"""Per-read ground truth + error profile from a PBSIM3 MAF, plus FASTQ quality.
usage: read_truth.py <maf.gz> <reads.fastq[.gz]> <genome_len> <out_prefix>
writes <out>.reads.tsv, <out>.coverage_1kb.tsv, <out>.summary.json"""
import sys, gzip, json, math
maf, fq, G, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
GENOME = sys.argv[5] if len(sys.argv) > 5 else "/home/abhiraj_arya/research/GNNome/data/CP017100.fasta"
GSEQ = "".join(l.strip() for l in open(GENOME) if not l.startswith(">")).upper()
fixed = []   # PBSIM3 builds can write shifted ref starts for '-' reads: verify every block against the genome
op = lambda p: gzip.open(p, "rt") if p.endswith(".gz") else open(p)
comp = str.maketrans("ACGTacgt", "TGCAtgca")
rows = {}; ref = None
for l in op(maf):
    if not l.startswith("s "): continue
    f = l.split()
    if f[1] == "ref": ref = f; continue
    rstart, rspan = int(ref[2]), int(ref[3]); rt, qt = ref[6], f[6]
    rtxt = rt.replace("-", "").upper()
    if GSEQ[rstart:rstart + len(rtxt)] != rtxt:
        k = GSEQ.find(rtxt, max(0, rstart - 5000), rstart + 5000 + len(rtxt))
        if k < 0: raise SystemExit(f"cannot place MAF block of {f[1]}")
        fixed.append(k - rstart); rstart = k
    rspan = len(rtxt)
    m = x = ins = dele = 0; hp_err = 0
    for i, (a, b) in enumerate(zip(rt, qt)):
        if a == b: m += 1
        elif a == "-": ins += 1
        elif b == "-": dele += 1
        else: x += 1
    cols = len(rt)
    rows[f[1]] = dict(read=f[1], strand=f[4], ref_start=rstart, ref_end=rstart + rspan, read_len=int(f[5]),
                      match=m, mismatch=x, ins=ins, dele=dele, identity=m / cols,
                      _qt=qt.replace("-", ""))
    ref = None
# FASTQ: mean phred, and orientation check against MAF
checked = okc = 0
with op(fq) as f:
    for h in f:
        s = f.readline().strip(); f.readline(); q = f.readline().strip()
        r = rows[h[1:].split()[0]]
        r["mean_q"] = sum(ord(c) - 33 for c in q) / len(q)
        if checked < 500:
            exp = r["_qt"] if r["strand"] == "+" else r["_qt"].translate(comp)[::-1]
            okc += (exp.upper() == s.upper()); checked += 1
for r in rows.values(): del r["_qt"]
cols = ["read","strand","ref_start","ref_end","read_len","match","mismatch","ins","dele","identity","mean_q"]
with open(out + ".reads.tsv", "w") as w:
    w.write("\t".join(cols) + "\n")
    for r in rows.values(): w.write("\t".join(str(round(r[c], 5)) if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
# true coverage in 1 kb bins (difference array)
d = [0] * (G + 1)
for r in rows.values(): d[r["ref_start"]] += 1; d[r["ref_end"]] -= 1
cov, c = [], 0
for i in range(G): c += d[i]; cov.append(c)
with open(out + ".coverage_1kb.tsv", "w") as w:
    w.write("bin_start\tmean_depth\tmin_depth\n")
    for b in range(0, G, 1000):
        seg = cov[b:b + 1000]; w.write(f"{b}\t{sum(seg)/len(seg):.2f}\t{min(seg)}\n")
L = sorted(r["read_len"] for r in rows.values()); tot = sum(L)
acc = 0; n50 = 0
for x in sorted(L, reverse=True):
    acc += x
    if acc >= tot / 2: n50 = x; break
ids = sorted(r["identity"] for r in rows.values())
TM = sum(r["match"] for r in rows.values()); TX = sum(r["mismatch"] for r in rows.values())
TI = sum(r["ins"] for r in rows.values()); TD = sum(r["dele"] for r in rows.values()); TC = TM + TX + TI + TD
pct = lambda a, p: a[min(len(a) - 1, int(p * len(a)))]
mean = lambda a: sum(a) / len(a)
sd = lambda a: math.sqrt(sum((v - mean(a)) ** 2 for v in a) / (len(a) - 1))
S = dict(reads=len(L), bases=tot, coverage=tot / G, len_mean=mean(L), len_sd=sd(L), len_min=L[0], len_max=L[-1], len_n50=n50,
         len_q25=pct(L, .25), len_median=pct(L, .5), len_q75=pct(L, .75),
         identity_aggregate=TM / TC, identity_read_mean=mean(ids), identity_read_sd=sd(ids),
         identity_p05=pct(ids, .05), identity_median=pct(ids, .5), identity_p95=pct(ids, .95),
         mismatch_rate=TX / TC, insertion_rate=TI / TC, deletion_rate=TD / TC,
         error_share_mismatch=TX / (TX + TI + TD), error_share_ins=TI / (TX + TI + TD), error_share_del=TD / (TX + TI + TD),
         mean_q_mean=mean([r["mean_q"] for r in rows.values()]),
         strand_plus=sum(r["strand"] == "+" for r in rows.values()),
         depth_mean=mean(cov), depth_min=min(cov), depth_max=max(cov),
         bp_depth_lt5=sum(c < 5 for c in cov), bp_depth_lt10=sum(c < 10 for c in cov), bp_depth0=sum(c == 0 for c in cov),
         orientation_check=f"{okc}/{checked} FASTQ reads equal MAF read sequence (revcomp for '-')",
         maf_ref_start_corrected=len(fixed), maf_ref_start_shift_max=max(fixed) if fixed else 0,
         maf_ref_start_shift_median=sorted(fixed)[len(fixed)//2] if fixed else 0)
json.dump(S, open(out + ".summary.json", "w"), indent=1)
print(json.dumps(S, indent=1))
