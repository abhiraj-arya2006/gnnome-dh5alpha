"""Quantify GNNome junction (cut-point) errors from the reads' true alignments.

GNNome builds contig = seg(n0)[:prefix0] + seg(n1)[:prefix1] + ... + seg(nk), prefix = max(0, len(seg) - GFA overlap).
The join u->v is exact when u_seg[prefix] sits at the same reference base as v_seg[0]. Using the PBSIM MAF
alignment of every read (base-level read->reference map) and the read's Raven trim offset, the correct prefix is
computed; discrepancy = used_prefix - correct_prefix (>0 duplicated bases = insertion in the contig, <0 = deletion).
The predicted discrepancies are compared with the indels QUAST observed within 50 bp of each junction.

usage: junction_check.py <gfa> <gnnome_run_dir> <maf.gz> <reads.fastq[.gz]> <used_snps> <out_prefix>
"""
import sys, gzip, pickle, json, bisect
from collections import Counter
gfa, run, maf, fq, snps, out = sys.argv[1:7]
rc = lambda s: s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
op = lambda p: gzip.open(p, "rt") if p.endswith(".gz") else open(p)
GSEQ = "".join(l.strip() for l in open("/home/abhiraj_arya/research/GNNome/data/CP017100.fasta") if not l.startswith(">")).upper()

seg, n2r, r2n, ol = {}, {}, {}, {}
idx = 0
for line in open(gfa):
    f = line.rstrip("\n").split("\t")
    if f[0] == "S":
        seg[idx] = f[2]; n2r[idx] = n2r[idx + 1] = f[1]; r2n[f[1]] = (idx, idx + 1); idx += 2
    elif f[0] == "L":
        a, o1, b, o2, cig = f[1:6]; o = int(cig[:-1]); A, B = r2n[a], r2n[b]
        sr = A[0] if o1 == "+" else A[1]; dr = B[0] if o2 == "+" else B[1]
        sv = B[1] if o2 == "+" else B[0]; dv = A[1] if o1 == "+" else A[0]
        ol[(sr, dr)] = o; ol[(sv, dv)] = o
walk = pickle.load(open(f"{run}/decode/0_walks.pkl", "rb"))[0]
need = {n2r[n] for n in walk}

# FASTQ reads (for trim offsets)
fqseq = {}
with op(fq) as f:
    for h in f:
        s = f.readline().strip(); f.readline(); f.readline()
        n = h[1:].split()[0]
        if n in need: fqseq[n] = s
# MAF: read (FASTQ orientation offset) -> reference position
rmap, strand = {}, {}
ref = None
for l in op(maf):
    if not l.startswith("s "): continue
    p = l.split()
    if p[1] == "ref": ref = p; continue
    name = p[1]
    if name in need:
        rpos = int(ref[2]); rtxt = ref[6].replace("-", "").upper()
        if GSEQ[rpos:rpos + len(rtxt)] != rtxt:          # PBSIM3 '-' strand start shift (see read_truth.py)
            rpos = GSEQ.find(rtxt, max(0, rpos - 5000), rpos + 5000 + len(rtxt)); assert rpos >= 0
        m = []
        for a, b in zip(ref[6], p[6]):
            if b != "-": m.append(rpos)
            if a != "-": rpos += 1
        L = len(m)
        strand[name] = p[4]
        rmap[name] = m if p[4] == "+" else m[::-1]     # index by FASTQ-orientation offset
    ref = None

def node_refpos(n):
    r = n2r[n]; s = seg[n - (n % 2)]; t5 = fqseq[r].find(s)
    assert t5 >= 0
    m = rmap[r][t5:t5 + len(s)]
    return m if n % 2 == 0 else m[::-1]                 # odd node = reverse complement

rows = []
for u, v in zip(walk[:-1], walk[1:]):
    used = max(0, len(seg[u - (u % 2)]) - ol[(u, v)])
    pu, pv = node_refpos(u), node_refpos(v)
    target = pv[0]
    inc = pu[-1] >= pu[0]
    if inc: i = bisect.bisect_left(pu, target)
    else:
        neg = [-x for x in pu]; i = bisect.bisect_left(neg, -target)
    rows.append(dict(src=n2r[u], dst=n2r[v], overlap=ol[(u, v)], used_prefix=used, correct_prefix=i, discrepancy=used - i))
# contig coordinates of junctions and observed indels nearby
pos = 0; jc = []
for r in rows: pos += r["used_prefix"]; jc.append(pos)
ev = []
for l in open(snps):
    f = l.rstrip("\n").split("\t")
    rb, cb, cp = f[3], f[4], int(f[5]) - 1
    if rb == ".": ev.append((cp, len(cb)))
    elif cb == ".": ev.append((cp, -len(rb)))
ev.sort(); evp = [e[0] for e in ev]
for r, c in zip(rows, jc):
    a, b = bisect.bisect_left(evp, c - 50), bisect.bisect_right(evp, c + 50)
    r["observed_net_indel_within_50bp"] = sum(e[1] for e in ev[a:b])
    r["contig_pos"] = c
d = [r["discrepancy"] for r in rows]
exact = [r for r in rows if r["discrepancy"] == r["observed_net_indel_within_50bp"]]
big = [r for r in rows if abs(r["discrepancy"]) > 5]
big_ok = [r for r in big if r["discrepancy"] == r["observed_net_indel_within_50bp"]]
S = dict(junctions=len(rows), discrepancy_zero=sum(x == 0 for x in d), abs_le5=sum(abs(x) <= 5 for x in d),
         abs_gt5=len(big), discrepancy_hist=dict(sorted(Counter(max(-30, min(30, x)) for x in d).items())),
         predicted_equals_observed_net_indel=f"{len(exact)}/{len(rows)}",
         predicted_gt5_equals_observed=f"{len(big_ok)}/{len(big)}",
         pearson_r=None)
import statistics as st
x = [r["discrepancy"] for r in rows]; y = [r["observed_net_indel_within_50bp"] for r in rows]
mx, my = st.mean(x), st.mean(y); sx = sum((a - mx) ** 2 for a in x) ** .5; sy = sum((b - my) ** 2 for b in y) ** .5
S["pearson_r"] = sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy) if sx and sy else None
json.dump(S, open(out + ".summary.json", "w"), indent=1)
with open(out + ".tsv", "w") as w:
    k = list(rows[0].keys()); w.write("\t".join(k) + "\n")
    for r in rows: w.write("\t".join(str(r[c]) for c in k) + "\n")
print(json.dumps(S, indent=1))
