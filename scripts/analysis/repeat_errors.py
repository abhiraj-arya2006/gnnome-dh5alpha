"""Task 9 - repeat-associated errors, investigated down to graph and read level.

For one assembly:
  1. Error events from QUAST used_snps (mismatch / indel <=5 bp / indel >5 bp, ref + contig coordinates)
     and recalculated breakpoints (task 8).
  2. Genomic context of every event: repeat family/class, annotated element, homopolymer (>=7), low complexity,
     true read depth, mean identity of the reads covering the site.
  3. Enrichment vs. the genome background with a circular-shift permutation test (1000 shifts).
  4. Read-level origin (unpolished assemblies only, where every base is copied from one read):
       inherited  - one read contains the exact contig window around the event -> the read carried the error
       junction   - no single read contains the window -> created where two reads were joined
  5. Graph level (GNNome only): the contig is rebuilt from the walk (reads[src][:prefix] ... + last read) and checked
     byte-for-byte against the assembly, giving exact junction coordinates, the two reads, the GFA overlap, the edge's
     ground-truth class and GNNome probability, and how many alternative/spurious out-edges the source node had.
  6. A dossier row for every large indel (>5 bp) and every breakpoint.

usage: repeat_errors.py <label> <assembly.fasta> <used_snps> <breakpoints.tsv|-> <truth reads.tsv> <reads.fastq[.gz]>
                        <gfa|-> <gnnome_run_dir|-> <edges.tsv|-> <unpolished 0|1> <outdir>
"""
import sys, gzip, json, bisect, random, pickle, os
from collections import Counter, defaultdict
import numpy as np

(label, asm_path, snps_path, bp_path, truth_path, reads_path, gfa_path, run_dir, edges_path,
 unpolished, outdir) = sys.argv[1:12]
unpolished = unpolished == "1"
R7 = "/home/abhiraj_arya/research/analysis/07_repeats"
G = 4583637
random.seed(1)
rc = lambda s: s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
GSEQ = "".join(l.strip() for l in open("/home/abhiraj_arya/research/GNNome/data/CP017100.fasta") if not l.startswith(">")).upper()

# ---------------- tracks ----------------
def bed(path, cols):
    out = []
    for l in open(path):
        f = l.rstrip("\n").split("\t"); out.append(tuple([int(f[1]), int(f[2])] + [f[c] for c in cols]))
    return out
repeats = bed(f"{R7}/repeats.bed", [3, 4, 5, 6])          # locus, family, class, copies
elements = bed(f"{R7}/elements.bed", [3, 4])              # class, name
hps = bed(f"{R7}/homopolymers.bed", [3])
dust = bed(f"{R7}/lowcomplexity.bed", [3])
def mask(ivs, pad=0):
    m = np.zeros(G, dtype=bool)
    for iv in ivs: m[max(0, iv[0] - pad):min(G, iv[1] + pad)] = True
    return m
M = {"repeat": mask(repeats), "homopolymer": mask(hps, 1), "lowcomplexity": mask(dust)}
for cls in ("rRNA", "IS_transposase", "Rhs", "phage_related", "unannotated"):
    M["repeat_" + cls] = mask([r for r in repeats if r[4] == cls])
rep_idx = np.full(G, -1, dtype=np.int32)
for i, r in enumerate(repeats): rep_idx[r[0]:r[1]] = i
el_idx = np.full(G, -1, dtype=np.int32)
for i, e in enumerate(elements): el_idx[e[0]:e[1]] = i

# ---------------- truth: reads, depth, local read identity ----------------
T = {}
for i, l in enumerate(open(truth_path)):
    if i == 0: continue
    f = l.rstrip("\n").split("\t")
    T[f[0]] = (f[1], int(f[2]), int(f[3]), float(f[9]))
diff = np.zeros(G + 1); idsum = np.zeros(G + 1)
for s, a, b, idy in T.values():
    diff[a] += 1; diff[b] -= 1; idsum[a] += idy; idsum[b] -= idy
depth = np.cumsum(diff)[:G]; mean_id = np.cumsum(idsum)[:G] / np.maximum(depth, 1)
starts = sorted((a, r) for r, (s, a, b, idy) in T.items()); st_pos = [x[0] for x in starts]
maxlen = max(b - a for s, a, b, idy in T.values())
def reads_at(p):
    j = bisect.bisect_right(st_pos, p); out = []
    for k in range(j - 1, -1, -1):
        a, r = starts[k]
        if a < p - maxlen: break
        if T[r][2] > p: out.append(r)
    return out

# ---------------- assembly + events ----------------
asm = "".join(l.strip() for l in open(asm_path) if not l.startswith(">")).upper()
events = []
for l in open(snps_path):
    f = l.rstrip("\n").split("\t")
    rp, rb, cb, cp = int(f[2]) - 1, f[3], f[4], int(f[5]) - 1
    if rb == "." : typ, ln = "ins", len(cb)
    elif cb == ".": typ, ln = "del", len(rb)
    else: typ, ln = "mismatch", len(rb)
    cat = "mismatch" if typ == "mismatch" else ("indel_le5" if ln <= 5 else "indel_gt5")
    events.append(dict(kind="small", cat=cat, typ=typ, len=ln, ref=rp, ctg=cp))
bps = []
if bp_path != "-" and os.path.exists(bp_path):
    seen = []
    for i, l in enumerate(open(bp_path)):
        if i == 0: continue
        f = l.rstrip("\n").split("\t")
        if f[2] == "origin": continue
        pos = int(f[3]); src = f[1]
        hit = next((s for s in seen if abs(s["ref"] - pos) < 2000), None)
        if hit: hit["sources"].add(src); continue
        d = dict(kind="breakpoint", cat="breakpoint", typ=f[2], len=abs(int(f[9])) if f[9] else 0, ref=pos,
                 ctg=int(f[5]), sources={src}); seen.append(d)
    bps = seen
allev = events + bps

# ---------------- GNNome junctions (exact) ----------------
junctions = []
if run_dir != "-":
    walks = pickle.load(open(f"{run_dir}/decode/0_walks.pkl", "rb"))
    seg, ol, r2n = {}, {}, {}
    idx = 0; n2r = {}
    for line in open(gfa_path):
        f = line.rstrip("\n").split("\t")
        if f[0] == "S":
            seg[idx] = f[2]; seg[idx + 1] = rc(f[2]); n2r[idx] = n2r[idx + 1] = f[1]; r2n[f[1]] = (idx, idx + 1); idx += 2
        elif f[0] == "L":
            a, o1, b, o2, cig = f[1:6]; o = int(cig[:-1])
            A, B = r2n[a], r2n[b]
            sr = A[0] if o1 == "+" else A[1]; dr = B[0] if o2 == "+" else B[1]
            sv = B[1] if o2 == "+" else B[0]; dv = A[1] if o1 == "+" else A[0]
            ol[(sr, dr)] = o; ol[(sv, dv)] = o
    E = {}
    if edges_path != "-":
        for i, l in enumerate(open(edges_path)):
            if i == 0: continue
            f = l.rstrip("\n").split("\t"); E[(int(f[1]), int(f[2]))] = (f[15], float(f[17]))
    outdeg = Counter(); spur_out = Counter()
    for (u, v), (why, p) in E.items():
        outdeg[u] += 1
        if why in ("no_genomic_overlap", "wrong_direction", "strand_conflict"): spur_out[u] += 1
    w = walks[0]; pieces = []; pos = 0
    for u, v in zip(w[:-1], w[1:]):
        pre = max(0, len(seg[u]) - ol[(u, v)])  # inference.py clamps negative prefixes to 0
        pieces.append(seg[u][:pre]); pos += pre
        cls, prob = E.get((u, v), ("?", float("nan")))
        junctions.append(dict(ctg=pos, src=n2r[u], dst=n2r[v], src_node=u, dst_node=v, overlap=ol[(u, v)], prefix=pre,
                              edge_class=cls, prob=prob, src_outdeg=outdeg[u], src_spurious_out=spur_out[u]))
    pieces.append(seg[w[-1]])
    rebuilt = "".join(pieces).upper()
    walk_rebuild_ok = rebuilt == asm
    jpos = [j["ctg"] for j in junctions]
else:
    walk_rebuild_ok = None; jpos = []
def nearest_junction(cp):
    if not jpos: return None, None
    k = bisect.bisect_left(jpos, cp); best = None
    for kk in (k - 1, k):
        if 0 <= kk < len(jpos) and (best is None or abs(jpos[kk] - cp) < abs(jpos[best] - cp)): best = kk
    return jpos[best] - cp, junctions[best]

# ---------------- read-level origin ----------------
reads = {}
if unpolished:
    need = set()
    sample_small = [e for e in events if e["cat"] != "indel_gt5"]
    random.shuffle(sample_small); sample_small = sample_small[:3000]
    test = sample_small + [e for e in events if e["cat"] == "indel_gt5"] + bps
    for e in test:
        for r in reads_at(min(G - 1, max(0, e["ref"]))): need.add(r)
    op = gzip.open if reads_path.endswith(".gz") else open
    with op(reads_path, "rt") as f:
        for h in f:
            s = f.readline().strip(); f.readline(); f.readline()
            n = h[1:].split()[0]
            if n in need: reads[n] = s.upper()
    FL = 30
    for e in test:
        a, b = e["ctg"] - FL, e["ctg"] + max(e["len"], 1) + FL
        if e["typ"] == "del": b = e["ctg"] + FL
        if a < 0 or b > len(asm): e["origin"] = "edge_of_contig"; continue
        win = asm[a:b]
        lo, hi = max(0, e["ref"] - 20000), min(G, e["ref"] + 20000)
        if win in GSEQ[lo:hi] or rc(win) in GSEQ[lo:hi]:
            # the joined sequence exists in the genome itself (e.g. collapse/expansion between repeat copies):
            # matching reads would not prove a read error, so this is reported separately
            e["origin"] = "genomic_window"; e["source_reads"] = []; continue
        srcs = [r for r in reads_at(e["ref"]) if r in reads and (win in reads[r] or win in rc(reads[r]))]
        e["origin"] = "inherited" if srcs else "junction"
        e["source_reads"] = srcs[:3]
        if srcs: e["source_identity"] = T[srcs[0]][3]

# ---------------- annotate ----------------
def ctx(e):
    p = min(G - 1, max(0, e["ref"]))
    ri, ei = rep_idx[p], el_idx[p]
    e["repeat_family"] = repeats[ri][3] if ri >= 0 else ""; e["repeat_class"] = repeats[ri][4] if ri >= 0 else ""
    e["repeat_copies"] = repeats[ri][5] if ri >= 0 else ""
    e["element"] = f"{elements[ei][2]}:{elements[ei][3]}" if ei >= 0 else ""
    e["homopolymer"] = bool(M["homopolymer"][p]); e["lowcomplexity"] = bool(M["lowcomplexity"][p])
    e["depth"] = int(depth[p]); e["reads_mean_identity"] = round(float(mean_id[p]), 5)
    if jpos:
        d, j = nearest_junction(e["ctg"]); e["junction_dist"] = d; e["junction"] = j
for e in allev: ctx(e)

# ---------------- enrichment (circular-shift permutation) ----------------
def enrich(evs, track, nperm=1000):
    if not evs: return None
    pos = np.array([min(G - 1, max(0, e["ref"])) for e in evs])
    obs = int(M[track][pos].sum()); exp_frac = float(M[track].mean())
    rng = np.random.default_rng(7); ge = 0
    for _ in range(nperm):
        if int(M[track][(pos + rng.integers(0, G)) % G].sum()) >= obs: ge += 1
    return dict(n=len(evs), in_track=obs, frac=obs / len(evs), genome_frac=exp_frac,
                fold=(obs / len(evs)) / exp_frac if exp_frac else None, p_perm_ge=(ge + 1) / (nperm + 1))
groups = defaultdict(list)
for e in allev: groups[e["cat"]].append(e)
enr = {cat: {tr: enrich(evs, tr) for tr in M} for cat, evs in groups.items()}

# ---------------- outputs ----------------
os.makedirs(outdir, exist_ok=True)
summary = dict(label=label, assembly_len=len(asm), events=dict(Counter(e["cat"] for e in allev)),
               gnnome_walk_rebuild_matches_assembly=walk_rebuild_ok, gnnome_junctions=len(junctions), enrichment=enr)
if unpolished:
    tested = [e for e in allev if "origin" in e]
    summary["read_origin"] = {cat: dict(Counter(e["origin"] for e in tested if e["cat"] == cat)) for cat in groups}
    inh = [e["source_identity"] for e in tested if e.get("origin") == "inherited" and e["cat"] != "breakpoint"]
    summary["inherited_source_read_identity_mean"] = float(np.mean(inh)) if inh else None
    summary["all_reads_identity_mean"] = float(np.mean([v[3] for v in T.values()]))
if jpos:
    for cat, evs in groups.items():
        d = [abs(e["junction_dist"]) for e in evs]
        summary.setdefault("junction_distance", {})[cat] = dict(
            n=len(d), within_50bp=sum(x <= 50 for x in d), within_500bp=sum(x <= 500 for x in d),
            expected_within_50bp_if_random=len(d) * min(1.0, 101 * len(jpos) / len(asm)))
    if unpolished:
        jt = [e for e in allev if e.get("origin") == "junction"]
        summary["junction_origin_events_within_50bp_of_walk_junction"] = f"{sum(abs(e['junction_dist']) <= 50 for e in jt)}/{len(jt)}"
json.dump(summary, open(f"{outdir}/{label}.summary.json", "w"), indent=1, default=str)
cols = ["kind", "cat", "typ", "len", "ref", "ctg", "repeat_family", "repeat_class", "repeat_copies", "element", "homopolymer",
        "lowcomplexity", "depth", "reads_mean_identity", "origin", "source_reads", "source_identity", "junction_dist"]
jcols = ["src", "dst", "overlap", "prefix", "edge_class", "prob", "src_outdeg", "src_spurious_out"]
with open(f"{outdir}/{label}.dossier.tsv", "w") as w:
    w.write("\t".join(cols + (["junction_" + c for c in jcols] if jpos else []) + ["sources"]) + "\n")
    for e in sorted([e for e in allev if e["cat"] in ("indel_gt5", "breakpoint")], key=lambda e: e["ref"]):
        row = [str(e.get(c, "")) for c in cols]
        if jpos: row += [str(e["junction"][c]) for c in jcols]
        row.append(",".join(sorted(e.get("sources", []))))
        w.write("\t".join(row) + "\n")
with open(f"{outdir}/{label}.events.tsv", "w") as w:
    w.write("cat\ttyp\tlen\tref\tctg\trepeat_class\thomopolymer\tlowcomplexity\tdepth\torigin\tjunction_dist\n")
    for e in allev:
        w.write("\t".join(str(e.get(c, "")) for c in ["cat", "typ", "len", "ref", "ctg", "repeat_class", "homopolymer",
                                                      "lowcomplexity", "depth", "origin", "junction_dist"]) + "\n")
print(json.dumps({k: v for k, v in summary.items() if k != "enrichment"}, indent=1, default=str))
for cat in enr:
    for tr in ("repeat", "homopolymer", "lowcomplexity", "repeat_rRNA", "repeat_IS_transposase"):
        r = enr[cat][tr]
        if r: print(f"  {cat:10s} {tr:22s} {r['in_track']:>6}/{r['n']:<7} {100*r['frac']:6.2f}% vs genome {100*r['genome_frac']:5.2f}%  fold {r['fold']:.2f}  p={r['p_perm_ge']:.4f}")
