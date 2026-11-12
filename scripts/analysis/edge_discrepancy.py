"""Cut-point error of every TRUE edge in a Raven graph.

For edge u->v the correct prefix (number of u bases before v starts) is computed from the base-level PBSIM
alignments (MAF, start-corrected). Raven's prefix = len(u_seg) - GFA overlap = Raven's edge length
(lhs_begin - rhs_begin of a minimizer-chain overlap). error = raven_prefix - correct_prefix.
Also marks, for each source node, its shortest-prefix true out-edge (what a transitive-reduction/"nearest
neighbour" layout such as Raven's own unitig uses) and the edges used by the GNNome walk.

usage: edge_discrepancy.py <gfa> <edges.tsv (graph_truth)> <maf.gz> <reads.fastq[.gz]> <gnnome_run_dir> <out_prefix>
"""
import sys, gzip, pickle, json, statistics as stt
import numpy as np
gfa, edges_tsv, maf, fq, run, out = sys.argv[1:7]
rc = lambda s: s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
op = lambda p: gzip.open(p, "rt") if p.endswith(".gz") else open(p)
GSEQ = "".join(l.strip() for l in open("/home/abhiraj_arya/research/GNNome/data/CP017100.fasta") if not l.startswith(">")).upper()

seg, n2r = {}, {}; idx = 0
for line in open(gfa):
    f = line.rstrip("\n").split("\t")
    if f[0] == "S": seg[f[1]] = f[2]; n2r[idx] = n2r[idx + 1] = f[1]; idx += 2
fqseq = {}
with op(fq) as f:
    for h in f:
        s = f.readline().strip(); f.readline(); f.readline(); n = h[1:].split()[0]
        if n in seg: fqseq[n] = s
rmap = {}; ref = None
for l in op(maf):
    if not l.startswith("s "): continue
    p = l.split()
    if p[1] == "ref": ref = p; continue
    if p[1] in seg:
        rpos = int(ref[2]); rtxt = ref[6].replace("-", "").upper()
        if GSEQ[rpos:rpos + len(rtxt)] != rtxt: rpos = GSEQ.find(rtxt, max(0, rpos - 5000), rpos + 5000 + len(rtxt))
        rt = np.frombuffer(ref[6].encode(), dtype=np.uint8); qt = np.frombuffer(p[6].encode(), dtype=np.uint8)
        refpos = rpos + np.cumsum(rt != 45) - (rt != 45)      # ref coordinate of each column (45 = '-')
        m = refpos[qt != 45].astype(np.int64)
        s = seg[p[1]]; t5 = fqseq[p[1]].find(s)
        mm = m if p[4] == "+" else m[::-1]                      # index by FASTQ offset
        rmap[p[1]] = mm[t5:t5 + len(s)]
    ref = None
def pos(n):
    m = rmap[n2r[n]]; return m if n % 2 == 0 else m[::-1]
walk = pickle.load(open(f"{run}/decode/0_walks.pkl", "rb"))[0]
wset = set(zip(walk[:-1], walk[1:]))
rows = []
for i, l in enumerate(open(edges_tsv)):
    if i == 0: continue
    f = l.rstrip("\n").split("\t")
    if f[15] != "true": continue
    u, v, ol = int(f[1]), int(f[2]), int(f[11])
    pu, pv = pos(u), pos(v); target = pv[0]
    if pu[-1] >= pu[0]: c = int(np.searchsorted(pu, target, "left"))
    else: c = int(np.searchsorted(-pu, -target, "left"))
    rp = len(seg[n2r[u]]) - ol
    rows.append((u, v, ol, rp, c, rp - c, float(f[17]), (u, v) in wset))
best = {}
for r in rows:
    if r[0] not in best or r[3] < best[r[0]][3]: best[r[0]] = r
nn = set((r[0], r[1]) for r in best.values())
def summ(sel):
    e = np.array([abs(r[5]) for r in sel]); p = np.array([r[3] for r in sel])
    return dict(n=len(sel), abs_err_median=float(np.median(e)), abs_err_mean=float(e.mean()), frac_abs_gt5=float((e > 5).mean()),
                prefix_median=float(np.median(p)))
bins = [0, 1000, 2000, 4000, 6000, 8000, 10000, 15000, 40000]
by_prefix = []
for a, b in zip(bins, bins[1:]):
    sel = [r for r in rows if a <= r[3] < b]
    if sel: by_prefix.append(dict(prefix_bin=f"{a}-{b}", **summ(sel)))
e_all = np.array([r[5] for r in rows]); p_all = np.array([r[3] for r in rows])
S = dict(true_edges=len(rows), all_true_edges=summ(rows), nearest_neighbour_edges=summ([r for r in rows if (r[0], r[1]) in nn]),
         gnnome_walk_edges=summ([r for r in rows if r[7]]), by_prefix_length=by_prefix,
         spearman_abs_err_vs_prefix=None)
from scipy.stats import spearmanr
rho, pv = spearmanr(np.abs(e_all), p_all); S["spearman_abs_err_vs_prefix"] = dict(rho=float(rho), p=float(pv))
json.dump(S, open(out + ".summary.json", "w"), indent=1)
with open(out + ".tsv", "w") as w:
    w.write("src\tdst\toverlap\traven_prefix\tcorrect_prefix\terror\tgnnome_prob\tin_walk\n")
    for r in rows: w.write("\t".join(str(x) for x in r) + "\n")
print(json.dumps(S, indent=1))
