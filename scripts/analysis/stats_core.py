"""Task 12 (part 1) - statistics on the Run 1 / Run 2 assemblies themselves.

1. GNNome vs Raven error rates: paired 100-kb windows (same genome positions), Wilcoxon signed-rank test and a
   genome block bootstrap (10,000 resamples of windows) for the difference in errors per 100 kb.
2. GNNome edge-score AUROC (true vs spurious): stratified bootstrap 95 % CI.
3. Junction cut-error model: Pearson r with Fisher-z 95 % CI; binomial test for indels >5 bp within 50 bp of a junction.
4. Repeat / homopolymer / low-complexity enrichment: all permutation p-values from task 9, Benjamini-Hochberg FDR.
5. Read selection: Mann-Whitney U and bootstrap CI for the difference in median per-read error (GNNome walk vs Raven layout).
"""
import json, csv, glob, math, pickle
import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score
A = "/home/abhiraj_arya/research/analysis"; G1 = "/home/abhiraj_arya/research/GNNome"; R2 = "/home/abhiraj_arya/research/run2_corrected"
N = 4583637; rng = np.random.default_rng(12); B = 10000
out = {}
tsv = lambda p: list(csv.DictReader(open(p), delimiter="\t"))

# 1. paired windows
def windows(path, cats, w=100000):
    pos = [int(e["ref"]) for e in tsv(path) if e["cat"] in cats]
    return np.histogram(pos, bins=np.arange(0, N + w, w))[0].astype(float)
comp = {}
for name, r, g, cats in (("run1_unpolished", "R1_Raven_p0", "R1_GNNome", ("mismatch", "indel_le5", "indel_gt5")),
                         ("run1_unpolished_indel_gt5", "R1_Raven_p0", "R1_GNNome", ("indel_gt5",)),
                         ("run2_unpolished", "R2_Raven_p0", "R2_GNNome", ("mismatch", "indel_le5", "indel_gt5")),
                         ("run2_unpolished_indel_gt5", "R2_Raven_p0", "R2_GNNome", ("indel_gt5",)),
                         ("run2_racon2", "R2_Raven_Racon2", "R2_GNNome_Racon2", ("mismatch", "indel_le5", "indel_gt5")),
                         ("run2_racon2_vs_raven_internal_p2", "R2_Raven_internal_p2", "R2_GNNome_Racon2", ("mismatch", "indel_le5", "indel_gt5"))):
    a = windows(f"{A}/09_repeat_errors/{r}.events.tsv", cats); b = windows(f"{A}/09_repeat_errors/{g}.events.tsv", cats)
    d = b - a; n = len(d)
    boot = np.array([d[rng.integers(0, n, n)].mean() for _ in range(B)])
    nz = d[d != 0]
    p = stats.wilcoxon(b, a).pvalue if len(nz) else 1.0
    comp[name] = dict(windows=n, raven_total=int(a.sum()), gnnome_total=int(b.sum()),
                      diff_per_100kb_mean=float(d.mean()), ci95=[float(np.quantile(boot, .025)), float(np.quantile(boot, .975))],
                      windows_gnnome_lower=int((d < 0).sum()), windows_gnnome_higher=int((d > 0).sum()), windows_equal=int((d == 0).sum()),
                      wilcoxon_p=float(p), ratio_gnnome_over_raven=float(b.sum() / a.sum()) if a.sum() else None)
out["gnnome_minus_raven_errors_per_100kb_window"] = comp

# 2. AUROC bootstrap
au = {}
for run in ("run1", "run2"):
    E = tsv(f"{A}/05_edge_predictions/{run}.edges.tsv")
    pos = np.array([float(e["logit"]) for e in E if e["reason"] == "true"])
    neg = np.array([float(e["logit"]) for e in E if e["reason"] in ("no_genomic_overlap", "wrong_direction", "strand_conflict")])
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]; s = np.r_[pos, neg]
    bs = []
    for _ in range(2000):
        pp = pos[rng.integers(0, len(pos), len(pos))]; nn = neg[rng.integers(0, len(neg), len(neg))]
        bs.append(roc_auc_score(np.r_[np.ones(len(pp)), np.zeros(len(nn))], np.r_[pp, nn]))
    au[run] = dict(auroc=float(roc_auc_score(y, s)), ci95=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))],
                   n_true=len(pos), n_spurious=len(neg),
                   note="only %d spurious edges: the CI is wide and the estimate is dominated by these few edges" % len(neg))
out["gnnome_auroc_true_vs_spurious"] = au

# 3. junction model
J = tsv(f"{A}/09_repeat_errors/R1_GNNome.junctions.tsv")
x = np.array([int(j["discrepancy"]) for j in J], float); y = np.array([int(j["observed_net_indel_within_50bp"]) for j in J], float)
r, p = stats.pearsonr(x, y); z = math.atanh(r); se = 1 / math.sqrt(len(x) - 3)
out["run1_junction_cut_error_model"] = dict(n=len(x), pearson_r=float(r), ci95=[math.tanh(z - 1.96 * se), math.tanh(z + 1.96 * se)], p=float(p))
jb = {}
for lab in ("R1_GNNome", "R2_GNNome"):
    s = json.load(open(f"{A}/09_repeat_errors/{lab}.summary.json"))["junction_distance"]["indel_gt5"]
    exp_frac = s["expected_within_50bp_if_random"] / s["n"]
    jb[lab] = dict(indel_gt5=s["n"], within_50bp=s["within_50bp"], expected=s["expected_within_50bp_if_random"],
                   binomial_p=float(stats.binom_test(s["within_50bp"], s["n"], exp_frac, alternative="greater")))
out["large_indels_at_gnnome_junctions"] = jb

# 4. enrichment with BH
tests = []
for f in sorted(glob.glob(f"{A}/09_repeat_errors/*.summary.json")):
    s = json.load(open(f))
    if "enrichment" not in s: continue
    for cat, tr in s["enrichment"].items():
        for track, v in tr.items():
            if v: tests.append((s["label"], cat, track, v["in_track"], v["n"], v["fold"], v["p_perm_ge"]))
p = np.array([t[-1] for t in tests]); m = len(p); order = np.argsort(p)
q = np.empty(m); prev = 1.0
for rank, i in reversed(list(enumerate(order, 1))):
    prev = min(prev, p[i] * m / rank); q[i] = prev
sig = [dict(assembly=t[0], events=t[1], track=t[2], in_track=t[3], n=t[4], fold=t[5], p=t[6], q=float(qq))
       for t, qq in zip(tests, q) if qq < 0.05]
out["enrichment_fdr"] = dict(tests=m, significant_q_lt_0p05=sig,
                             min_q=float(q.min()), note="one-sided circular-shift permutation tests (enrichment), BH across all tests")

# 5. read selection
rs = {}
for run, gfa, asm, wd in (("run1", f"{G1}/data/graph_1.gfa", f"{G1}/data/CP017100_assembly.fasta", f"{G1}/data/raven"),
                          ("run2", f"{R2}/02_raven_p0/graph_1.gfa", f"{R2}/02_raven_p0/raven_p0_assembly.fasta", f"{R2}/04_gnnome/raven")):
    T = {x["read"]: 1 - float(x["identity"]) for x in tsv(f"{A}/03_reads/{run}.reads.tsv")}
    rc = lambda s: s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    seg = {}; n2r = {}; idx = 0
    for line in open(gfa):
        f = line.rstrip("\n").split("\t")
        if f[0] == "S": seg[idx] = f[2]; seg[idx + 1] = rc(f[2]); n2r[idx] = n2r[idx + 1] = f[1]; idx += 2
    ctg = "".join(l.strip() for l in open(asm) if not l.startswith(">"))
    rav = {n2r[n] for n, s in seg.items() if ctg.find(s[:1000]) >= 0 and ctg.find(s[:1000], ctg.find(s[:1000]) + 1) < 0}
    gnn = {n2r[n] for n in pickle.load(open(f"{wd}/decode/0_walks.pkl", "rb"))[0]}
    g = np.array([T[r] for r in gnn]); rr = np.array([T[r] for r in rav])
    bd = [np.median(g[rng.integers(0, len(g), len(g))]) - np.median(rr[rng.integers(0, len(rr), len(rr))]) for _ in range(5000)]
    rs[run] = dict(gnnome_reads=len(g), raven_reads=len(rr), median_err_gnnome_pct=float(np.median(g) * 100), median_err_raven_pct=float(np.median(rr) * 100),
                   diff_median_pct=float((np.median(g) - np.median(rr)) * 100), ci95_pct=[float(np.quantile(bd, .025) * 100), float(np.quantile(bd, .975) * 100)],
                   mannwhitney_p=float(stats.mannwhitneyu(g, rr).pvalue))
out["read_selection"] = rs
json.dump(out, open(f"{A}/12_statistics/core_statistics.json", "w"), indent=1)
print(json.dumps(out, indent=1))
