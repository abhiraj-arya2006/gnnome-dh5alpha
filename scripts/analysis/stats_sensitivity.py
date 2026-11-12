"""Task 6 + 12 (part 2): sensitivity experiment summary, paired GNNome-vs-Raven statistics, figures F8/F9.

Pairing: within a replicate both assemblers use the same reads and the same Raven graph, so GNNome - Raven
differences are paired. Per condition: mean, SD and t-based 95 % CI of the paired difference (n = seeds).
Across all replicates: Wilcoxon signed-rank. Effect of read accuracy on the difference: OLS slope (accuracy series).
usage: stats_sensitivity.py <collected.tsv>
"""
import sys, csv, json, math
import numpy as np
from scipy import stats
sys.path.insert(0, "/home/abhiraj_arya/research/analysis/scripts")
import figures as F
import matplotlib.pyplot as plt
A = "/home/abhiraj_arya/research/analysis"
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
num = lambda x: float(x) if x not in ("", None) else float("nan")
for r in rows:
    for k in list(r):
        if k not in ("replicate", "cond"): r[k] = num(r[k])
    r["raven_racon2_err"] = r["raven_racon2_mm"] + r["raven_racon2_indel"]
    r["gnnome_racon2_err"] = r["gnnome_racon2_mm"] + r["gnnome_racon2_indel"]
    r["d_p0_identity_pp"] = 100 * (r["gnnome_p0_identity"] - r["raven_p0_identity"])
    r["d_racon2_err"] = r["gnnome_racon2_err"] - r["raven_racon2_err"]
    r["d_p0_localmis"] = r["gnnome_p0m80_localmis"] - r["raven_p0m80_localmis"]
    r["d_source_err_pp"] = 100 * (r["gnnome_source_read_error"] - r["raven_source_read_error"])
ORDER = ["acc085_d30_L9k", "acc090_d30_L9k", "acc095_d30_L9k", "acc099_d30_L9k", "acc085_d60_L9k", "acc085_d30_L18k", "acc099_d60_L18k", "run1reads_d30"]
NAME = {"acc085_d30_L9k": "acc 0.85 · 30× · 9 kb (= Run 2 settings)", "acc090_d30_L9k": "acc 0.90 · 30× · 9 kb", "acc095_d30_L9k": "acc 0.95 · 30× · 9 kb",
        "acc099_d30_L9k": "acc 0.99* · 30× · 9 kb", "acc085_d60_L9k": "acc 0.85 · 60× · 9 kb", "acc085_d30_L18k": "acc 0.85 · 30× · 18.5 kb",
        "acc099_d60_L18k": "acc 0.99* · 60× · 18.5 kb", "run1reads_d30": "Run 1 reads, 50 % subsample (~30×)"}
conds = [c for c in ORDER if any(r["cond"] == c for r in rows)]
def tci(x):
    x = np.array([v for v in x if not math.isnan(v)]); n = len(x)
    if n < 2: return (float(x.mean()) if n else float("nan"), float("nan"), float("nan"), n)
    m, s = x.mean(), x.std(ddof=1); h = stats.t.ppf(.975, n - 1) * s / math.sqrt(n)
    return float(m), float(m - h), float(m + h), n
summary = {}
metrics = ["all_reads_error_mean", "raven_source_read_error", "gnnome_source_read_error", "read_identity", "read_coverage", "read_len_n50", "graph_edges", "edges_per_node", "edges_spurious", "edges_containment",
           "gnn_auroc", "gnn_logit_true_median", "walk_spurious_or_backward", "raven_p0_identity", "gnnome_p0_identity",
           "raven_p0m80_localmis", "gnnome_p0m80_localmis", "raven_racon2_err", "gnnome_racon2_err", "raven_racon2_mis", "gnnome_racon2_mis",
           "raven_racon2_gf", "gnnome_racon2_gf", "raven_p0_contigs", "gnnome_p0_contigs", "d_p0_identity_pp", "d_racon2_err", "d_p0_localmis", "d_source_err_pp"]
for c in conds:
    R = [r for r in rows if r["cond"] == c]
    summary[c] = {m: dict(zip(("mean", "ci_lo", "ci_hi", "n"), tci([r[m] for r in R]))) for m in metrics}
    for d in ("d_p0_identity_pp", "d_racon2_err"):
        x = [r[d] for r in R]; summary[c][d]["gnnome_better_in"] = f"{sum((v > 0) if d == 'd_p0_identity_pp' else (v < 0) for v in x)}/{len(x)}"
tests = {}
for d, better in (("d_p0_identity_pp", 1), ("d_racon2_err", -1), ("d_source_err_pp", -1)):
    x = np.array([r[d] for r in rows if not math.isnan(r[d])])
    w = stats.wilcoxon(x) if np.any(x != 0) else None
    tests[d] = dict(n=len(x), mean=float(x.mean()), median=float(np.median(x)), gnnome_better=int(np.sum(np.sign(x) == better)),
                    wilcoxon_p=float(w.pvalue) if w else 1.0)
acc = [r for r in rows if r["cond"] in ORDER[:4]]
if len(acc) > 3:
    lr = stats.linregress([r["read_identity"] * 100 for r in acc], [r["d_p0_identity_pp"] for r in acc])
    tests["slope_identity_diff_vs_read_identity"] = dict(slope_pp_per_pp=lr.slope, p=lr.pvalue, r=lr.rvalue, n=len(acc))
json.dump(dict(conditions=summary, pooled_tests=tests, n_replicates=len(rows)), open(f"{A}/12_statistics/sensitivity_statistics.json", "w"), indent=1)
with open(f"{A}/06_sensitivity/condition_summary.tsv", "w") as w:
    w.write("condition\tn\t" + "\t".join(f"{m}_mean\t{m}_ci95" for m in metrics) + "\n")
    for c in conds:
        s = summary[c]; w.write(f"{c}\t{s['read_identity']['n']}\t" + "\t".join(
            f"{s[m]['mean']:.6g}\t[{s[m]['ci_lo']:.4g}, {s[m]['ci_hi']:.4g}]" for m in metrics) + "\n")
print(json.dumps(tests, indent=1))
for c in conds:
    s = summary[c]
    print(f"{c:18s} n={s['read_identity']['n']} readID={100*s['read_identity']['mean']:.2f}% spur={s['edges_spurious']['mean']:.1f} "
          f"logitT={s['gnn_logit_true_median']['mean']:.1f} p0 id R/G={100*s['raven_p0_identity']['mean']:.3f}/{100*s['gnnome_p0_identity']['mean']:.3f} "
          f"dID={s['d_p0_identity_pp']['mean']:+.3f}pp [{s['d_p0_identity_pp']['ci_lo']:+.3f},{s['d_p0_identity_pp']['ci_hi']:+.3f}] ({s['d_p0_identity_pp']['gnnome_better_in']}) "
          f"racon2 err R/G={s['raven_racon2_err']['mean']:.2f}/{s['gnnome_racon2_err']['mean']:.2f} dErr={s['d_racon2_err']['mean']:+.2f} ({s['d_racon2_err']['gnnome_better_in']})")

# ------------------------------------------------------------- F8
fig, axs = plt.subplots(2, 3, figsize=(14, 6.8)); fig.subplots_adjust(wspace=0.62, hspace=0.55); ax = axs.ravel()
series = [r for r in rows if r["cond"] in ORDER[:4] or r["cond"] == "run1reads_d30"]
other = [r for r in rows if r["cond"] in ("acc085_d60_L9k", "acc085_d30_L18k", "acc099_d60_L18k")]
xs = lambda R: [100 * r["read_identity"] for r in R]
for key, c, mk, lab in (("raven", F.RAVEN, "o", "Raven"), ("gnnome", F.GNN, "D", "GNNome")):
    ax[0].scatter(xs(series), [100 * r[f"{key}_p0_identity"] for r in series], color=c, marker=mk, s=26, edgecolor=F.SURF, linewidth=0.6, label=lab, zorder=3)
    ax[1].scatter(xs(series), [100 * r[f"{key}_source_read_error"] for r in series], color=c, marker=mk, s=26, edgecolor=F.SURF, linewidth=0.6, label=lab, zorder=3)
    ax[2].scatter(xs(series), [r[f"{key}_racon2_err"] + 0.01 for r in series], color=c, marker=mk, s=26, edgecolor=F.SURF, linewidth=0.6, label=lab, zorder=3)
ax[1].scatter(xs(series), [100 * r["all_reads_error_mean"] for r in series], color=F.MUTED, marker="_", s=60, label="all reads (mean)", zorder=2)
ax[0].plot([84, 100], [84, 100], color=F.MUTED, lw=0.8, ls="--"); ax[0].text(86.5, 85.2, "y = x", fontsize=6.5, color=F.INK2)
ax[0].set(xlabel="Read identity (%)", ylabel="Unpolished assembly identity (%)", title="Unpolished assembly accuracy"); ax[0].legend(loc="upper left")
ax[1].set_yscale("log"); ax[1].set(xlabel="Read identity (%)", ylabel="Source-read error, length-weighted (%)", title="Accuracy of the reads copied"); ax[1].legend(loc="lower left", fontsize=6.5)
ax[2].set_yscale("log"); ax[2].set(xlabel="Read identity (%)", ylabel="Errors per 100 kbp + 0.01 (log)", title="After Racon ×2"); ax[2].legend(loc="upper right")
ax[3].scatter(xs(series + other), [r["gnn_logit_true_median"] for r in series + other], color=F.INK2, s=24, zorder=3)
ax[3].axhline(0, color=F.MUTED, lw=0.8, ls="--"); ax[3].text(86, 0.8, "p = 0.5", fontsize=6.5, color=F.INK2)
ax[3].set(xlabel="Read identity (%)", ylabel="Median GNNome logit (true edges)", title="GNNome score calibration")
ax[4].scatter(xs(series + other), [r["edges_spurious"] for r in series + other], color=F.SPUR, s=24, zorder=3, label="spurious edges in graph")
ax[4].scatter(xs(series + other), [r["walk_spurious_or_backward"] for r in series + other], color=F.INK, marker="x", s=24, zorder=4, label="… used by GNNome walk")
ax[4].set(xlabel="Read identity (%)", ylabel="Edges", title="Spurious overlaps (ground truth)"); ax[4].legend(loc="upper left", fontsize=6.5)
labs = [NAME[c].split(" (")[0] for c in conds]
for i, c in enumerate(conds):
    v = [r["d_p0_identity_pp"] for r in rows if r["cond"] == c]
    ax[5].scatter(v, [i] * len(v), color=F.MUTED, s=16, zorder=2)
    m = summary[c]["d_p0_identity_pp"]
    if not math.isnan(m["ci_lo"]): ax[5].plot([m["ci_lo"], m["ci_hi"]], [i, i], color=F.GNN, lw=2)
    ax[5].plot([m["mean"]], [i], marker="D", ms=6, color=F.GNN, markeredgecolor=F.SURF, zorder=3)
ax[5].axvline(0, color=F.MUTED, lw=0.8); ax[5].set_yticks(range(len(conds))); ax[5].set_yticklabels(labs, fontsize=6.5); ax[5].invert_yaxis()
ax[5].set(xlabel="GNNome − Raven unpolished identity (pp)", title="Paired difference (mean, 95 % CI)")
for a, l in zip(ax, "abcdef"): F.panel(a, l)
F.save(fig, "F8_sensitivity")

# ------------------------------------------------------------- F9 forest (polished)
fig, ax = plt.subplots(1, 2, figsize=(11, 0.45 * len(conds) + 1.4)); fig.subplots_adjust(wspace=0.08)
for j, (d, xl) in enumerate((("d_p0_identity_pp", "Unpolished identity: GNNome − Raven (pp)"), ("d_racon2_err", "After Racon ×2 errors/100 kbp: GNNome − Raven"))):
    a = ax[j]
    for i, c in enumerate(conds):
        v = [r[d] for r in rows if r["cond"] == c]; m = summary[c][d]
        a.scatter(v, [i] * len(v), color=F.MUTED, s=14, zorder=2)
        if not math.isnan(m["ci_lo"]): a.plot([m["ci_lo"], m["ci_hi"]], [i, i], color=F.GNN, lw=2)
        a.plot([m["mean"]], [i], marker="D", ms=6, color=F.GNN, markeredgecolor=F.SURF, zorder=3)
    a.axvline(0, color=F.INK2, lw=0.8); a.set_yticks(range(len(conds))); a.invert_yaxis(); a.set_xlabel(xl)
    a.set_yticklabels([NAME[c] for c in conds] if j == 0 else [], fontsize=7)
    a.text(0.98 if j == 0 else 0.02, 0.02, "GNNome better →" if j == 0 else "← GNNome better", transform=a.transAxes, ha="right" if j == 0 else "left", va="bottom", fontsize=7, color=F.INK2)
fig.suptitle("Paired GNNome − Raven differences per condition (orange: mean ± 95 % CI; grey: individual seeds)", x=0.12, ha="left", y=1.02, fontsize=9, fontweight="bold")
F.save(fig, "F9_forest")
