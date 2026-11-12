"""Publication figures (task 11). Static PNG (300 dpi) + PDF, light surface.

Colour roles (validated with validate_palette.py, a port of the dataviz skill's validator):
  Raven #2a78d6, GNNome #eb6834, third series #1baf7a  -> all-pairs PASS (scatter-safe)
  Run 1 #1baf7a, Run 2 #6250d6                          -> all-pairs PASS (used only where Raven blue is absent)
  edge classes: true = recessive ink, containment #eda100, spurious #e34948 (+ direct labels; contrast relief)
Identity is never colour-alone: every multi-series panel has a legend and/or direct labels.
usage: figures.py [fig ...]   (default: all available)
"""
import sys, json, csv, gzip, pickle, os, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

A = "/home/abhiraj_arya/research/analysis"; OUT = f"{A}/11_figures"; R2 = "/home/abhiraj_arya/research/run2_corrected"
G1 = "/home/abhiraj_arya/research/GNNome"
RAVEN, GNN, THIRD = "#2a78d6", "#eb6834", "#1baf7a"
RUN1, RUN2 = "#1baf7a", "#6250d6"
CONT, SPUR = "#eda100", "#e34948"
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.labelsize": 8, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5, "axes.axisbelow": True,
    "lines.linewidth": 2.0, "legend.frameon": False, "legend.fontsize": 7.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "pdf.fonttype": 42})
def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=300, bbox_inches="tight"); fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight"); plt.close(fig)
    print("wrote", name)
def tsv(p): return list(csv.DictReader(open(p), delimiter="\t"))
def label_end(ax, x, y, text, color):
    ax.annotate(text, (x, y), xytext=(4, 0), textcoords="offset points", va="center", fontsize=7.5, color=INK2)
def panel(ax, letter):
    ax.text(-0.14, 1.06, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", color=INK, va="bottom")

# ---------------------------------------------------------------- F1 reads
def fig1():
    r = {k: tsv(f"{A}/03_reads/{k}.reads.tsv") for k in ("run1", "run2")}
    fig, ax = plt.subplots(1, 4, figsize=(12, 2.8)); fig.subplots_adjust(wspace=0.45)
    for k, c, lab in (("run1", RUN1, "Run 1"), ("run2", RUN2, "Run 2")):
        L = np.array([int(x["read_len"]) for x in r[k]]) / 1000
        ax[0].hist(L, bins=np.arange(0, 75, 1), histtype="step", color=c, lw=1.6, label=lab)
        I = np.array([float(x["identity"]) for x in r[k]]) * 100
        ax[1].hist(I, bins=np.arange(60, 100.25, 0.25), histtype="step", color=c, lw=1.6, label=lab)
    ax[0].set(xlabel="Read length (kb)", ylabel="Reads", title="Read length"); ax[0].legend(loc="upper right")
    ax[1].set(xlabel="Per-read identity to reference (%)", ylabel="Reads", title="Read accuracy"); ax[1].set_yscale("log")
    ax[1].legend(loc="upper left")
    # error vs position
    prof = {}
    for line in open(f"{A}/03_reads/error_position_profile.txt"):
        if line.startswith("Run"): cur = line.split()[0]; prof[cur] = []
        elif "from 5-prime" in line:
            f = line.split(); prof[cur].append((f[0], float(f[3]), float(f[6])))
    for k, c, lab in (("Run1", RUN1, "Run 1"), ("Run2", RUN2, "Run 2")):
        y = np.array([p[1] for p in prof[k]]); y = y / y[-1]
        ax[2].plot(range(len(y)), y, color=c, marker="o", ms=4, label=lab)
    ax[2].set_xticks(range(len(prof["Run1"]))); ax[2].set_xticklabels([p[0].replace("-end", "+") for p in prof["Run1"]], rotation=45, ha="right")
    ax[2].set(xlabel="Distance from read 5' end (bp)", ylabel="Error rate / rate >10 kb", title="Errors along the read"); ax[2].legend()
    for k, c, lab in (("run1", RUN1, "Run 1"), ("run2", RUN2, "Run 2")):
        cv = tsv(f"{A}/03_reads/{k}.coverage_1kb.tsv")
        x = np.array([int(z["bin_start"]) for z in cv]) / 1e6; y = np.array([float(z["mean_depth"]) for z in cv])
        k10 = np.convolve(y, np.ones(10) / 10, mode="same")
        ax[3].plot(x, k10, color=c, lw=1.2, label=lab)
    ax[3].set(xlabel="Genome position (Mb)", ylabel="True depth (10-kb mean)", title="True coverage", ylim=(0, None)); ax[3].legend(loc="lower center", ncol=2)
    for a, l in zip(ax, "abcd"): panel(a, l)
    save(fig, "F1_reads")

# ---------------------------------------------------------------- F2 graphs
def fig2():
    S = {k: json.load(open(p)) for k, p in (("Run 1", f"{A}/05_edge_predictions/run1.summary.json"),
                                            ("Run 2", f"{A}/05_edge_predictions/run2.summary.json"),
                                            ("K-12 control", f"{A}/04_graphs/k12_control.summary.json"))}
    fig, ax = plt.subplots(1, 4, figsize=(12, 2.7)); fig.subplots_adjust(wspace=0.42)
    for a, (k, s) in zip(ax[:3], S.items()):
        h = {int(x): y for x, y in s["outdeg_hist"].items()}; xs = sorted(h)
        a.bar(xs, [h[x] for x in xs], width=0.8, color=INK2 if k != "K-12 control" else MUTED)
        c = s["complexity"]
        a.set(xlabel="Out-degree", ylabel="Nodes", title=f"{k}: out-degree")
        a.text(0.02, 0.97, f"{c['nodes']:,} nodes\n{c['edges']:,} edges\n{c['edges_per_node']:.1f} edges/node\n{100*c['frac_transitive']:.0f}% transitive",
               transform=a.transAxes, ha="left", va="top", fontsize=7, color=INK2)
        a.set_ylim(0, max(h.values()) * 1.55)
    cls = [("true", INK2), ("containment", CONT), ("spurious", SPUR)]
    for i, k in enumerate(("Run 1", "Run 2")):
        t = S[k]["truth"]; vals = [t["edges_true"], t["edges_containment"], t["edges_spurious"]]
        for j, ((n, col), v) in enumerate(zip(cls, vals)):
            ax[3].bar(i * 4 + j, v, color=col, width=0.85)
            ax[3].text(i * 4 + j, v * 1.25, f"{v:,}", ha="center", fontsize=6.5, color=INK2)
    ax[3].set_yscale("log"); ax[3].set_xticks([1, 5]); ax[3].set_xticklabels(["Run 1", "Run 2"]); ax[3].set_ylim(1, 4e5)
    ax[3].set(ylabel="Directed edges (log)", title="Ground-truth edge classes")
    ax[3].legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c) for _, c in cls], labels=[n for n, _ in cls], loc="upper right", fontsize=7)
    for a, l in zip(ax, "abcd"): panel(a, l)
    save(fig, "F2_graph_complexity")

# ---------------------------------------------------------------- F3 GNNome predictions
def fig3():
    from sklearn.metrics import roc_curve
    fig, ax = plt.subplots(1, 3, figsize=(10, 2.8)); fig.subplots_adjust(wspace=0.35)
    for a, run, lab in ((ax[0], "run1", "Run 1"), (ax[1], "run2", "Run 2")):
        E = tsv(f"{A}/05_edge_predictions/{run}.edges.tsv")
        lg = {c: np.array([float(e["logit"]) for e in E if (e["reason"] == "true" if c == "true" else
                         e["reason"].startswith("containment") if c == "containment" else e["reason"] in ("no_genomic_overlap", "wrong_direction", "strand_conflict"))])
              for c in ("true", "containment", "spurious")}
        lo = min(v.min() for v in lg.values() if len(v)); hi = max(v.max() for v in lg.values() if len(v))
        bins = np.linspace(lo, hi, 60)
        a.hist(lg["true"], bins=bins, color=INK2, alpha=0.85, label=f"true (n={len(lg['true']):,})")
        a.hist(lg["containment"], bins=bins, color=CONT, label=f"containment (n={len(lg['containment'])})")
        for x in lg["spurious"]: a.axvline(x, color=SPUR, lw=1.4, ymax=0.25)
        a.plot([], [], color=SPUR, lw=1.4, label=f"spurious (n={len(lg['spurious'])}, ticks)")
        a.set_yscale("log"); a.axvline(0, color=MUTED, lw=0.8, ls="--")
        a.text(0.0, 1.0, " p = 0.5", transform=a.get_xaxis_transform(), fontsize=6.5, color=MUTED, va="top")
        a.set(xlabel="GNNome edge logit", ylabel="Edges (log)", title=f"{lab}: GNNome edge scores"); a.legend(loc="upper left", fontsize=6.5)
        y = [1] * len(lg["true"]) + [0] * len(lg["spurious"]); s = list(lg["true"]) + list(lg["spurious"])
        fpr, tpr, _ = roc_curve(y, s); auc = json.load(open(f"{A}/05_edge_predictions/{run}.summary.json"))["truth"]["auroc_true_vs_spurious"]
        ax[2].plot(fpr, tpr, color=RUN1 if run == "run1" else RUN2, drawstyle="steps-post", label=f"{lab}  AUROC {auc:.4f}")
    ax[2].plot([0, 1], [0, 1], color=MUTED, lw=0.8, ls="--")
    ax[2].set(xlabel="False-positive rate (spurious edges)", ylabel="True-positive rate (true edges)", title="True vs spurious edges", xlim=(-0.02, 1.02), ylim=(0, 1.02))
    ax[2].legend(loc="lower right")
    for a, l in zip(ax, "abc"): panel(a, l)
    save(fig, "F3_gnnome_edge_predictions")

# ---------------------------------------------------------------- F4 junction mechanism
def fig4():
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.0)); fig.subplots_adjust(wspace=0.42)
    J = tsv(f"{A}/09_repeat_errors/R1_GNNome.junctions.tsv")
    x = np.array([int(j["discrepancy"]) for j in J]); y = np.array([int(j["observed_net_indel_within_50bp"]) for j in J])
    rng = np.random.default_rng(0); jit = lambda n: rng.uniform(-0.25, 0.25, n)
    ax[0].scatter(x + jit(len(x)), y + jit(len(y)), s=9, color=GNN, edgecolor=SURF, linewidth=0.4)
    lim = (min(x.min(), y.min()) - 3, max(x.max(), y.max()) + 3); ax[0].plot(lim, lim, color=MUTED, lw=0.8, ls="--")
    r = json.load(open(f"{A}/09_repeat_errors/R1_GNNome.junctions.summary.json"))
    ax[0].text(0.03, 0.97, f"394 GNNome junctions, Run 1\nPearson r = {r['pearson_r']:.3f}\nexact match {r['predicted_equals_observed_net_indel']}",
               transform=ax[0].transAxes, va="top", fontsize=7, color=INK2)
    ax[0].set(xlabel="Predicted cut error from read truth (bp)", ylabel="Observed net indel at junction (bp)", title="Run 1: junction indel = cut error")
    printed_sets = {}
    for a, run, gfa, lab in ((ax[1], "R1", f"{G1}/data/graph_1.gfa", "Run 1"), (ax[2], "R2", f"{R2}/02_raven_p0/graph_1.gfa", "Run 2")):
        r2n = {}; idx = 0; pr = set()
        for line in open(gfa):
            f = line.rstrip("\n").split("\t")
            if f[0] == "S": r2n[f[1]] = (idx, idx + 1); idx += 2
            elif f[0] == "L":
                aa, o1, b, o2 = f[1:5]; A_, B_ = r2n[aa], r2n[b]; pr.add((A_[0] if o1 == "+" else A_[1], B_[0] if o2 == "+" else B_[1]))
        rows = tsv(f"{A}/09_repeat_errors/{run}_edge_cut_error.tsv")
        e_w = np.array([abs(int(q["error"])) for q in rows if (int(q["src"]), int(q["dst"])) in pr])
        e_t = np.array([abs(int(q["error"])) for q in rows if (int(q["src"]), int(q["dst"])) not in pr])
        cap = 400; bins = np.concatenate([np.arange(0, 21), [30, 50, 100, 200, cap]])
        for e, c, n in ((e_w, RAVEN, "written in GFA (Raven's own length)"), (e_t, GNN, "reverse-complement twin (length reused)")):
            cdf = [np.mean(np.minimum(e, cap) <= b) for b in bins]
            a.plot(bins, cdf, color=c, drawstyle="steps-post", label=n)
        a.set_xscale("symlog", linthresh=5); a.set_xlim(0, cap)
        a.set(xlabel="|cut error| of true edges (bp)", ylabel="Cumulative fraction", title=f"{lab}: overlap-length accuracy", ylim=(0, 1.02))
        a.axvline(5, color=MUTED, lw=0.8, ls=":")
        a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), fontsize=6.5, ncol=1)
    for a, l in zip(ax, "abc"): panel(a, l)
    save(fig, "F4_junction_mechanism")

# ---------------------------------------------------------------- F5 read selection
def fig5():
    fig, ax = plt.subplots(1, 3, figsize=(12.5, 2.9)); fig.subplots_adjust(wspace=0.55)
    def layout(gfa, asm, walk_dir):
        rc = lambda s: s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
        seg = {}; n2r = {}; idx = 0
        for line in open(gfa):
            f = line.rstrip("\n").split("\t")
            if f[0] == "S": seg[idx] = f[2]; seg[idx + 1] = rc(f[2]); n2r[idx] = n2r[idx + 1] = f[1]; idx += 2
        ctg = "".join(l.strip() for l in open(asm) if not l.startswith(">"))
        rav = set()
        for n, s in seg.items():
            k = ctg.find(s[:1000])
            if k >= 0 and ctg.find(s[:1000], k + 1) < 0: rav.add(n2r[n])
        w = pickle.load(open(f"{walk_dir}/decode/0_walks.pkl", "rb"))[0]
        return rav, {n2r[n] for n in w}
    for a, run, gfa, asm, wd, lab in ((ax[0], "run1", f"{G1}/data/graph_1.gfa", f"{G1}/data/CP017100_assembly.fasta", f"{G1}/data/raven", "Run 1"),
                                       (ax[1], "run2", f"{R2}/02_raven_p0/graph_1.gfa", f"{R2}/02_raven_p0/raven_p0_assembly.fasta", f"{R2}/04_gnnome/raven", "Run 2")):
        T = {x["read"]: 100 * (1 - float(x["identity"])) for x in tsv(f"{A}/03_reads/{run}.reads.tsv")}
        rav, gnn = layout(gfa, asm, wd)
        data = [list(T.values()), [T[r] for r in rav], [T[r] for r in gnn]]
        bp = a.boxplot(data, vert=False, widths=0.55, showfliers=False, patch_artist=True, medianprops=dict(color=INK, lw=1.4))
        for p, c in zip(bp["boxes"], (MUTED, RAVEN, GNN)): p.set(facecolor=c, alpha=0.9, edgecolor=c)
        a.set_yticks([1, 2, 3]); a.set_yticklabels([f"all reads\n(n={len(data[0]):,})", f"Raven layout\n(n={len(data[1]):,})", f"GNNome walk\n(n={len(data[2]):,})"])
        a.set(xlabel="Per-read error (%)", title=f"{lab}: source-read quality")
        if run == "run1": a.set_xscale("log")
    # expected vs observed (Run 1)
    lines = open(f"{A}/09_repeat_errors/R1_contributing_read_quality.txt").read().splitlines()
    vals = {}
    for l in lines[:2]:
        f = l.split(); name = f[0]
        exp = float(l.split("expected inherited errors=")[1].split()[0].replace(",", "")); obs = float(l.split(")=")[1].replace(",", ""))
        vals[name] = (exp, obs)
    for i, (n, c) in enumerate((("Raven", RAVEN), ("GNNome", GNN))):
        e, o = vals[n]
        ax[2].bar(i * 3, e, color=c, alpha=0.45, width=0.9); ax[2].bar(i * 3 + 1, o, color=c, width=0.9)
        ax[2].text(i * 3, e + 120, f"{e:,.0f}", ha="center", fontsize=7, color=INK2); ax[2].text(i * 3 + 1, o + 120, f"{o:,.0f}", ha="center", fontsize=7, color=INK2)
    ax[2].set_xticks([0, 1, 3, 4]); ax[2].set_xticklabels(["exp.", "obs.", "exp.", "obs."], fontsize=7)
    ax[2].text(0.5, -0.2, "Raven", transform=ax[2].get_xaxis_transform(), ha="center", fontsize=7.5, color=INK2); ax[2].text(3.5, -0.2, "GNNome", transform=ax[2].get_xaxis_transform(), ha="center", fontsize=7.5, color=INK2)
    ax[2].set(ylabel="Small errors in unpolished contig", title="Run 1: expected vs observed errors")
    for a, l in zip(ax, "abc"): panel(a, l)
    save(fig, "F5_read_selection")

# ---------------------------------------------------------------- F6 genome map
def fig6():
    Gl = 4583637 / 1e6
    fig, ax = plt.subplots(6, 1, figsize=(11, 7.2), sharex=True, gridspec_kw=dict(height_ratios=[1.1, 1, 1, 1, 0.8, 0.8], hspace=0.75))
    rep = [l.rstrip("\n").split("\t") for l in open(f"{A}/07_repeats/repeats.bed")]
    classes = ["rRNA", "IS_transposase", "Rhs", "phage_related", "tRNA", "unannotated"]
    for i, c in enumerate(classes):
        for r in rep:
            if r[5] == c: ax[0].add_patch(plt.Rectangle((int(r[1]) / 1e6, i - 0.35), max((int(r[2]) - int(r[1])) / 1e6, 0.004), 0.7, color=INK2))
    ax[0].set_yticks(range(len(classes))); ax[0].set_yticklabels(["rRNA", "IS", "Rhs", "phage", "tRNA", "other"], fontsize=6.5)
    ax[0].set_ylim(-0.6, len(classes) - 0.4); ax[0].grid(False); ax[0].set_title("Repeat families (self-alignment ≥300 bp, ≥85 %; classes from GenBank)")
    def density(path, cats, w=20000):
        pos = [int(e["ref"]) for e in tsv(path) if e["cat"] in cats]
        h, b = np.histogram(pos, bins=np.arange(0, 4583637 + w, w)); return b[:-1] / 1e6, h / (w / 1e5)
    for a, run, title, cats in ((ax[1], "R1", "Run 1 unpolished: small errors per 100 kb", ("mismatch", "indel_le5")),
                                (ax[2], "R2", "Run 2 unpolished: small errors per 100 kb", ("mismatch", "indel_le5")),
                                (ax[3], "R2", "Run 2 polished (Racon ×2): residual errors per 100 kb", ("mismatch", "indel_le5", "indel_gt5"))):
        for lab, c, name in (("Raven", RAVEN, "Raven_p0" if "unpol" in title else "Raven_Racon2"), ("GNNome", GNN, "GNNome" if "unpol" in title else "GNNome_Racon2")):
            x, y = density(f"{A}/09_repeat_errors/{run}_{name}.events.tsv", cats)
            a.plot(x, y, color=c, lw=1.0, label=lab)
        a.set_title(title); a.set_ylim(0, None); a.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, fontsize=6.5)
    def ticks(a, items, title):
        for y, (lab, pts, c, mk) in enumerate(items):
            a.scatter([p / 1e6 for p in pts], [y] * len(pts), marker=mk, s=22, color=c, edgecolor=SURF, linewidth=0.4)
        a.set_yticks(range(len(items))); a.set_yticklabels([i[0] for i in items], fontsize=6.5); a.set_ylim(-0.7, len(items) - 0.3); a.set_title(title)
    big = lambda run, n: [int(e["ref"]) for e in tsv(f"{A}/09_repeat_errors/{run}_{n}.events.tsv") if e["cat"] in ("indel_gt5", "breakpoint")]
    ticks(ax[4], [("Raven", big("R1", "Raven_p0"), RAVEN, "o"), ("GNNome", big("R1", "GNNome"), GNN, "D")], "Run 1 unpolished: indels >5 bp and breakpoints")
    sp = []
    for run in ("run1", "run2"):
        for e in tsv(f"{A}/05_edge_predictions/{run}.edges.tsv"):
            if e["reason"] == "no_genomic_overlap": sp += [int(e["src_start"]), int(e["dst_start"])]
    ticks(ax[5], [("spurious overlaps", sp, SPUR, "v"), ("Raven tRNA collapse", [591774], RAVEN, "o")], "Graph-level repeat events")
    ax[5].annotate("Rhs copies", (3.63, 0), xytext=(3.9, 0.45), fontsize=6.5, color=INK2, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
    ax[5].set_xlabel("Genome position (Mb)"); ax[5].set_xlim(0, Gl)
    save(fig, "F6_genome_map")

FIGS = {"1": fig1, "2": fig2, "3": fig3, "4": fig4, "5": fig5, "6": fig6}

# ---------------------------------------------------------------- F7 polishing
def fig7():
    def q(p):
        d = {r[0]: r[1] for r in csv.reader(open(p), delimiter="\t") if len(r) > 1}
        return float(d["# mismatches per 100 kbp"]) + float(d["# indels per 100 kbp"])
    P1 = f"{A}/10_polishing/run1"
    run1 = {"Raven": [q(f"{G1}/data/classical_raven_quast_results/report.tsv"), q(f"{P1}/raven/quast_round1/report.tsv"), q(f"{P1}/raven/quast_round2/report.tsv")],
            "GNNome": [q(f"{G1}/data/raven/quast_results/report.tsv"), q(f"{P1}/gnnome/quast_round1/report.tsv"), q(f"{P1}/gnnome/quast_round2/report.tsv")]}
    run2 = {"Raven": [q(f"{R2}/05_quast/minid80_raven_p0/report.tsv"), q(f"{A}/10_polishing/run2/quast_raven_round1/report.tsv"), q(f"{R2}/07_polish/quast/raven_racon2/report.tsv")],
            "GNNome": [q(f"{R2}/05_quast/minid80_gnnome_p0/report.tsv"), q(f"{A}/10_polishing/run2/quast_gnnome_round1/report.tsv"), q(f"{R2}/07_polish/quast/gnnome_racon2/report.tsv")]}
    p2 = q(f"{R2}/07_polish/quast/raven_internal_p2/report.tsv")
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.1)); fig.subplots_adjust(wspace=0.62)
    for a, data, lab in ((ax[0], run1, "Run 1 (99.84 % reads, 60×)"), (ax[1], run2, "Run 2 (86.01 % reads, 30×)")):
        fmt = lambda y: f"{y:,.0f}" if y >= 100 else f"{y:.2f}"
        for (n, v), c, mk in zip(data.items(), (RAVEN, GNN), ("o", "D")):
            a.plot([0, 1, 2], v, color=c, marker=mk, ms=6, markeredgecolor=SURF, markeredgewidth=0.8,
                   label=f"{n}: " + " → ".join(fmt(y) for y in v))
        a.set_yscale("symlog", linthresh=0.1); a.set_xticks([0, 1, 2]); a.set_xticklabels(["unpolished", "Racon ×1", "Racon ×2"])
        a.set(ylabel="Errors per 100 kbp (QUAST)", title=lab); a.legend(loc="upper right", fontsize=7); a.set_xlim(-0.3, 2.5)
    ax[1].scatter([2.15], [p2], marker="s", s=36, color=RAVEN, facecolor=SURF, linewidth=1.4, zorder=3, label=f"Raven -p 2 (built-in): {p2:.2f}")
    ax[1].legend(loc="upper right", fontsize=7)
    ax[1].text(0.03, 0.04, "unpolished: QUAST\n--min-identity 80", transform=ax[1].transAxes, fontsize=6.3, color=INK2, ha="left")
    ax[0].text(0.03, 0.05, "after Racon ×2 both are\nsequence-identical (0 errors)", transform=ax[0].transAxes, fontsize=6.5, color=INK2, ha="left")
    S = json.load(open(f"{A}/10_polishing/run2/residual_errors.json"))
    cats = [("in homopolymer ≥4", "indels_in_homopolymer_ge4"), ("shared with other\nassembly", "shared_with_other"), ("in repeat", "in_repeat"),
            ("near GNNome\njunction (±50 bp)", "within_50bp_of_gnnome_junction")]
    for j, (n, c) in enumerate((("raven", RAVEN), ("gnnome", GNN))):
        tot = S[n]["events"]; vals = [100 * S[n][k] / tot for _, k in cats]
        ax[2].barh(np.arange(len(cats)) + (0.2 if j == 0 else -0.2), vals, height=0.38, color=c, label=n.replace("gnnome", "GNNome").replace("raven", "Raven"))
        for i, v in enumerate(vals): ax[2].text(v + 0.8, i + (0.2 if j == 0 else -0.2), f"{v:.1f}%", va="center", fontsize=6.5, color=INK2)
    hb = 100 * S["genome_bases_in_homopolymer_ge4"]
    ax[2].plot([hb, hb], [-0.45, 0.45], color=INK2, ls=":", lw=1.2)
    ax[2].annotate(f"genome: {hb:.1f}%", (hb, -0.45), xytext=(2, 2), textcoords="offset points", fontsize=6.3, color=INK2)
    ax[2].set_yticks(range(len(cats))); ax[2].set_yticklabels([c for c, _ in cats], fontsize=7); ax[2].invert_yaxis()
    ax[2].set(xlabel="% of residual errors after Racon ×2 (Run 2)", title="What polishing leaves behind", xlim=(0, 60)); ax[2].legend(loc="lower right")
    for a, l in zip(ax, "abc"): panel(a, l)
    save(fig, "F7_polishing")
FIGS["7"] = fig7

if __name__ == "__main__":
    for k in (sys.argv[1:] or FIGS): FIGS[k]()
