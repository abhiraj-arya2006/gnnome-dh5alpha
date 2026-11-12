"""Collect one row of metrics per sensitivity replicate (task 6) into a TSV.
Recomputes read truth and edge labels with the final script versions so every replicate is analysed identically.
usage: collect_sensitivity.py <runs_dir> <out.tsv>
"""
import sys, os, json, glob, csv, subprocess
runs, out = sys.argv[1:3]
PY = "/home/abhiraj_arya/miniconda3/envs/gnnome/bin/python"; S = "/home/abhiraj_arya/research/analysis/scripts"
def quast(p):
    try: return {r[0]: r[1] for r in csv.reader(open(p), delimiter="\t") if len(r) > 1}
    except FileNotFoundError: return {}
def num(x):
    try: return float(x)
    except (TypeError, ValueError): return None
rows = []
for d in sorted(glob.glob(f"{runs}/*_s*")):
    if not os.path.exists(f"{d}/DONE"): continue
    p = dict(kv.split("=") for kv in open(f"{d}/params.txt").read().split())
    subprocess.run([PY, "-I", f"{S}/read_truth.py", f"{d}/sim/sim_0001.maf.gz", f"{d}/sim/sim_0001.fq.gz", "4583637", f"{d}/sim/truth"],
                   capture_output=True, check=True)
    subprocess.run([PY, "-I", f"{S}/graph_truth.py", f"{d}/raven/graph_1.gfa", f"{d}/gnnome/raven", f"{d}/sim/truth.reads.tsv",
                    f"{d}/gnnome/edges"], capture_output=True, check=True)
    rt = json.load(open(f"{d}/sim/truth.summary.json")); gt = json.load(open(f"{d}/gnnome/edges.summary.json"))
    cx, tr = gt["complexity"], gt["truth"]
    ident = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(f"{d}/identity_unpolished.tsv")}
    row = dict(replicate=os.path.basename(d), cond=p["cond"], accuracy=float(p["acc"]), depth=float(p["depth"]),
               len_mean=float(p["len_mean"]), seed=int(p["seed"]),
               read_identity=rt["identity_aggregate"], read_coverage=rt["coverage"], read_len_n50=rt["len_n50"], reads=rt["reads"],
               graph_nodes=cx["nodes"], graph_edges=cx["edges"], edges_per_node=cx["edges_per_node"], frac_transitive=cx["frac_transitive"],
               weak_components=cx["weak_components"], edges_spurious=tr["edges_spurious"], edges_containment=tr["edges_containment"],
               gnn_auroc=tr["auroc_true_vs_spurious"], gnn_ap=tr["ap_true_vs_spurious"], gnn_prob_true_median=tr["prob_true_median"],
               gnn_logit_true_median=tr["logit_true_median"], walk_edges=tr["walk_edges"],
               walk_spurious_or_backward=len(tr["walk_spurious_or_backward_edges"]), walk_nontrue=tr["walk_nontrue_edges"],
               walk_genome_cov=tr["walk_genome_covered_frac"])
    # assembly length / contig count from the FASTA (identity_unpolished.tsv kept only the last PAF record's length)
    fa = {"raven": f"{d}/raven/raven_p0.fasta", "gnnome": f"{d}/gnnome/raven/assembly/0_assembly.fasta"}
    seqs = {}
    for a, p in fa.items():
        s, cur = [], []
        for l in open(p):
            if l.startswith(">"):
                if cur: s.append("".join(cur)); cur = []
            else: cur.append(l.strip())
        if cur: s.append("".join(cur))
        seqs[a] = s
    # read selection: length-weighted per-read error of the reads each assembler copied (exact-substring layout)
    T = {}
    for i, l in enumerate(open(f"{d}/sim/truth.reads.tsv")):
        if i: f = l.split("\t"); T[f[0]] = 1 - float(f[9])
    rcx = lambda s: s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    segs = {}; n2r = {}; idx = 0; olx = {}; r2n = {}
    for line in open(f"{d}/raven/graph_1.gfa"):
        f = line.rstrip("\n").split("\t")
        if f[0] == "S": segs[idx] = f[2]; segs[idx + 1] = rcx(f[2]); n2r[idx] = n2r[idx + 1] = f[1]; r2n[f[1]] = (idx, idx + 1); idx += 2
        elif f[0] == "L":
            a_, o1, b_, o2, cig = f[1:6]; o = int(cig[:-1]); A_, B_ = r2n[a_], r2n[b_]
            olx[(A_[0] if o1 == "+" else A_[1], B_[0] if o2 == "+" else B_[1])] = o; olx[(B_[1] if o2 == "+" else B_[0], A_[1] if o1 == "+" else A_[0])] = o
    import pickle
    w = pickle.load(open(f"{d}/gnnome/raven/decode/0_walks.pkl", "rb"))[0]
    gc = [(n2r[u], max(0, len(segs[u]) - olx[(u, v)])) for u, v in zip(w[:-1], w[1:])] + [(n2r[w[-1]], len(segs[w[-1]]))]
    ctg = "".join(seqs["raven"]); hits = []
    for n, s in segs.items():
        k = ctg.find(s[:1000])
        if k >= 0 and ctg.find(s[:1000], k + 1) < 0: hits.append((k, n))
    hits.sort(); rl = [(n2r[n], (hits[i + 1][0] if i + 1 < len(hits) else len(ctg)) - k) for i, (k, n) in enumerate(hits)]
    for a, contrib in (("gnnome", gc), ("raven", rl)):
        L = sum(c for _, c in contrib)
        row[f"{a}_source_read_error"] = sum(c * T[r] for r, c in contrib) / L if L else None
        row[f"{a}_source_reads"] = len(contrib)
    row["all_reads_error_mean"] = sum(T.values()) / len(T)
    for a in ("raven", "gnnome"):
        row[f"{a}_p0_identity"] = float(ident[a][3]); row[f"{a}_p0_len"] = sum(len(x) for x in seqs[a]); row[f"{a}_p0_contigs"] = len(seqs[a])
        row[f"{a}_p0_refcov"] = int(ident[a][2])
        q0 = quast(f"{d}/quast_p0_minid80/{a}/report.tsv"); q2 = quast(f"{d}/quast_racon2/{a}/report.tsv")
        for k, n in (("# contigs", "contigs"), ("Genome fraction (%)", "gf"), ("# mismatches per 100 kbp", "mm"), ("# indels per 100 kbp", "indel"),
                     ("# misassemblies", "mis"), ("# local misassemblies", "localmis"), ("NGA50", "nga50"), ("Unaligned length", "unaligned")):
            row[f"{a}_p0m80_{n}"] = num(q0.get(k)); row[f"{a}_racon2_{n}"] = num(q2.get(k))
    rows.append(row)
with open(out, "w") as w:
    k = list(rows[0].keys()); w.write("\t".join(k) + "\n")
    for r in rows: w.write("\t".join("" if r[c] is None else str(r[c]) for c in k) + "\n")
print(f"{len(rows)} replicates -> {out}")
