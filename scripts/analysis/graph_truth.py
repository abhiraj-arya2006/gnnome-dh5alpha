"""Graph complexity + ground-truth edge labels + GNNome prediction analysis.

usage: graph_truth.py <gfa> <gnnome_run_dir (…/raven)> <reads.tsv or -> <out_prefix> [genome_len]

Node convention (identical to graph_parser.only_from_gfa): the i-th S line gives
node 2i (read in FASTQ orientation, which is the GFA orientation) and node 2i+1
(reverse complement). Every L line gives two directed edges.

Truth (needs PBSIM read positions from read_truth.py):
  node genome strand = read strand ('+' = FASTQ read equals forward genome) for
  even nodes, flipped for odd nodes. Edge u->v is TRUE when both nodes lie on the
  same genome strand, their true genomic intervals intersect, and v continues u in
  the direction of travel (tolerance TOL bp). Otherwise FALSE with a reason:
  strand_conflict, no_genomic_overlap (reads from different loci, i.e. a
  repeat-induced/spurious overlap), or wrong_direction.
"""
import sys, json, pickle, math, statistics as st
from collections import Counter, defaultdict
import torch, dgl, networkx as nx
from sklearn.metrics import roc_auc_score, average_precision_score

gfa, run, truth_path, out = sys.argv[1:5]
G = int(sys.argv[5]) if len(sys.argv) > 5 else 4583637
TOL = 0  # strict geometry: noisy-read coordinates do not justify a tolerance

# ---- parse GFA exactly like graph_parser ----
r2n, n2r, seglen, edges_ol = {}, {}, {}, {}
idx = 0
for line in open(gfa):
    f = line.rstrip("\n").split("\t")
    if f[0] == "S":
        r2n[f[1]] = (idx, idx + 1); n2r[idx] = n2r[idx + 1] = f[1]; seglen[f[1]] = len(f[2]); idx += 2
    elif f[0] == "L":
        a, o1, b, o2, cig = f[1:6]
        ol = int(cig[:-1])
        if ol == 0: continue
        A, B = r2n[a], r2n[b]
        sr = A[0] if o1 == "+" else A[1]; dr = B[0] if o2 == "+" else B[1]
        sv = B[1] if o2 == "+" else B[0]; dv = A[1] if o1 == "+" else A[0]
        edges_ol[(sr, dr)] = ol; edges_ol[(sv, dv)] = ol

g = dgl.load_graphs(f"{run}/processed/0.dgl")[0][0]
u, v = [t.tolist() for t in g.edges()]
E = list(zip(u, v))
assert set(E) == set(edges_ol), "DGL edges do not match the GFA"
logit = torch.load(f"{run}/decode/0_predicts.pt").squeeze().tolist()
assert len(logit) == len(E)
prob = [1 / (1 + math.exp(-x)) for x in logit]
walks = pickle.load(open(f"{run}/decode/0_walks.pkl", "rb"))
walk_edges = set()
for w in walks:
    walk_edges |= set(zip(w[:-1], w[1:]))

# ---- complexity ----
D = nx.DiGraph(); D.add_nodes_from(range(g.num_nodes())); D.add_edges_from(E)
outd = [D.out_degree(n) for n in D]; ind = [D.in_degree(n) for n in D]
wcc = sorted((len(c) for c in nx.weakly_connected_components(D)), reverse=True)
trans = sum(1 for (a, c) in E if set(D.successors(a)) & set(D.predecessors(c)))
ols = list(edges_ol.values())
cx = dict(gfa=gfa, segments=len(r2n), links=len(E) // 2, nodes=g.num_nodes(), edges=len(E),
          edges_per_node=len(E) / g.num_nodes(),
          out_degree_mean=st.mean(outd), out_degree_max=max(outd),
          frac_nodes_outdeg_gt1=sum(d > 1 for d in outd) / len(outd),
          frac_nodes_outdeg_gt3=sum(d > 3 for d in outd) / len(outd),
          nodes_outdeg0=sum(d == 0 for d in outd), nodes_indeg0=sum(d == 0 for d in ind),
          weak_components=len(wcc), largest_component=wcc[0], isolated_nodes=sum(1 for c in wcc if c == 1),
          transitive_edges=trans, frac_transitive=trans / len(E),
          cyclomatic_number=len(E) - g.num_nodes() + len(wcc),
          overlap_len_median=st.median(ols), overlap_len_mean=st.mean(ols), overlap_len_min=min(ols),
          segment_len_median=st.median(seglen.values()))
outdeg_hist = Counter(outd)

# ---- truth ----
res = dict(complexity=cx, outdeg_hist=dict(sorted(outdeg_hist.items())))
rows = []
if truth_path != "-":
    T = {}
    for i, l in enumerate(open(truth_path)):
        if i == 0: continue
        f = l.split("\t"); T[f[0]] = (f[1], int(f[2]), int(f[3]))
    def node_pos(n):
        s, a, b = T[n2r[n]]
        if n % 2: s = "-" if s == "+" else "+"
        return s, a, b
    reason = Counter()
    for k, (a_, b_) in enumerate(E):
        su, au, bu = node_pos(a_); sv, av, bv = node_pos(b_)
        inter = min(bu, bv) - max(au, av)
        fwd = (av > au and bv > bu) if su == "+" else (av < au and bv < bu)
        if su != sv: lab, why = 0, "strand_conflict"
        elif inter <= 0: lab, why = 0, "no_genomic_overlap"
        elif fwd: lab, why = 1, "true"
        elif av <= au and bv >= bu: lab, why = 2, "containment_src_in_dst"
        elif av >= au and bv <= bu: lab, why = 2, "containment_dst_in_src"
        else: lab, why = 0, "wrong_direction"
        reason[why] += 1
        gap = max(au, av) - min(bu, bv) if inter <= 0 else 0
        rows.append((k, a_, b_, n2r[a_], n2r[b_], su, au, bu, sv, av, bv, edges_ol[(a_, b_)], inter, gap,
                     lab, why, logit[k], prob[k], int((a_, b_) in walk_edges)))
    # primary metric: true (1) vs spurious/false (0); containment edges (2) excluded
    prim = [r for r in rows if r[14] in (0, 1)]
    y = [r[14] for r in prim]; s = [r[16] for r in prim]
    # secondary metric: containment counted as negative
    y2 = [1 if r[14] == 1 else 0 for r in rows]; s2 = [r[16] for r in rows]
    pos = [r[17] for r in rows if r[14] == 1]; neg = [r[17] for r in rows if r[14] != 1]
    tp = sum(1 for r in rows if r[17] >= 0.5 and r[14] == 1); fp = sum(1 for r in rows if r[17] >= 0.5 and r[14] != 1)
    fn = sum(1 for r in rows if r[17] < 0.5 and r[14] == 1); tn = sum(1 for r in rows if r[17] < 0.5 and r[14] != 1)
    walk_rows = [r for r in rows if r[18]]
    walk_false = [r for r in walk_rows if r[14] != 1]
    # genome covered by the walk's reads (true positions)
    cov = bytearray(G)
    for w in walks:
        for n in w:
            _, a, b = node_pos(n)
            cov[a:b] = b"\x01" * (b - a)
    # false edges that are "attractive": scored above the best true alternative from the same source
    best_true = defaultdict(lambda: -1e9)
    for r in rows:
        if r[14] == 1: best_true[r[1]] = max(best_true[r[1]], r[16])
    false_beats_true = sum(1 for r in rows if r[14] == 0 and r[1] in best_true and r[16] > best_true[r[1]])
    ok = lambda yy: 0 < sum(yy) < len(yy)
    res["truth"] = dict(
        tolerance_bp=TOL, edges_true=sum(r[14] == 1 for r in rows), edges_spurious=sum(r[14] == 0 for r in rows),
        edges_containment=sum(r[14] == 2 for r in rows), frac_spurious=sum(r[14] == 0 for r in rows) / len(rows),
        edge_classes=dict(reason), nodes_with_spurious_out_edge=len({r[1] for r in rows if r[14] == 0}),
        auroc_true_vs_spurious=roc_auc_score(y, s) if ok(y) else None,
        ap_true_vs_spurious=average_precision_score(y, s) if ok(y) else None,
        ap_spurious_class=average_precision_score([1 - t for t in y], [-x for x in s]) if ok(y) else None,
        auroc_true_vs_spurious_or_containment=roc_auc_score(y2, s2) if ok(y2) else None,
        prob_true_median=st.median(pos), prob_nontrue_median=st.median(neg) if neg else None,
        prob_true_p05=sorted(pos)[int(.05 * len(pos))], prob_true_p95=sorted(pos)[int(.95 * len(pos))],
        logit_true_median=st.median(r[16] for r in rows if r[14] == 1),
        confusion_at_0p5=dict(tp=tp, fp=fp, fn=fn, tn=tn),
        false_edges_scored_above_best_true_sibling=false_beats_true,
        walks=len(walks), walk_nodes=sum(len(w) for w in walks), walk_edges=len(walk_rows),
        walk_edge_classes=dict(Counter(r[15] for r in walk_rows)),
        walk_nontrue_edges=len(walk_false),
        walk_spurious_or_backward_edges=[dict(src=r[3], dst=r[4], why=r[15], src_pos=f"{r[5]}{r[6]}-{r[7]}", dst_pos=f"{r[8]}{r[9]}-{r[10]}",
                                      overlap=r[11], prob=round(r[17], 4)) for r in walk_false if r[14] == 0],
        walk_genome_covered_bp=sum(cov), walk_genome_covered_frac=sum(cov) / G)
    with open(out + ".edges.tsv", "w") as w:
        w.write("edge_id\tsrc\tdst\tsrc_read\tdst_read\tsrc_strand\tsrc_start\tsrc_end\tdst_strand\tdst_start\tdst_end\t"
                "gfa_overlap\tgenomic_intersection\tgenomic_gap\tlabel\treason\tlogit\tprob\tin_walk\n")
        for r in rows:
            w.write("\t".join(str(round(x, 5)) if isinstance(x, float) else str(x) for x in r) + "\n")
json.dump(res, open(out + ".summary.json", "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != "outdeg_hist"}, indent=1)[:6000])
