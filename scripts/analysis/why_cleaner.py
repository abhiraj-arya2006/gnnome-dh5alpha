"""Why does GNNome copy cleaner reads? Tests on the real inference graphs:
(1) overlap_similarity (a model input) tracks the true error of the two reads;
(2) the model's edge score rises with overlap_similarity;
(3) at every walk decision, the chosen edge vs. the alternative true edges."""
import sys, json, dgl, torch, numpy as np, csv
from scipy.stats import spearmanr, wilcoxon
dglp, edges_tsv, reads_tsv, out = sys.argv[1:5]
g = dgl.load_graphs(dglp)[0][0]
sim = g.edata['overlap_similarity'].numpy(); oll = g.edata['overlap_length'].numpy()
E = list(csv.DictReader(open(edges_tsv), delimiter="\t"))
assert len(E) == g.num_edges()
src, dst = [x.numpy() for x in g.edges()]
assert all(int(e["src"]) == src[i] and int(e["dst"]) == dst[i] for i, e in enumerate(E))
err = {r["read"]: 1 - float(r["identity"]) for r in csv.DictReader(open(reads_tsv), delimiter="\t")}
logit = np.array([float(e["logit"]) for e in E]); lab = np.array([int(e["label"]) for e in E]); inw = np.array([int(e["in_walk"]) for e in E])
pair_err = np.array([err[e["src_read"]] + err[e["dst_read"]] for e in E])
dst_err = np.array([err[e["dst_read"]] for e in E])
R = {}
t = lab == 1
R["spearman_sim_vs_pair_read_error_trueedges"] = spearmanr(sim[t], pair_err[t])[0]
R["spearman_logit_vs_sim_trueedges"] = spearmanr(logit[t], sim[t])[0]
R["spearman_logit_vs_overlap_len_trueedges"] = spearmanr(logit[t], oll[t])[0]
R["spearman_logit_vs_dst_read_error_trueedges"] = spearmanr(logit[t], dst_err[t])[0]
# walk decisions: at each walk node, chosen edge vs other TRUE out-edges of that node
out_by = {}
for i in range(len(E)):
    out_by.setdefault(src[i], []).append(i)
d_err, d_sim, d_len, n = [], [], [], 0
for i in np.where(inw == 1)[0]:
    alts = [j for j in out_by[src[i]] if j != i and lab[j] == 1]
    if not alts: continue
    n += 1
    d_err.append(dst_err[i] - np.mean([dst_err[j] for j in alts]))
    d_sim.append(sim[i] - np.mean([sim[j] for j in alts]))
    d_len.append(oll[i] - np.mean([oll[j] for j in alts]))
d_err, d_sim, d_len = map(np.array, (d_err, d_sim, d_len))
R["walk_decisions_with_true_alternatives"] = n
R["chosen_next_read_cleaner_than_alternatives_frac"] = float((d_err < 0).mean())
R["chosen_next_read_error_minus_alt_median"] = float(np.median(d_err))
R["chosen_higher_similarity_frac"] = float((d_sim > 0).mean())
R["chosen_longer_overlap_frac"] = float((d_len > 0).mean())
R["wilcoxon_p_err"] = float(wilcoxon(d_err)[1])
R = {k: (float(v) if not isinstance(v, int) else v) for k, v in R.items()}
json.dump(R, open(out, "w"), indent=1); print(json.dumps(R, indent=1))
