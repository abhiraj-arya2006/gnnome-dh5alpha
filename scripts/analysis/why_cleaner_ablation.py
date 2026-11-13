"""Causal test: re-score with the pretrained model, (a) real features, (b) overlap_similarity replaced by its mean.
At every node with >=2 TRUE out-edges, does the top-scoring true edge lead to a cleaner read than the other true edges?
Run with cwd = GNNome repo. usage: why_cleaner_ablation.py <0.dgl> <edges.tsv> <reads.tsv> <out.json>"""
import sys, json, csv, numpy as np, torch, dgl
sys.path.insert(0, "."); import models
from configs.hyperparameters import get_hyperparameters
from utils.data_utils import preprocess_graph, add_positional_encoding
dglp, edges_tsv, reads_tsv, out = sys.argv[1:5]
h = get_hyperparameters()
model = models.SymGatedGCNModel(h['node_features'], h['edge_features'], h['dim_latent'], h['hidden_ne_features'], h['num_gnn_layers'],
                                h['hidden_edge_scores'], h['normalization'], dropout=None)
model.load_state_dict(torch.load("weights/weights.pt", map_location="cpu")); model.eval()
g = add_positional_encoding(preprocess_graph(dgl.load_graphs(dglp)[0][0]))
def score(e):
    pi = g.ndata['in_deg'].unsqueeze(1).float(); pi = (pi - pi.mean()) / pi.std()
    po = g.ndata['out_deg'].unsqueeze(1).float(); po = (po - po.mean()) / po.std()
    with torch.no_grad(): return model(g, torch.cat((pi, po), 1), e).squeeze(-1).numpy()
e_real = g.edata['e'].clone(); s_real = score(e_real)
e_abl = e_real.clone(); e_abl[:, 1] = e_real[:, 1].mean(); s_abl = score(e_abl)
E = list(csv.DictReader(open(edges_tsv), delimiter="\t"))
logit_tsv = np.array([float(x["logit"]) for x in E])
err = {r["read"]: 1 - float(r["identity"]) for r in csv.DictReader(open(reads_tsv), delimiter="\t")}
lab = np.array([int(x["label"]) for x in E]); derr = np.array([err[x["dst_read"]] for x in E])
src = g.edges()[0].numpy(); by = {}
for i in np.where(lab == 1)[0]: by.setdefault(src[i], []).append(i)
def test(s):
    win, d = 0, []
    nodes = [v for v in by.values() if len(v) >= 2]
    for v in nodes:
        top = max(v, key=lambda j: s[j]); others = [derr[j] for j in v if j != top]
        d.append(derr[top] - np.mean(others)); win += derr[top] < min(others)
    d = np.array(d)
    return dict(nodes=len(nodes), top_edge_cleaner_than_mean_alt_frac=float((d < 0).mean()),
                top_edge_cleanest_of_all_frac=float(win / len(nodes)), median_err_diff=float(np.median(d)),
                chance_cleanest=float(np.mean([1 / len(v) for v in nodes])))
R = dict(model_reproduces_saved_logits_maxabsdiff=float(np.abs(s_real - logit_tsv).max()),
         real=test(s_real), similarity_ablated=test(s_abl),
         rank_corr_scores_real_vs_ablated=float(np.corrcoef(np.argsort(np.argsort(s_real)), np.argsort(np.argsort(s_abl)))[0, 1]))
json.dump(R, open(out, "w"), indent=1); print(json.dumps(R, indent=1))
