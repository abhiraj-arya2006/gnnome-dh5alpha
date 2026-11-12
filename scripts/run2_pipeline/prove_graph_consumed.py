"""Prove which GFA a GNNome run consumed.
Rebuilds node ids and directed edges straight from the GFA (same rules as
graph_parser.only_from_gfa) and compares them with the DGL graph that
inference.py loaded, plus the node sequences GNNome used for the assembly."""
import sys, hashlib, pickle, dgl, torch
gfa, run_dir = sys.argv[1], sys.argv[2]           # run_dir = <data>/raven
rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()

r2n, seqs, exp_edges, S, L = {}, {}, set(), 0, 0
idx = 0
for line in open(gfa):
    f = line.split()
    if f[0] == "S":
        r2n[f[1]] = (idx, idx + 1); seqs[idx] = f[2]; seqs[idx + 1] = rc(f[2]); idx += 2; S += 1
    elif f[0] == "L":
        L += 1
        a, o1, b, o2, cig = f[1:6]
        if int(cig[:-1]) == 0: continue
        A, B = r2n[a], r2n[b]
        sr = A[0] if o1 == "+" else A[1]; dr = B[0] if o2 == "+" else B[1]
        sv = B[1] if o2 == "+" else B[0]; dv = A[1] if o1 == "+" else A[0]
        exp_edges |= {(sr, dr), (sv, dv)}

g = dgl.load_graphs(f"{run_dir}/processed/0.dgl")[0][0]
u, v = g.edges()
got_edges = set(zip(u.tolist(), v.tolist()))
reads = pickle.load(open(f"{run_dir}/info/0_reads.pkl", "rb"))
r2n_pkl = pickle.load(open(f"{run_dir}/info/0_read_to_node.pkl", "rb"))
pred = torch.load(f"{run_dir}/decode/0_predicts.pt")

seq_ok = sum(1 for n, s in seqs.items() if reads.get(n) == s)
print(f"GFA            : {gfa}")
print(f"GFA sha256     : {sha(gfa)}")
print(f"GFA counts     : S={S}  L={L}  -> expected nodes 2*S={2*S}, expected edges 2*L={2*L}")
print(f"GNNome run dir : {run_dir}")
print(f"0.dgl sha256   : {sha(f'{run_dir}/processed/0.dgl')}")
print(f"0.dgl counts   : nodes={g.num_nodes()}  edges={g.num_edges()}")
print(f"0_predicts.pt  : {tuple(pred.shape)} edge scores (one per DGL edge)")
print(f"node ids       : read_to_node.pkl == ids derived from GFA order: {r2n_pkl == r2n}")
print(f"node sequences : {seq_ok}/{len(seqs)} nodes in 0_reads.pkl equal the GFA segment (or its reverse complement)")
print(f"edges          : {len(exp_edges & got_edges)} of {len(exp_edges)} GFA-derived edges present in 0.dgl; "
      f"{len(got_edges - exp_edges)} DGL edges not explained by the GFA")
ok = (g.num_nodes() == 2*S and g.num_edges() == 2*L == len(exp_edges) and got_edges == exp_edges
      and seq_ok == len(seqs) and r2n_pkl == r2n and pred.shape[0] == g.num_edges())
print(f"VERDICT        : {'MATCH - this GNNome run consumed exactly this GFA' if ok else 'NO MATCH'}")
