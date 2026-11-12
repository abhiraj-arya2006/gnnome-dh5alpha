"""Task 10 - what polishing fixes and what survives.

For a pair of polished assemblies of the same run (Raven-derived, GNNome-derived) the residual QUAST events are
characterised by:
  homopolymer context (length of the longest single-base run touching the event, ref coordinates),
  shared vs. unique (same type within 2 bp in the other assembly -> consensus-limited site, not assembler-specific),
  repeats / low complexity / true depth, and
  distance to the genomic positions of the GNNome walk junctions (do twin-edge cut errors survive polishing?).
usage: polish_analysis.py <label> <raven_used_snps> <gnnome_used_snps> <reads.tsv> <gnnome_edges.tsv> <out.json>
"""
import sys, json, bisect
import numpy as np
label, rs, gs, truth, edges, out = sys.argv[1:7]
G = "".join(l.strip() for l in open("/home/abhiraj_arya/research/GNNome/data/CP017100.fasta") if not l.startswith(">")).upper()
N = len(G)
R7 = "/home/abhiraj_arya/research/analysis/07_repeats"
# homopolymer run length at each base
run = np.ones(N, dtype=np.int32); i = 0
while i < N:
    j = i
    while j + 1 < N and G[j + 1] == G[i]: j += 1
    run[i:j + 1] = j - i + 1; i = j + 1
hp_ctx = np.maximum(np.maximum(run, np.roll(run, 1)), np.roll(run, -1))
rep = np.zeros(N, bool)
for l in open(f"{R7}/repeats.bed"): f = l.split("\t"); rep[int(f[1]):int(f[2])] = True
low = np.zeros(N, bool)
for l in open(f"{R7}/lowcomplexity.bed"): f = l.split("\t"); low[int(f[1]):int(f[2])] = True
diff = np.zeros(N + 1)
for k, l in enumerate(open(truth)):
    if k: f = l.split("\t"); diff[int(f[2])] += 1; diff[int(f[3])] -= 1
depth = np.cumsum(diff)[:N]
def events(p):
    ev = []
    for l in open(p):
        f = l.rstrip("\n").split("\t"); rp = int(f[2]) - 1
        t = "ins" if f[3] == "." else ("del" if f[4] == "." else "mismatch")
        ev.append((rp, t, len(f[4]) if t == "ins" else len(f[3])))
    return ev
# GNNome walk junction genomic positions (first base of the destination read)
jg = []
for k, l in enumerate(open(edges)):
    if k == 0: continue
    f = l.rstrip("\n").split("\t")
    if f[18] == "1": jg.append(int(f[9]) if f[8] == "+" else int(f[10]))
jg.sort()
def near(p, lst, d):
    k = bisect.bisect_left(lst, p - d); return k < len(lst) and lst[k] <= p + d
S = dict(label=label, genome_bases_in_homopolymer_ge4=float((hp_ctx >= 4).mean()), genome_repeat_frac=float(rep.mean()))
evs = {"raven": events(rs), "gnnome": events(gs)}
for a, b in (("raven", "gnnome"), ("gnnome", "raven")):
    E = evs[a]; other = sorted((p, t) for p, t, _ in evs[b]); op = [x[0] for x in other]
    def shared(p, t):
        k = bisect.bisect_left(op, p - 2)
        while k < len(op) and op[k] <= p + 2:
            if other[k][1] == t: return True
            k += 1
        return False
    ind = [e for e in E if e[1] != "mismatch"]
    sh = [shared(p, t) for p, t, _ in E]
    S[a] = dict(events=len(E), mismatches=sum(e[1] == "mismatch" for e in E), indels=len(ind),
                insertions=sum(e[1] == "ins" for e in E), deletions=sum(e[1] == "del" for e in E),
                indel_len_1=sum(e[2] == 1 for e in ind), indel_gt5=sum(e[2] > 5 for e in ind),
                indels_in_homopolymer_ge4=sum(hp_ctx[p] >= 4 for p, _, _ in ind),
                frac_indels_in_homopolymer_ge4=float(np.mean([hp_ctx[p] >= 4 for p, _, _ in ind])) if ind else None,
                fold_homopolymer_ge4=(float(np.mean([hp_ctx[p] >= 4 for p, _, _ in ind])) / float((hp_ctx >= 4).mean())) if ind else None,
                in_repeat=int(sum(rep[p] for p, _, _ in E)), in_lowcomplexity=int(sum(low[p] for p, _, _ in E)),
                shared_with_other=int(sum(sh)), unique=int(len(E) - sum(sh)),
                mean_true_depth_at_events=float(np.mean([depth[p] for p, _, _ in E])) if E else None,
                genome_mean_depth=float(depth.mean()),
                within_50bp_of_gnnome_junction=int(sum(near(p, jg, 50) for p, _, _ in E)),
                expected_within_50bp_if_random=len(E) * len(jg) * 101 / N)
json.dump(S, open(out, "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print(json.dumps(S, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
