"""Re-verify the key claim of every task from the files on disk. Prints PASS/FAIL per check; exit 1 on any FAIL."""
import json, csv, hashlib, os, subprocess, glob, pickle, gzip
A = "/home/abhiraj_arya/research/analysis"; G1 = "/home/abhiraj_arya/research/GNNome"; R2 = "/home/abhiraj_arya/research/run2_corrected"
E1 = "/home/abhiraj_arya/research/run1_evidence"
res = []
def check(task, name, ok, detail=""):
    res.append((task, name, bool(ok), detail)); print(f"[{'PASS' if ok else 'FAIL'}] {task:>2} {name}" + (f"  ({detail})" if detail else ""))
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
J = lambda p: json.load(open(p))
tsv = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
seq = lambda p: "".join(l.strip() for l in open(p) if not l.startswith(">")).upper()

# 1 inventory
r = subprocess.run("cd %s && sha256sum -c --quiet hashes/MANIFEST.sha256" % E1, shell=True, capture_output=True)
check(1, "Run 1 evidence manifest verifies", r.returncode == 0)
r = subprocess.run("cd %s && sha256sum -c --quiet MANIFEST.sha256" % R2, shell=True, capture_output=True)
check(1, "Run 2 evidence manifest verifies", r.returncode == 0)
inv = tsv(f"{A}/01_inventory/assemblies.tsv")
bad = [x["label"] for x in inv if sha("/home/abhiraj_arya/research/" + x["assembly"]) != x["sha256"]]
check(1, "all 9 inventoried assemblies still match their recorded SHA256", not bad and len(inv) == 9, ",".join(bad))
# 2 conditions
c = {x["factor"]: x for x in tsv(f"{A}/02_conditions/conditions.tsv")}
check(2, "conditions table has Run 1 and Run 2 columns for every factor", len(c) >= 15 and all(x["run1"] and x["run2"] for x in c.values()))
# 3 reads
r1, r2 = J(f"{A}/03_reads/run1.summary.json"), J(f"{A}/03_reads/run2.summary.json")
check(3, "read counts 14,884 / 14,874", r1["reads"] == 14884 and r2["reads"] == 14874)
check(3, "coverage 60.00x / 30.00x", round(r1["coverage"], 2) == 60.0 and round(r2["coverage"], 2) == 30.0)
check(3, "FASTQ orientation matches MAF (500/500 both)", r1["orientation_check"].startswith("500/500") and r2["orientation_check"].startswith("500/500"))
check(3, "PBSIM MAF '-' strand start shift: 5,224 in Run 1, 0 in Run 2", r1["maf_ref_start_corrected"] == 5224 and r2["maf_ref_start_corrected"] == 0)
check(3, "identical error-type shares (6:55:39 model) in both runs",
      abs(r1["error_share_ins"] - r2["error_share_ins"]) < 0.002 and abs(r1["error_share_mismatch"] - r2["error_share_mismatch"]) < 0.002)
# 4 graphs
g1, g2 = J(f"{A}/05_edge_predictions/run1.summary.json"), J(f"{A}/05_edge_predictions/run2.summary.json")
check(4, "Run 1 graph = 2*S nodes / 2*L edges (5,794 / 74,728)", g1["complexity"]["nodes"] == 5794 and g1["complexity"]["edges"] == 74728)
check(4, "Run 2 graph = 2*S nodes / 2*L edges (1,524 / 5,936)", g2["complexity"]["nodes"] == 1524 and g2["complexity"]["edges"] == 5936)
# 5 edge predictions
t1, t2 = g1["truth"], g2["truth"]
check(5, "edge classes sum to all edges (both runs)", t1["edges_true"] + t1["edges_spurious"] + t1["edges_containment"] == 74728 and
      t2["edges_true"] + t2["edges_spurious"] + t2["edges_containment"] == 5936)
check(5, "GNNome walks use 0 spurious/backward edges (both runs)", not t1["walk_spurious_or_backward_edges"] and not t2["walk_spurious_or_backward_edges"])
check(5, "Run 2 GNNome probabilities collapsed (median p(true) < 1e-6)", t2["prob_true_median"] < 1e-6, f"{t2['prob_true_median']:.2e}")
# 6 sensitivity
S = tsv(f"{A}/06_sensitivity/collected.tsv") if os.path.exists(f"{A}/06_sensitivity/collected.tsv") else []
check(6, "all 21 main replicates collected (7 conditions x 3 seeds), no failures", sum(not x["cond"].startswith("run1reads") for x in S) == 21 and not os.path.exists(f"{A}/06_sensitivity/failures.txt"), f"{sum(not x['cond'].startswith('run1reads') for x in S)} main + {sum(x['cond'].startswith('run1reads') for x in S)} Run 1-subsample")
if S:
    ok = all(float(x["gnnome_source_read_error"]) < float(x["raven_source_read_error"]) for x in S)
    check(6, "GNNome copies more accurate reads than Raven in every replicate", ok)
# 7 repeats
gb = "".join(l[10:].replace(" ", "").strip() for l in open(f"{A}/07_repeats/CP017100.1.gb").read().split("ORIGIN")[1].splitlines() if l.strip() and l.strip() != "//").upper()
check(7, "GenBank sequence identical to CP017100.fasta", gb == seq(f"{G1}/data/CP017100.fasta"))
rr = sorted((int(l.split("\t")[1]), int(l.split("\t")[2])) for l in open(f"{A}/07_repeats/elements.bed") if l.split("\t")[3] == "rRNA")
clusters = 0; last = -10**9
for a, b in rr:
    if a - last > 10000: clusters += 1
    last = b
check(7, "7 rRNA operons annotated (E. coli expectation)", clusters == 7, f"{clusters} clusters")
check(7, "Rhs family has 5 copies", any(x["class"] == "Rhs" and x["copies"] == "5" for x in tsv(f"{A}/07_repeats/repeat_families.tsv")))
# 8 breakpoints
def cnt(lab, src):
    return sorted(x["type"] for x in tsv(f"{A}/08_breakpoints/{lab}.tsv") if x["source"] == src and x["type"] != "origin")
check(8, "old Run 2 GNNome has the K-12-vs-DH5a breakpoint pattern (nucmer)", cnt("K12_reference", "nucmer") == cnt("R2old_GNNome_K12graph", "nucmer"),
      f"{len(cnt('K12_reference', 'nucmer'))} breakpoints")
check(8, "all polished assemblies: 0 breakpoints in every source", all(not [x for x in tsv(f"{A}/08_breakpoints/{l}.tsv") if x["type"] != "origin"]
                                                                       for l in ("R2_Raven_Racon2", "R2_GNNome_Racon2", "R2_Raven_internal_p2")))
# 9 repeat errors
for lab in ("R1_GNNome", "R2_GNNome"):
    check(9, f"{lab}: contig rebuilt from walk is byte-identical to the assembly", J(f"{A}/09_repeat_errors/{lab}.summary.json")["gnnome_walk_rebuild_matches_assembly"] is True)
jc = J(f"{A}/09_repeat_errors/R1_GNNome.junctions.summary.json")
check(9, "Run 1 junction cut-error model r > 0.99", jc["pearson_r"] > 0.99, f"r={jc['pearson_r']:.4f}")
o = J(f"{A}/09_repeat_errors/R1_GNNome.summary.json")["read_origin"]["indel_gt5"]
check(9, "Run 1: all 35 GNNome indels >5 bp are junction-origin", o == {"junction": 35})
txt = open(f"{A}/09_repeat_errors/R2_breakpoint_mechanisms.txt").read()
check(9, "Run 2: 15/15 junction-origin GNNome breakpoints on twin edges", "15/15 junction-origin GNNome breakpoints sit on twin edges" in txt)
check(9, "Raven tRNA-Met collapse traced to edge S1_7704->S1_3959 (cut error -391)", "true cut error -391 <== spans deletion" in open(f"{A}/09_repeat_errors/R1_raven_tRNA_collapse_edge.txt").read())
sp = open(f"{A}/09_repeat_errors/spurious_edges_explained.txt").read()
check(9, "all 4 Run 2 spurious edges join Rhs (F005) copies", sp.split("===== run2")[1].count("F005|Rhs") == 8)
fdr = J(f"{A}/12_statistics/core_statistics.json")["enrichment_fdr"]["significant_q_lt_0p05"]
check(9, "no repeat-class enrichment survives FDR (only homopolymer/low-complexity)", all(x["track"] in ("homopolymer", "lowcomplexity") for x in fdr))
# 10 polishing
a, b = seq(f"{A}/10_polishing/run1/raven/round2.fasta"), seq(f"{A}/10_polishing/run1/gnnome/round2.fasta")
check(10, "Run 1 polished Raven == GNNome, exact substring of reference", a == b and a in seq(f"{G1}/data/CP017100.fasta"))
pr = J(f"{A}/10_polishing/run2/residual_errors.json")
check(10, "Run 2 residual errors: >40 % in homopolymers >=4 (both)", pr["raven"]["frac_indels_in_homopolymer_ge4"] > 0.4 and pr["gnnome"]["frac_indels_in_homopolymer_ge4"] > 0.4)
# 11 figures
figs = [f"F{i}_" for i in range(1, 10)]
have = [f for f in figs if glob.glob(f"{A}/11_figures/{f}*.png") and glob.glob(f"{A}/11_figures/{f}*.pdf")]
check(11, "figures F1-F9 exist as PNG and PDF", len(have) == 9, f"{len(have)}/9")
# 12 statistics
cs = J(f"{A}/12_statistics/core_statistics.json")["gnnome_minus_raven_errors_per_100kb_window"]
check(12, "Run 1 unpolished: GNNome fewer errors in 46/46 windows", cs["run1_unpolished"]["windows_gnnome_lower"] == 46)
check(12, "Run 2 polished: no significant difference (CI spans 0)", cs["run2_racon2"]["ci95"][0] < 0 < cs["run2_racon2"]["ci95"][1])
check(12, "sensitivity statistics written", os.path.exists(f"{A}/12_statistics/sensitivity_statistics.json"))
n_fail = sum(not r[2] for r in res)
print(f"\n{len(res) - n_fail}/{len(res)} checks passed")
json.dump([dict(task=t, check=n, passed=p, detail=d) for t, n, p, d in res], open(f"{A}/12_statistics/verification.json", "w"), indent=1)
raise SystemExit(1 if n_fail else 0)
