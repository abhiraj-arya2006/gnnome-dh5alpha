# GNNome vs Raven on *E. coli* DH5α — Run 1 & Run 2 Evidence Report

**Prepared:** 2026-11-12 · **Machine:** WSL2 (Linux 6.18), 4 CPUs, 11 GB RAM · **Reference:** `CP017100.1` *E. coli* K-12 NEB 5-alpha (DH5α), 4,583,637 bp, SHA256 `5e8be5a50a977206b334dfb2ec1489fbb65d7a88ee77aa63bb72b1ff15945d75`

All paths are under `/home/abhiraj_arya/research/`. `R2/` = `run2_corrected/`, `R1E/` = `run1_evidence/`, `G/` = `GNNome/`.

---

## 0. Direct answers to the five evidence requests

### Q1 — Run 2 read hash, read count and PBSIM3 seed

| Item | Value | Evidence |
|---|---|---|
| Reads file used in Run 2 | `G/data/dh5alpha/simulated_reads.fastq.gz` (named in the original Run 2 command, `~/.bash_history` line 160); read-only copy at `R2/00_reads/simulated_reads.fastq.gz` | `R2/00_reads/history_evidence.txt` |
| File SHA256 | `f5f25cc4da2201e12d4dc74024c468fe5b9b6c49447e16d9193251807815e38e` | `R2/00_reads/reads.sha256` |
| Content SHA256 (decompressed FASTQ) | `5a6ba8d273b356fa8f3bfb9f2406c993199c1f85dd2a5a94b08f5360826b69dc` — identical for `simulated_0001.fq.gz` and `dh5alpha_reads.fastq.gz` (those differ only in gzip header) | `R2/00_reads/reads.sha256` |
| Read count | **14,874** (also `loaded 14874 sequences` in the Raven log) | `R2/00_reads/read_stats.txt`, `R2/02_raven_p0/raven.log` |
| Total bases | **137,509,430 bp** | `R2/00_reads/read_stats.txt` |
| Coverage | **30.00×** (137,509,430 / 4,583,637) | `R2/00_reads/read_stats.txt` |
| Read length | mean 9,245 bp, N50 12,942 bp, min 109, max 72,272 | `R2/00_reads/read_stats.txt` |
| Read accuracy | PBSIM's own statistic: mean 0.8513 (SD 0.0416); alignment identity from the MAF (matching columns / all columns): **86.01 %** | `R2/08_proofs/pbsim3_regeneration.log`, `R2/00_reads/simulation_metadata.txt` |
| **PBSIM3 seed** | **`1783636376`** (= Unix time 2026-07-09 22:32:56 UTC, the same second `simulated_0001.ref` was written) | `R2/08_proofs/pbsim3_seed_and_command.txt` |
| **PBSIM3 command (recovered)** | `./vendor/pbsim3/src/pbsim --strategy wgs --method qshmm --qshmm ./vendor/pbsim3/data/QSHMM-RSII.model --depth 30 --genome data/dh5alpha/CP017100.fasta --prefix data/dh5alpha/simulated --seed 1783636376` | same |
| PBSIM3 build | `G/vendor/pbsim3/src/pbsim`, SHA256 `5ce6a5e9446528e19a0cdab397b0832369a32030bbb8050d2fe24214196c79d5`, source commit `35ce0debb8b9df8a45ac9bb85420686776e9e721` | same |

**How the seed was proven.** The original PBSIM3 command was not logged. PBSIM3 uses `seed = time(NULL)` unless `--seed` is given (`pbsim.cpp:253`) and prints it only to stderr, which was not saved. Every Unix second from 22:30:36 to 22:34:37 UTC was tried as a seed (`R2/scripts/find_pbsim_seed.sh`); only `1783636376` reproduces the first read. A **full 30× regeneration** with that seed is **byte-identical** to the original:

| File | Regenerated content SHA256 | Original content SHA256 |
|---|---|---|
| `simulated_0001.fq.gz` | `5a6ba8d2…826b69dc` | `5a6ba8d2…826b69dc` ✅ |
| `simulated_0001.maf.gz` | `48653c61…f0e437` | `48653c61…f0e437` ✅ |
| `simulated_0001.ref` | byte-identical | ✅ |

So the seed and every content-affecting PBSIM3 parameter are proven. Only `--prefix` is inferred from the file names; it does not affect content.

### Q2 — Exact Raven graph-output command and the `print_graphs` commit

| Item | Value |
|---|---|
| Binary | `G/vendor/raven-1.8.1/build/bin/raven`, Raven 1.8.1, SHA256 `2ba8260b04c6f19ec77da5e33d1b2cf6e63f0ff76e162bdbc20a7f6cc3f3cbdd`, built 2026-07-08 16:14:15 UTC, Release, `racon_enable_cuda=OFF` |
| Source commit | **`3207c8d9819174cd28d0a250a9c9a3156811e890`** — *"Print graphs in CSV and GFA"* (2022-05-03), remote `https://github.com/lbcb-sci/raven.git`, 0 uncommitted changes, no source file newer than the binary |
| Graph-output mechanism | `RavenLib/src/assemble.cc:871–872`: `PrintCsv(graph, "graph_1.csv", …)` and `PrintGfa(graph, "graph_1.gfa")` at the start of `Assemble()`. **There is no flag.** The raw overlap graph (before transitive reduction and polishing) is always written to the current working directory. |
| **Exact command** | `cd /home/abhiraj_arya/research/run2_corrected/02_raven_p0 && /usr/bin/time -v /home/abhiraj_arya/research/GNNome/vendor/raven-1.8.1/build/bin/raven -t 4 -p0 /home/abhiraj_arya/research/run2_corrected/00_reads/simulated_reads.fastq.gz > raven_p0_assembly.fasta 2> raven.log` |
| Result | `exit=0`, 23.7 s, 0.6 GB RAM; `graph_1.gfa` written 2026-11-12 18:01:42 UTC |
| Output graph | `R2/02_raven_p0/graph_1.gfa`, **SHA256 `6d5df9820ffddae793538c091eb8ac98013a890ab8ed1bd402b454ed4c8a5483`, S = 762, L = 2,968** |
| Conda Raven ruled out | Behavioural test on 3,000 reads: the fork wrote `graph_1.gfa` + `graph_1.csv`; Conda Raven 1.8.1 (`c0e9d370…2b22`) wrote **no** graph (`R2/01_raven_build/behavior_test.txt`). Conda's `--graphical-fragment-assembly` writes only the collapsed final graph (`S=1`), unusable by GNNome. |
| Determinism | A second, separate invocation (`-p 2`) wrote a graph with the **same SHA256** `6d5df982…5483` (`R2/07_polish/determinism_check.txt`). |

Evidence: `R2/08_proofs/raven_graph_command.txt`, `R2/02_raven_p0/command.txt`, `R2/02_raven_p0/raven.log`, `R2/01_raven_build/`.

### Q3 — Hash of the checkpoint used

| Item | Value |
|---|---|
| Checkpoint | `G/weights/weights.pt` (original pretrained weights; 940,485 B; mtime 2026-07-07 09:21:07 UTC) |
| **SHA256** | **`240076f02c737f6d37fa07b6830f587418b0698496d637337d1378fa3f6740ff`** |
| Run 2 | passed explicitly (`--model weights/weights.pt`); log: *"Loading model parameters from: weights/weights.pt" → "Computing the scores with the model…" → "Model used: weights/weights.pt"* (`R2/04_gnnome/inference.log`, `R2/04_gnnome/weights.sha256`) |
| Run 1 | same file; a sandboxed re-inference with it reproduces the Run 1 assembly and walks **byte-for-byte** (`R1E/checks/reinference.sha256`) |
| Training | **None.** No training command was run in either run. |

### Q4 — Proof that `inference.py` consumed the new graph, not an old one

`inference.py` loads the DGL graph `<data>/raven/processed/0.dgl` and the node data in `<data>/raven/info/`, both produced from the GFA by `create_inference_graphs.py`. The proof script `R2/scripts/prove_graph_consumed.py` rebuilds node IDs and every directed edge **directly from the GFA** (same rules as `graph_parser.only_from_gfa`: 2 nodes per `S` line, 2 edges per `L` line). It then compares them with the graph and node sequences GNNome actually used.

| | **Corrected Run 2** | Run 1 | Old Run 2 (invalid) |
|---|---|---|---|
| GFA | `R2/02_raven_p0/graph_1.gfa` | `G/data/graph_1.gfa` | `G/graph_1.gfa` (K-12) |
| GFA SHA256 | `6d5df982…8a5483` | `b1209948…100270` | `62176d9d…dcdfc` |
| GFA S / L | **762 / 2,968** | 2,897 / 37,364 | 2,561 / 27,451 |
| Expected nodes / edges (2S / 2L) | **1,524 / 5,936** | 5,794 / 74,728 | 5,122 / 54,902 |
| `0.dgl` nodes / edges | **1,524 / 5,936** ✅ | 5,794 / 74,728 ✅ | 5,122 / 54,902 |
| `0.dgl` SHA256 | `fa0859d6…fbde` | `247ec2cd…0ce4` | `569148b8…c1ff9` |
| `0_predicts.pt` (scores per edge) | **5,936** ✅ | 74,728 ✅ | 54,902 |
| Node IDs equal GFA order | **True** | True | True |
| Node sequences equal GFA segments | **1,524 / 1,524** | 5,794 / 5,794 | 5,122 / 5,122 |
| GFA edges found in `0.dgl` | **5,936 / 5,936; 0 unexplained** | 74,728 / 74,728; 0 unexplained | 54,902 / 54,902; 0 unexplained |
| `inference.log` | `nb_original_nodes: 1524` | — | — |
| **Verdict** | **Consumed exactly the new DH5α graph** | Consumed the Run 1 DH5α graph | Consumed the **K-12** graph |

**Negative controls:** the corrected Run 2 GNNome graph vs the K-12 GFA has 0/5,122 matching sequences, 8/54,902 edges and 5,928 unexplained edges → **NO MATCH**. The old Run 2 graph vs the new DH5α GFA has 0/1,524 sequences → **NO MATCH**.

Order of events: `graph_1.gfa` written 18:01:42 UTC → preprocessing 18:02 → inference 18:03 (2026-11-12). The fresh output directory `R2/04_gnnome/` had no cached `decode/0_predicts.pt`, so scores were computed by the model, not reloaded.

Evidence: `R2/08_proofs/graph_consumed_run2.txt`, `R2/08_proofs/graph_consumed_run1_and_old_run2.txt`, `R2/04_gnnome/{commands.txt,preprocess.log,inference.log,pkl_vs_gfa.txt}`.

### Q5 — Full QUAST reports and command lines, same QUAST version

All evaluations used **QUAST 5.3.0** (`/home/abhiraj_arya/miniconda3/envs/gnnome/bin/quast.py`), the **same reference file** `G/data/CP017100.fasta` (4,583,637 bp), and **default settings**: `MODE: default, min contig length 500, min alignment length 65, min alignment IDY 95.0, ambiguity one, min local misassembly 200, min extensive misassembly 1000`. Each log has 0 errors and 1 warning: *"Python locale settings can't be changed"* (an environment notice).

**The three Run 2 assemblies** (full `report.txt` in Appendix A; complete report folders with `report.html/.tsv/.pdf`, `contigs_reports/` and `icarus.html` at the paths shown):

| Assembly | Command line (verbatim from `quast.log`, run from `G/`) | Report folder |
|---|---|---|
| Raven + Racon ×2 | `quast.py /home/abhiraj_arya/research/run2_corrected/07_polish/raven_racon2/round2.fasta -r data/CP017100.fasta -o /home/abhiraj_arya/research/run2_corrected/07_polish/quast/raven_racon2` | `R2/07_polish/quast/raven_racon2/` |
| GNNome + Racon ×2 | `quast.py /home/abhiraj_arya/research/run2_corrected/07_polish/gnnome_racon2/round2.fasta -r data/CP017100.fasta -o /home/abhiraj_arya/research/run2_corrected/07_polish/quast/gnnome_racon2` | `R2/07_polish/quast/gnnome_racon2/` |
| Raven `-p 2` (built-in polishing) | `quast.py /home/abhiraj_arya/research/run2_corrected/07_polish/raven_internal_p2/raven_p2_assembly.fasta -r data/CP017100.fasta -o /home/abhiraj_arya/research/run2_corrected/07_polish/quast/raven_internal_p2` | `R2/07_polish/quast/raven_internal_p2/` |

**Also recorded** (same version, same reference):

| Assembly | Command line | Note |
|---|---|---|
| Run 2 Raven `-p0` (unpolished) | `quast.py …/02_raven_p0/raven_p0_assembly.fasta -r data/CP017100.fasta -o …/05_quast/raven_p0` | default settings; alignment failed (identity < 95 %) — basic stats only |
| Run 2 GNNome (unpolished) | `quast.py …/04_gnnome/raven/assembly/0_assembly.fasta -r data/CP017100.fasta -o …/05_quast/gnnome_p0` | as above |
| Run 2 Raven `-p0`, GNNome | `quast.py --min-identity 80 … -o …/05_quast/minid80_{raven_p0,gnnome_p0}` | **supplementary, non-default** |
| Run 1 Raven `-p0` | `quast.py data/CP017100_assembly.fasta -r data/CP017100.fasta -o data/classical_raven_quast_results` | original, 2026-07-09 |
| Run 1 GNNome | `quast.py data/raven/assembly/0_assembly.fasta -r data/CP017100.fasta -o data/raven/quast_results` | original, 2026-07-09 |

**Thread count.** QUAST picks its thread count automatically: the Run 1 runs (July) used 3 threads and the Run 2 runs used 1. Two controls show this does not change any metric:
1. The July QUAST of the old Run 2 Raven assembly (3 threads) and today's QUAST of the byte-identical Raven `-p 2` assembly (1 thread) give identical `report.tsv` values.
2. Re-scoring Raven + Racon ×2 with `-t 3` gives the same values as the 1-thread run (`R2/08_proofs/quast_threads3_rerun_raven_racon2/`).

Evidence: `R2/08_proofs/quast_settings_proof.txt` (command, version, start time, reference and full parameter line for every QUAST run).

---

## 1. Executive summary

| Question | Answer |
|---|---|
| Is Run 1 valid? | **Yes — fully verified.** Genuine DH5α graph, correct reads, original checkpoint, GitHub copies identical to local files. |
| Is the *old* Run 2 GNNome result valid? | **No.** It consumed the stale **K-12** `graph_1.gfa` (proven edge-by-edge, Q4). It must not be reported. |
| Is the *old* Run 2 Raven baseline valid? | **Yes.** Our fresh Raven `-p 2` run reproduced it **byte-for-byte**. |
| Is the *corrected* Run 2 valid? | **Yes.** Proven reads and seed, verified Raven fork, new DH5α graph consumed exactly by GNNome, original checkpoint, full logs. |
| Was a new GNNome inference needed for Run 1? | **No.** A sandboxed re-inference reproduced the Run 1 assembly byte-for-byte. |
| Was any training done? | **No.** Only the pretrained `weights/weights.pt` was used. |

### Who performed better

| Setting | Winner | Margin |
|---|---|---|
| **Run 1** — 60× high-accuracy reads, unpolished (`-p0`) | **GNNome** | **44 % fewer mismatches** (394 vs 708; 8.60 vs 15.45 /100 kbp) and **41 % fewer indels** (4,327 vs 7,375; 94.45 vs 160.98 /100 kbp); otherwise identical (1 contig, 99.930 % genome fraction, 0 misassemblies) |
| **Run 2** — 30× noisy reads (86 %), unpolished (`-p0`) | **GNNome (slightly)**, supplementary only | Identity to reference 88.26 % vs 88.02 %; unaligned 3,827 vs 23,614 bp at QUAST `--min-identity 80` (non-default setting) |
| **Run 2** — after identical external Racon ×2 | **Tie; Raven marginally ahead** | Mismatches 34 vs 35, indels 1,312 vs 1,368 (0.74 vs 0.76 and 28.66 vs 29.88 /100 kbp) across 4.58 Mb; same genome fraction (99.867 %), 0 misassemblies for both |
| Polishing method (Run 2) | **External minimap2 + Racon ×2** beats Raven's built-in `-p 2` | Raven `-p 2`: 50 mismatches, 1,614 indels (1.09 / 35.25 per 100 kbp) |

**Bottom line:** the pretrained GNNome model gives a clear accuracy gain on the accurate-read graph (Run 1). On noisy reads (Run 2), which are unlike its HiFi training data, its advantage is small before polishing and disappears after polishing, because Racon consensus corrects nearly all base errors on both paths.

> Run 1 and Run 2 used **different read sets** (60× near-perfect vs 30× noisy). Compare GNNome vs Raven **within** each run, never numbers across runs.

---

## 2. Evidence locations

| Folder | Contents |
|---|---|
| `R1E/github_d3ea188/results/ecoli-dh5a/full/` | 99 files extracted from `development-team2` commit `d3ea188` (read-only) |
| `R1E/local_run1/` | Read-only copy of local Run 1 files with original July 9 timestamps (112 files) |
| `R1E/{checks,hashes,scripts,reinference}/` | Run 1 check outputs, SHA256 lists, checker script, sandboxed re-inference |
| `R1E/hashes/MANIFEST.sha256` | 211-file manifest (verified OK) |
| `R2/00_reads/` | Verified Run 2 reads + PBSIM metadata (read-only copies) |
| `R2/01_raven_build/` | Raven build verification + behavioural test |
| `R2/02_raven_p0/` | Fresh graph, unpolished Raven assembly, command, full log |
| `R2/03_graph_validation/` | Graph hashes, counts, read and reference checks |
| `R2/04_gnnome/` | GNNome preprocessing + inference inputs, logs, outputs |
| `R2/05_quast/` | QUAST on unpolished assemblies (default + supplementary `--min-identity 80`) |
| `R2/06_compare/` | Comparison tables, identity tables, strain controls |
| `R2/07_polish/` | Racon ×2 polishing (both assemblies), Raven `-p 2`, QUAST of polished assemblies |
| `R2/08_proofs/` | PBSIM3 seed + regeneration, Raven command, graph-consumption proofs, QUAST settings proof |
| `R2/scripts/` | `check_gfa_reads.py`, `compare_quast.py`, `polish2.sh`, `find_pbsim_seed.sh`, `prove_graph_consumed.py` |
| `R2/MANIFEST.sha256` | 472-file manifest (verified OK) |

Original files in `GNNome/` and `development-team2/` were **not modified**: July timestamps intact, and `git status` shows no changes under `data/` or `weights/`. All seed searches, the PBSIM regeneration and the re-inference wrote only to temporary or new evidence folders.

---

## 3. Run 1 — recovery and verification

### 3.1 Recovery from GitHub
- Commit **`d3ea1889f6341390efdb282c169a1843c1a06a56`** — `development-team2`, 2026-07-10 15:48:35 +0530, author `abhiraj-arya2006`, *"results: add remaining K12 example files and DH5α full baseline results"*.
- On GitHub: reachable from `origin/main`, `origin/abhiraj/baseline-results`, `origin/akshat/fedrann-reproduction`.
- Extracted with `git archive` (working tree untouched).

### 3.2 Inputs

| Item | Path | Value |
|---|---|---|
| Reads (60× PBSIM3) | `G/data/CP017100_sim.fastq` | **14,884 reads, 275,035,078 bp, 60.00×**; SHA256 `5d3b823e71fb4f4882e838d4dea4bdc57e1c1d6026f8b9f4b31a2418e358ac00` |
| Reads (gz original) | `G/data/CP017100_sim_0001.fq.gz` | SHA256 `dba9f62ecce42433760cc92007be11e9df056eefca7121479d1860581aa56ad1` |
| Reference | `G/data/CP017100.fasta` | SHA256 `5e8be5a5…5d75` (identical to `data/dh5alpha/CP017100.fasta` and both PBSIM `.ref` files) |
| Checkpoint | `G/weights/weights.pt` | SHA256 `240076f0…40ff` |

### 3.3 Graph provenance

| Graph | SHA256 | S / L | Verdict |
|---|---|---|---|
| **Run 1** `G/data/graph_1.gfa` (2026-07-09 21:05:26 UTC, 61,909,148 B) | `b12099486f6d317024e7794759a3732c5321fe73e5236c51771b9a7027100270` | 2,897 / 37,364 | **DH5α** |
| GitHub `d3ea188` copy | `b12099486f6d…100270` | — | **identical** to local |
| K-12 `G/graph_1.gfa` (2026-07-09 22:15:34 UTC) | `62176d9d91749b6f8fdd5794b05da1a5070ec7e5e23b9963018ab7ee642dcdfc` | 2,561 / 27,451 | K-12 |

Sequence-level read check (200 random segments, seed 42; read names prove nothing because every dataset uses `S1_N` names):

| Graph | vs 60× DH5α reads | vs K-12 reads | vs Run 2 30× reads |
|---|---|---|---|
| Run 1 `data/graph_1.gfa` | **200/200** | 0/200 | 0/200 |
| K-12 `graph_1.gfa` | 0/200 | **200/200** | 0/200 |

- Raven header `>Utg5800 LN:i:4581274 RC:i:2897 XO:i:0`: `RC:i:2897` matches the 2,897 graph segments.
- Only the `print_graphs` fork writes `graph_1.gfa` (Q2), so the Run 1 graph came from the fork.

### 3.4 GNNome inference
- Edge-level proof (Q4): Run 1 `0.dgl` = **5,794 nodes / 74,728 edges = 2 × (2,897 / 37,364)**; all nodes, sequences and edges match `data/graph_1.gfa`.
- **Sandboxed re-inference** (copy of `processed/` + `info/` only, `--model weights/weights.pt`): `0_assembly.fasta` and `0_walks.pkl` **byte-identical** to the originals; edge scores differ by at most 0.0011 (CPU floating-point noise).

### 3.5 Assemblies

| Assembly | Path | Length | SHA256 | = GitHub copy |
|---|---|---|---|---|
| GNNome | `G/data/raven/assembly/0_assembly.fasta` | 4,581,340 bp | `70b1e14862e14f1ca22bccb0722b8755541d2a59cb657bbeb81c6d7422e0cb37` | ✅ |
| Raven `-p0` | `G/data/CP017100_assembly.fasta` | 4,581,274 bp | `07fd962adc10701349609caed01b8ca0841b451034fc6bb79c78db2c06844b43` | ✅ |

### 3.6 QUAST (original reports, preserved)
QUAST 5.3.0, default settings (3 threads), same reference. Local reports are identical to the GitHub copies (`diff -rq` empty).
- GNNome: `G/data/raven/quast_results/` (2026-07-09 21:06:20 UTC)
- Raven: `G/data/classical_raven_quast_results/` (2026-07-09 21:12:11 UTC)

| Metric | Raven `-p0` | GNNome |
|---|---:|---:|
| # contigs | 1 | 1 |
| Total length | 4,581,274 | 4,581,340 |
| Genome fraction (%) | 99.930 | 99.930 |
| Duplication ratio | 1.000 | 1.000 |
| # misassemblies | 0 | 0 |
| Unaligned length | 0 | 0 |
| NGA50 | 4,581,274 | 4,581,340 |
| # mismatches (per 100 kbp) | 708 (15.45) | **394 (8.60)** |
| # indels (per 100 kbp) | 7,375 (160.98) | **4,327 (94.45)** |
| Indels ≤5 bp / >5 bp | 7,372 / 3 | 4,292 / 35 |
| Total indel length | 8,315 | 5,003 |

### 3.7 Run 1 — verified vs uncertain
**Verified:** reads, graph, graph consumed by GNNome (edge-level), checkpoint (byte-identical re-inference), both assemblies, QUAST version/settings/reference, GitHub ↔ local identity.

**Uncertain:**
- The original Raven command and log were never saved. `-p0` is supported by the fork graph output and the unpolished error rate; `-t 8` cannot be verified (thread count does not change the graph, see Q2 determinism).
- No original GNNome command or log was saved; `--data data --asm raven --out data/raven` is inferred from output paths.
- The Run 1 PBSIM3 command and seed were not recovered (not requested; same method as Q1 could be applied).
- `production.md` / `reproduction_proof.md` name the wrong read file (`dh5alpha_reads.fastq.gz`, created 89 min *after* the assembly). Real input: `data/CP017100_sim.fastq`.

---

## 4. Run 2 — corrected re-run

### 4.1 Reads, Raven build and fresh graph
See **Q1** (reads, seed, PBSIM3 command) and **Q2** (Raven build, exact command, graph SHA256 and S/L counts).

Raven `-p0` assembly: `>Utg1530 LN:i:4685099 RC:i:760`, SHA256 `b616d36b437960b013cc3ab775d56d54480e158d9e4e7725e8e328e2fe52c4d6`.

### 4.2 Graph validation

| Check | Result |
|---|---|
| Differs from K-12 graph (`62176d9d…`) and Run 1 graph (`b1209948…`) | ✅ |
| Segments vs Run 2 reads (200 random) | **200/200** identical sequences |
| Segments vs Run 1 reads / K-12 reads | 0/200 / 0/200 |
| Same graph from a second Raven invocation | ✅ identical SHA256 |

### 4.3 GNNome inference (no training)
- Preprocessing: `python create_inference_graphs.py --reads /home/abhiraj_arya/research/run2_corrected/00_reads/simulated_reads.fastq.gz --gfa /home/abhiraj_arya/research/run2_corrected/02_raven_p0/graph_1.gfa --asm raven --out /home/abhiraj_arya/research/run2_corrected/04_gnnome` — exit 0.
- Inference: `python inference.py --data /home/abhiraj_arya/research/run2_corrected/04_gnnome --asm raven --out /home/abhiraj_arya/research/run2_corrected/04_gnnome/raven --model weights/weights.pt` — exit 0.
- Graph consumed: proven edge-by-edge (Q4). Checkpoint: Q3.
- Outputs: `0.dgl` `fa0859d6ecc3b55dff2a969628c84999774163f60891662aa0208f03ee39fbde`, `0_predicts.pt` `c2e90422fb1c8ced73e9d66896a992b610a94e04eafd87aae7900bd39f6dd398`, **`0_assembly.fasta` 4,680,964 bp, SHA256 `fa0498763b2ab7885fc943ee94bb104cad41b8fb3c758530835eed40c70b81a4`**.

### 4.4 Unpolished evaluation — why default QUAST could not score it
Default QUAST reported both unpolished Run 2 assemblies as 100 % unaligned (*"Failed aligning the contigs for all the assemblies. Only basic stats are going to be evaluated."*). QUAST's aligner found the full-length alignment (`de:f:0.0997`, about 10 % divergence), but QUAST discards alignments below **95 % identity**. Raven `-p0` has no consensus step, so an assembly built from 86 %-accurate reads stays about 88 % accurate:

| Assembly | Identity vs CP017100 | Identity vs K-12 | Ref covered (CP017100) |
|---|---:|---:|---:|
| Run 2 Raven `-p0` | 88.02 % | 84.82 % | 4,472,255 bp |
| Run 2 GNNome | **88.26 %** | 86.80 % | **4,541,820 bp** |
| Run 1 Raven `-p0` (contrast) | 99.80 % | 96.07 % | 4,580,439 bp |

(minimap2 2.28 `-cx asm20`; `R2/06_compare/identity_unpolished.txt`.) Both are closer to DH5α than to K-12. The ~100 kb excess length comes from insertion errors in the reads.

**Supplementary** (non-default `--min-identity 80`; not comparable with Run 1):

| Metric | Raven `-p0` | GNNome |
|---|---:|---:|
| Genome fraction (%) | 99.375 | **99.786** |
| Local misassemblies | 15 | 19 |
| Unaligned length | 23,614 | **3,827** |
| NGA50 | 4,661,485 | **4,677,137** |
| Mismatches /100 kbp | 2,481.05 | **2,278.59** |
| Indels /100 kbp | 7,741.78 | **7,437.99** |

### 4.5 Polishing
Tools in a separate env (`polish`): **Racon 1.5.0**, **minimap2 2.28-r1209**.

- **Fair comparison:** identical 2 rounds on each unpolished assembly — `minimap2 -t 4 -x map-pb <draft> <reads>` → `racon -t 4 -m 3 -x -5 -g -4 -w 500` (Raven's own scoring values). Script: `R2/scripts/polish2.sh`; commands and logs in each output folder.
- **Reference run:** the same Raven fork with built-in polishing, `raven -t 4 -p 2 <reads>` — exit 0, 1 min 35 s (`R2/07_polish/raven_internal_p2/`).

| Assembly | Length | Identity vs CP017100 | SHA256 |
|---|---:|---:|---|
| Raven + Racon ×2 | 4,578,562 | 99.97 % | `50b2c09eaa096932cd671db27744b6fa2c58efb0cf72739f887c7a9a5c14796c` |
| GNNome + Racon ×2 | 4,578,651 | 99.97 % | `62ad524b1ea557825b80764871c86e05569a6dc4d409917cbb02110a4d1c0f96` |
| Raven `-p 2` (built-in) | 4,578,923 | 99.96 % | `7e6cbec861d420a8d5cac8a545c6b130a42e4c2c9344150643ae419f4e9cece0` |

Identity after round 1 → round 2: Raven 99.88 → 99.97 %, GNNome 99.84 → 99.97 %.

**Reproduction of the old Run 2 baseline:** the Raven `-p 2` assembly is **byte-identical** to the old Run 2 Conda output `G/data/dh5alpha/raven/output/dh5alpha_asm.fasta` (`7e6cbec8…`), and all its QUAST metrics match the old `G/data/dh5alpha/classical_raven_quast_results/`. The old Run 2 Raven baseline is therefore valid and reproducible.

### 4.6 Polished QUAST results (QUAST 5.3.0, default settings, `CP017100.fasta`)

| Metric | Raven + Racon ×2 | GNNome + Racon ×2 | Raven `-p 2` |
|---|---:|---:|---:|
| # contigs | 1 | 1 | 1 |
| Total length | 4,578,562 | 4,578,651 | 4,578,923 |
| Genome fraction (%) | 99.867 | 99.867 | 99.867 |
| Duplication ratio | 1.000 | 1.000 | 1.000 |
| # misassemblies / local | 0 / 0 | 0 / 0 | 0 / 0 |
| Unaligned length | 0 | 0 | 0 |
| NGA50 | 4,578,534 | 4,578,623 | 4,578,890 |
| # mismatches (per 100 kbp) | **34 (0.74)** | 35 (0.76) | 50 (1.09) |
| # indels (per 100 kbp) | **1,312 (28.66)** | 1,368 (29.88) | 1,614 (35.25) |
| Indels ≤5 bp / >5 bp | 1,312 / 0 | 1,367 / 1 | 1,614 / 0 |
| Total indel length | 1,488 | 1,543 | 1,827 |

### 4.7 Strain control (polished assemblies vs K-12 `NC_000913.3`, QUAST `--fast`)

| | Raven + Racon ×2 | GNNome + Racon ×2 | Raven `-p 2` |
|---|---:|---:|---:|
| Misassemblies vs K-12 | 26 | 26 | 26 |
| Unaligned vs K-12 (bp) | 38,614 | 38,614 | 38,619 |
| Genome fraction vs K-12 (%) | 97.622 | 97.606 | 97.589 |

Against CP017100 the same assemblies have **0** misassemblies and **0** unaligned bp → all are **DH5α**, not K-12.

---

## 5. Results that must NOT be reported

| Result | Location | Reason |
|---|---|---|
| Old Run 2 GNNome assembly + QUAST | `G/data/dh5alpha/raven/assembly/0_assembly.fasta`, `G/data/dh5alpha/quast_gnnome_raven/`, `development-team2/results/ecoli-dh5a/full/quast_report/` (after commit `a5f329e`) | Consumed the **K-12** graph (proven edge-by-edge, Q4). Its 23 "misassemblies" vs DH5α are K-12/DH5α strain differences, not GNNome errors. |
| Unpolished Run 2 default QUAST | `R2/05_quast/raven_p0`, `R2/05_quast/gnnome_p0` | Kept as evidence only; alignment failed (identity < 95 %), so there are no accuracy metrics. |
| Failed earlier attempt | `G/data/dh5alpha_corrected/raven_graph/` | Wrong relative path; empty `assembly.fasta`. |

---

## 6. Caveats

1. **Read-type mismatch.** The pretrained GNNome model was trained on accurate (HiFi-like) data. Run 2 reads are 86 %-accurate (PBSIM3 QSHMM-RSII model), so Run 2 tests GNNome outside its training regime.
2. **Raven thread count.** Run 1 used `-t 8` (per docs); Run 2 used `-t 4` (4-CPU machine). Two separate Run 2 Raven invocations wrote byte-identical graphs.
3. **QUAST thread count.** 3 threads (Run 1) vs 1 (Run 2); proven not to change any metric (Q5).
4. **Single replicate.** Each run uses one simulated read set. The Run 2 polished differences (1 mismatch, 56 indels) are small enough that a different simulation seed could change their direction.
5. **External vs built-in polishing.** Raven `-p 2` aligns reads with its own internal mapper; the fair comparison uses minimap2 + Racon for both. Report Raven `-p 2` as a reference point only.

---

## 7. Recommended documentation fixes
- In `production.md` / `reproduction_proof.md`: replace `data/dh5alpha/dh5alpha_reads.fastq.gz` with the real Run 1 input `data/CP017100_sim.fastq`, and mark the Run 1 Raven and GNNome commands as **reconstructed, not logged**.
- Record the recovered Run 2 PBSIM3 command and seed (Q1).
- In `development-team2/results/ecoli-dh5a/full/`: replace the K-12-contaminated GNNome result (commit `a5f329e`) with the corrected Run 2 outputs, and restore or link the Run 1 results from `d3ea188`.

---

## 8. Verification
```bash
cd ~/research/run1_evidence && sha256sum -c --quiet hashes/MANIFEST.sha256 && echo RUN1_OK
```
```bash
cd ~/research/run2_corrected && sha256sum -c --quiet MANIFEST.sha256 && echo RUN2_OK
```
Re-run the graph-consumption proof (in the `gnnome` environment):
```bash
cd ~/research/GNNome && python ~/research/run2_corrected/scripts/prove_graph_consumed.py ~/research/run2_corrected/02_raven_p0/graph_1.gfa ~/research/run2_corrected/04_gnnome/raven
```

---

## Appendix A — Full QUAST `report.txt` for the three Run 2 assemblies

All three: QUAST 5.3.0, run from `/home/abhiraj_arya/research/GNNome`, reference `data/CP017100.fasta`, `MODE: default, threads: 1, min contig length: 500, min alignment length: 65, min alignment IDY: 95.0, ambiguity: one, min local misassembly length: 200, min extensive misassembly length: 1000`. *"All statistics are based on contigs of size >= 500 bp, unless otherwise noted."*

| Metric | Raven + Racon ×2 (`round2`) | GNNome + Racon ×2 (`round2`) | Raven `-p 2` (`raven_p2_assembly`) |
|---|---:|---:|---:|
| QUAST started (UTC) | 2026-11-12 18:20:45 | 2026-11-12 18:20:50 | 2026-11-12 18:20:54 |
| # contigs (>= 0 bp) | 1 | 1 | 1 |
| # contigs (>= 1000 bp) | 1 | 1 | 1 |
| # contigs (>= 5000 bp) | 1 | 1 | 1 |
| # contigs (>= 10000 bp) | 1 | 1 | 1 |
| # contigs (>= 25000 bp) | 1 | 1 | 1 |
| # contigs (>= 50000 bp) | 1 | 1 | 1 |
| Total length (>= 0 bp) | 4578562 | 4578651 | 4578923 |
| Total length (>= 1000 bp) | 4578562 | 4578651 | 4578923 |
| Total length (>= 5000 bp) | 4578562 | 4578651 | 4578923 |
| Total length (>= 10000 bp) | 4578562 | 4578651 | 4578923 |
| Total length (>= 25000 bp) | 4578562 | 4578651 | 4578923 |
| Total length (>= 50000 bp) | 4578562 | 4578651 | 4578923 |
| # contigs | 1 | 1 | 1 |
| Largest contig | 4578562 | 4578651 | 4578923 |
| Total length | 4578562 | 4578651 | 4578923 |
| Reference length | 4583637 | 4583637 | 4583637 |
| GC (%) | 50.77 | 50.77 | 50.77 |
| Reference GC (%) | 50.77 | 50.77 | 50.77 |
| N50 | 4578562 | 4578651 | 4578923 |
| NG50 | 4578562 | 4578651 | 4578923 |
| N90 | 4578562 | 4578651 | 4578923 |
| NG90 | 4578562 | 4578651 | 4578923 |
| auN | 4578562.0 | 4578651.0 | 4578923.0 |
| auNG | 4573492.6 | 4573670.4 | 4574213.8 |
| L50 / LG50 / L90 / LG90 | 1 / 1 / 1 / 1 | 1 / 1 / 1 / 1 | 1 / 1 / 1 / 1 |
| # misassemblies | 0 | 0 | 0 |
| # misassembled contigs | 0 | 0 | 0 |
| Misassembled contigs length | 0 | 0 | 0 |
| # local misassemblies | 0 | 0 | 0 |
| # scaffold gap ext. mis. | 0 | 0 | 0 |
| # scaffold gap loc. mis. | 0 | 0 | 0 |
| # unaligned mis. contigs | 0 | 0 | 0 |
| # unaligned contigs | 0 + 0 part | 0 + 0 part | 0 + 0 part |
| Unaligned length | 0 | 0 | 0 |
| Genome fraction (%) | 99.867 | 99.867 | 99.867 |
| Duplication ratio | 1.000 | 1.000 | 1.000 |
| # N's per 100 kbp | 0.00 | 0.00 | 0.00 |
| # mismatches per 100 kbp | 0.74 | 0.76 | 1.09 |
| # indels per 100 kbp | 28.66 | 29.88 | 35.25 |
| Largest alignment | 4578534 | 4578623 | 4578890 |
| Total aligned length | 4578534 | 4578623 | 4578890 |
| NA50 | 4578534 | 4578623 | 4578890 |
| NGA50 | 4578534 | 4578623 | 4578890 |
| NA90 | 4578534 | 4578623 | 4578890 |
| NGA90 | 4578534 | 4578623 | 4578890 |
| auNA | 4578506.0 | 4578595.0 | 4578857.0 |
| auNGA | 4573436.7 | 4573614.5 | 4574147.9 |
| LA50 / LGA50 / LA90 / LGA90 | 1 / 1 / 1 / 1 | 1 / 1 / 1 / 1 | 1 / 1 / 1 / 1 |
| **From `misassemblies_report.txt`:** # mismatches | 34 | 35 | 50 |
| # indels | 1312 | 1368 | 1614 |
| # indels (<= 5 bp) | 1312 | 1367 | 1614 |
| # indels (> 5 bp) | 0 | 1 | 0 |
| Indels length | 1488 | 1543 | 1827 |

The original files are `R2/07_polish/quast/{raven_racon2,gnnome_racon2,raven_internal_p2}/report.txt` (plus `report.tsv`, `report.html`, `report.pdf`, `icarus.html`, `quast.log`, `contigs_reports/`).
