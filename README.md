# GNNome on *E. coli* DH5α — verified evaluation and mechanism analysis

Evaluation of the **pretrained GNNome** graph-neural-network assembler (SymGatedGCN, no retraining) against **Raven** on *E. coli* DH5α (CP017100.1), using simulated long reads, with ground-truth tracing of every graph edge and every assembly error.

**Start here:** [`docs/GNNome_DH5alpha_results_update.md`](docs/GNNome_DH5alpha_results_update.md), the full results with all figures.
Provenance and integrity evidence: [`docs/run1_run2_evidence_report.md`](docs/run1_run2_evidence_report.md).

## Key results

| Finding | Evidence |
|---|---|
| The earlier DH5α GNNome result (manuscript v28) was built from a stale **K-12** graph; its "repeat-linked failure" is the IS3-mediated ~210 kb inversion separating K-12 from DH5α | Edge-by-edge graph match; identical breakpoint set in the reference-vs-reference alignment |
| Unpolished, **GNNome is more accurate than Raven**: Run 1 has 44 % fewer mismatches and 41 % fewer indels; better in 21/23 replicates (p = 1.3×10⁻⁵) | `results/12_statistics/` |
| The cause is **read selection**: GNNome copies more accurate reads (23/23 replicates), which predicts the error counts to within 3–5 % | `results/09_repeat_errors/`, Fig. 5, Fig. 8 |
| GNNome's large indels come from a **parser bug**: it reuses Raven's one-directional overlap length for the reverse-complement edge (r = 0.998; 35/35 indels; 15/15 junction misassemblies) | Fig. 4 |
| GNNome **never used a repeat-induced false overlap** (2 runs + 23 replicates) | Fig. 3, Fig. 6 |
| **Polishing equalizes** the assemblers (Run 1: sequence-identical, error-free) | Fig. 7 |

## Repository layout

```
docs/                  results write-up and provenance/evidence report
figures/               F1–F9 (PNG 300 dpi + PDF)
scripts/analysis/      analysis pipeline (tasks 1–12), figures, statistics, verify_all.py
scripts/run2_pipeline/ helper scripts for the corrected Run 2 (seed search, graph proof, polishing)
results/               small result tables / JSON per task (01_inventory … 12_statistics)
evidence/              hash manifests, command logs and proofs for Run 1 and Run 2
DATA_MANIFEST.tsv      SHA256 of every large input/output kept out of git (reads, graphs, assemblies)
environment.md         software versions
```

## What is *not* in the repo
Reads, graphs, assemblies and the 24 sensitivity-replicate run folders (~8 GB) are excluded. Each can be identified by its SHA256 in `DATA_MANIFEST.tsv` and regenerated:
- Run 2 reads: PBSIM3 seed `1783636376` (command in `evidence/run2/pbsim3_seed_and_command.txt`; regeneration is byte-identical).
- Sensitivity replicates: `scripts/analysis/run_replicate.sh` with `results/06_sensitivity/jobs.txt`.

## Reproducing
Scripts use absolute paths of the original workstation (`/home/abhiraj_arya/research/...`: `GNNome/`, `run2_corrected/`, `analysis/`). Adjust these paths, then run, for example:
```bash
python scripts/analysis/verify_all.py      # re-checks every claim (35/35 PASS on the original machine)
python scripts/analysis/figures.py         # regenerates F1–F7
python scripts/analysis/stats_sensitivity.py results/06_sensitivity/collected.tsv   # F8, F9 + statistics
```

## Limitations
Simulated reads only; one genome; the available PBSIM3 build caps simulated accuracy at ~96 % (the HiFi regime is covered by the Run 1 reads); the parser-fix benefit is predicted, not yet run.
