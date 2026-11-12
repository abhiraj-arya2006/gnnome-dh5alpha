#!/usr/bin/env bash
# PBSIM3 uses seed = time(NULL) when --seed is not given. Try every Unix second in a window
# around the Run 2 simulation and keep the seed whose first read equals the original first read.
set -euo pipefail
G=$HOME/research/GNNome; TMP=$(mktemp -d)
T=$(zcat $G/data/dh5alpha/simulated_0001.fq.gz | sed -n 2p)
for s in $(seq "$1" "$2"); do
  $G/vendor/pbsim3/src/pbsim --strategy wgs --method qshmm --qshmm $G/vendor/pbsim3/data/QSHMM-RSII.model \
    --depth 0.01 --genome $G/data/dh5alpha/CP017100.fasta --prefix $TMP/s --seed $s 2>/dev/null
  [ "$(zcat $TMP/s_0001.fq.gz | sed -n 2p)" = "$T" ] && { echo "MATCH seed=$s ($(date -u -d @$s))"; rm -rf $TMP; exit 0; }
done
echo "no match"; rm -rf $TMP; exit 1
