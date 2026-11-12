#!/usr/bin/env bash
set -euo pipefail
in=$1; out=$2; reads=$3; P=$HOME/miniconda3/envs/polish/bin
mkdir -p "$out"; cp "$in" "$out/round0.fasta"
for i in 1 2; do
  echo "round $i: minimap2 -t 4 -x map-pb round$((i-1)).fasta reads | racon -t 4 -m 3 -x -5 -g -4 -w 500" | tee -a "$out/commands.txt"
  $P/minimap2 -t 4 -x map-pb "$out/round$((i-1)).fasta" "$reads" > "$out/round$i.paf" 2>> "$out/polish.log"
  $P/racon -t 4 -m 3 -x -5 -g -4 -w 500 "$reads" "$out/round$i.paf" "$out/round$((i-1)).fasta" > "$out/round$i.fasta" 2>> "$out/polish.log"
done
echo "done: $out" | tee -a "$out/commands.txt"
