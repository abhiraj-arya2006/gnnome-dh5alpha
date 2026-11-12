#!/usr/bin/env bash
# Same protocol as Run 2 (run2_corrected/scripts/polish2.sh): minimap2 -x map-pb + racon -m 3 -x -5 -g -4 -w 500, 2 rounds.
set -euo pipefail
G=/home/abhiraj_arya/research/GNNome; O=/home/abhiraj_arya/research/analysis/10_polishing/run1; PB=/home/abhiraj_arya/miniconda3/envs/polish/bin
READS=$G/data/CP017100_sim.fastq; QUAST=/home/abhiraj_arya/miniconda3/envs/gnnome/bin/quast.py
for x in raven:$G/data/CP017100_assembly.fasta gnnome:$G/data/raven/assembly/0_assembly.fasta; do
  o=$O/${x%%:*}; mkdir -p $o; cp ${x#*:} $o/round0.fasta
  for i in 1 2; do
    echo "round $i: minimap2 -t 2 -x map-pb | racon -t 2 -m 3 -x -5 -g -4 -w 500" >> $o/commands.txt
    $PB/minimap2 -t 2 -x map-pb $o/round$((i-1)).fasta $READS > $o/round$i.paf 2>> $o/polish.log
    $PB/racon -t 2 -m 3 -x -5 -g -4 -w 500 $READS $o/round$i.paf $o/round$((i-1)).fasta > $o/round$i.fasta 2>> $o/polish.log
    rm -f $o/round$i.paf
    (cd $G && $QUAST -t 2 $o/round$i.fasta -r data/CP017100.fasta -o $o/quast_round$i > /dev/null 2>&1)
  done
done
echo DONE > $O/DONE
