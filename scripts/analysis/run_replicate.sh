#!/usr/bin/env bash
# One sensitivity replicate: PBSIM3 -> Raven fork (-p0, graph) -> GNNome (pretrained, no training)
# -> truth/edge analysis -> identity + QUAST (unpolished: --min-identity 80) -> Racon x2 both -> QUAST default.
# usage: run_replicate.sh <cond> <accuracy> <depth> <len_mean> <len_sd> <seed> <threads>
set -euo pipefail
cond=$1 acc=$2 depth=$3 lm=$4 lsd=$5 seed=$6 T=$7
G=$HOME/research/GNNome; A=$HOME/research/analysis; PY=$HOME/miniconda3/envs/gnnome/bin/python
QUAST=$HOME/miniconda3/envs/gnnome/bin/quast.py; MM=$HOME/miniconda3/envs/gnnome/bin/minimap2; PB=$HOME/miniconda3/envs/polish/bin
REF=$G/data/CP017100.fasta; RAVEN=$G/vendor/raven-1.8.1/build/bin/raven
D=$A/06_sensitivity/runs/${cond}_s${seed}; rm -rf "$D"; mkdir -p "$D"/{sim,raven,gnnome,polish}; cd "$D"
export OMP_NUM_THREADS=$T MKL_NUM_THREADS=$T
log(){ echo "[$(date -u +%H:%M:%S)] $*" >> "$D/pipeline.log"; }
echo "cond=$cond acc=$acc depth=$depth len_mean=$lm len_sd=$lsd seed=$seed threads=$T" > params.txt
log "pbsim"; "$G/vendor/pbsim3/src/pbsim" --strategy wgs --method qshmm --qshmm "$G/vendor/pbsim3/data/QSHMM-RSII.model" \
  --depth "$depth" --accuracy-mean "$acc" --length-mean "$lm" --length-sd "$lsd" --genome "$REF" --prefix sim/sim --seed "$seed" 2> sim/pbsim.log
READS=$D/sim/sim_0001.fq.gz
log "raven"; (cd raven && "$RAVEN" -t "$T" -p0 --disable-checkpoints "$READS" > raven_p0.fasta 2> raven.log)
test -s raven/graph_1.gfa && test -s raven/raven_p0.fasta
log "gnnome"; (cd "$G" && "$PY" create_inference_graphs.py --reads "$READS" --gfa "$D/raven/graph_1.gfa" --asm raven --out "$D/gnnome" > "$D/gnnome/preprocess.log" 2>&1 \
  && "$PY" inference.py --data "$D/gnnome" --asm raven --out "$D/gnnome/raven" --model weights/weights.pt > "$D/gnnome/inference.log" 2>&1)
test -s gnnome/raven/assembly/0_assembly.fasta
log "truth"; "$PY" -I "$A/scripts/read_truth.py" sim/sim_0001.maf.gz "$READS" 4583637 sim/truth > /dev/null
"$PY" -I "$A/scripts/graph_truth.py" raven/graph_1.gfa gnnome/raven sim/truth.reads.tsv gnnome/edges > /dev/null 2>&1
log "identity"; for x in raven:raven/raven_p0.fasta gnnome:gnnome/raven/assembly/0_assembly.fasta; do
  "$MM" -t "$T" -cx asm20 --secondary=no "$REF" "${x#*:}" 2>/dev/null | cut -f1-12 | awk -v n=${x%%:*} '{al+=$11; m+=$10; tc+=$9-$8; q=$2} END{printf "%s\t%d\t%d\t%.6f\n", n, q, tc, (al?m/al:0)}'; done > identity_unpolished.tsv
log "quast unpolished"; for x in raven:raven/raven_p0.fasta gnnome:gnnome/raven/assembly/0_assembly.fasta; do
  "$QUAST" -t "$T" --min-identity 80 "${x#*:}" -r "$REF" -o "quast_p0_minid80/${x%%:*}" > /dev/null 2>&1; done
log "polish"; for x in raven:raven/raven_p0.fasta gnnome:gnnome/raven/assembly/0_assembly.fasta; do
  o=polish/${x%%:*}; mkdir -p $o; cp "${x#*:}" $o/round0.fasta
  for i in 1 2; do "$PB/minimap2" -t "$T" -x map-pb $o/round$((i-1)).fasta "$READS" > $o/round$i.paf 2>> $o/polish.log
    "$PB/racon" -t "$T" -m 3 -x -5 -g -4 -w 500 "$READS" $o/round$i.paf $o/round$((i-1)).fasta > $o/round$i.fasta 2>> $o/polish.log; rm -f $o/round$i.paf; done
  test -s $o/round2.fasta; done
log "quast polished"; for x in raven gnnome; do "$QUAST" -t "$T" polish/$x/round2.fasta -r "$REF" -o quast_racon2/$x > /dev/null 2>&1; done
rm -f gnnome/raven/info/0_reads.pkl   # large, regenerable
log "DONE"; touch DONE
