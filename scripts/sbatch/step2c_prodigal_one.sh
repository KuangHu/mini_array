#!/bin/bash
#SBATCH --job-name=step2c_pr
#SBATCH --output=/home/kuangh/logs/step2c_pr_%A_%a.out
#SBATCH --error=/home/kuangh/logs/step2c_pr_%A_%a.err
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --qos=standard
#SBATCH --partition=standard
#
# Array task: runs prodigal on one staged FASTA.
# Expects env vars GENOMES_FILE (newline-separated "guid<TAB>mode") and
# OUT_DIR (destination for .faa / .gff).
set -euo pipefail
GENOMES_FILE="${GENOMES_FILE:?GENOMES_FILE unset}"
OUT_DIR="${OUT_DIR:?OUT_DIR unset}"
P=/groups/rubin/projects/kuang/out/CRISPR_array_related

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$GENOMES_FILE")
GUID=$(echo "$LINE" | cut -f1)
MODE=$(echo "$LINE" | cut -f2)
if [ -z "$GUID" ] || [ -z "$MODE" ]; then
  echo "no work (empty line)"; exit 0
fi

FAA="$OUT_DIR/$GUID.faa"
GFF="$OUT_DIR/$GUID.gff"
if [ -s "$FAA" ]; then
  echo "already done: $GUID"
  exit 0
fi

mkdir -p "$OUT_DIR"
"$P/envs/identify/bin/prodigal" \
  -i "$P/input/staged/$GUID.fna" \
  -p "$MODE" \
  -a "$FAA" -f gff -o "$GFF" -q
echo "ok $GUID"
