#!/bin/bash
#SBATCH --job-name=cctyper
#SBATCH --output=/home/kuangh/logs/cctyper_%j.out
#SBATCH --error=/home/kuangh/logs/cctyper_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --qos=standard
#SBATCH --partition=standard
#
# Run CRISPRCasTyper on one staged FASTA.
#
# Usage (submitted from a login node):
#   sbatch scripts/sbatch/cctyper_one.sh <GENOME_UID> <PRODIGAL_MODE>
#     PRODIGAL_MODE ∈ {single, meta}
#       - single for isolates
#       - meta   for MAG / phage / plasmid / unknown
#
# Parameters mirror step1.md §8 and config/step1.yaml (cas_annotation section).

set -euo pipefail

G="${1:?Usage: sbatch cctyper_one.sh <GENOME_UID> <single|meta>}"
MODE="${2:?Specify prodigal mode: single or meta}"
case "$MODE" in
  single|meta) ;;
  *) echo "PRODIGAL_MODE must be single or meta, got $MODE" >&2; exit 2 ;;
esac
P="${PROJECT_ROOT:-/groups/rubin/projects/kuang/out/CRISPR_array_related}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

mkdir -p "$P/raw/cctyper" "$P/logs"

# cctyper fails if the output dir already exists; make sure we clean a half run.
OUT="$P/raw/cctyper/$G"
if [ -d "$OUT" ]; then
  echo "existing $OUT — cctyper will refuse to overwrite; remove or move before resubmitting" >&2
  exit 3
fi

conda run --no-capture-output -p "$P/envs/driver" \
  cctyper "$P/input/staged/$G.fna" "$OUT" \
  --threads "$THREADS" \
  --prodigal "$MODE" \
  --ccd 10000 \
  --skip_blast \
  --no_plot \
  > "$P/logs/cctyper.$G.log" 2>&1
