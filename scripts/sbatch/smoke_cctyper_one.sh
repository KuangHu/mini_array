#!/bin/bash
#SBATCH --job-name=smk_cctyper
#SBATCH --output=/home/kuangh/logs/smk_cctyper_%j.out
#SBATCH --error=/home/kuangh/logs/smk_cctyper_%j.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --qos=standard
#SBATCH --partition=standard
set -euo pipefail
G="${1:?Usage: sbatch smoke_cctyper_one.sh <GENOME_UID> <single|meta>}"
MODE="${2:?Specify prodigal mode: single or meta}"
case "$MODE" in single|meta) ;; *) echo "bad MODE: $MODE" >&2; exit 2;; esac
P="${PROJECT_ROOT:-/groups/rubin/projects/kuang/out/CRISPR_array_related}"
SMK="$P/smoke_test"
T="${SLURM_CPUS_PER_TASK:-4}"
mkdir -p "$SMK/raw/cctyper" "$SMK/logs"
OUT="$SMK/raw/cctyper/$G"
if [ -d "$OUT" ]; then
  echo "existing $OUT — remove before re-running" >&2; exit 3
fi
conda run --no-capture-output -p "$P/envs/driver" \
  cctyper "$SMK/input/staged/$G.fna" "$OUT" \
  --threads "$T" --prodigal "$MODE" --ccd 10000 --skip_blast --no_plot \
  > "$SMK/logs/cctyper.$G.log" 2>&1
