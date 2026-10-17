#!/bin/bash
#SBATCH --job-name=step2b
#SBATCH --output=/home/kuangh/logs/step2b_%j.out
#SBATCH --error=/home/kuangh/logs/step2b_%j.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --qos=standard
#SBATCH --partition=standard
set -euo pipefail
P="${PROJECT_ROOT:-/groups/rubin/projects/kuang/out/CRISPR_array_related}"
CONFIG="${CONFIG:-/home/kuangh/tools/CRISPR_array_related/config/step2.yaml}"
mkdir -p "$P/results/step2b_pilot" "$P/logs"
conda run --no-capture-output -p "$P/envs/driver" \
  python /home/kuangh/tools/CRISPR_array_related/scripts/step2b.py \
  --config "$CONFIG" \
  > "$P/logs/step2b_explainer.$SLURM_JOB_ID.log" 2>&1
