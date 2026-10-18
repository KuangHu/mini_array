#!/bin/bash
#SBATCH --job-name=step2c_prep
#SBATCH --output=/home/kuangh/logs/step2c_prep_%j.out
#SBATCH --error=/home/kuangh/logs/step2c_prep_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --qos=standard
#SBATCH --partition=standard

set -euo pipefail
P="${PROJECT_ROOT:-/groups/rubin/projects/kuang/out/CRISPR_array_related}"
CONFIG="${CONFIG:-/home/kuangh/tools/CRISPR_array_related/config/step2.yaml}"
mkdir -p "$P/logs"

# Catalog preparation needs only dedup + structural review.
# Optional protein analyses require a separate explicit proteomes run.
conda run --no-capture-output -p "$P/envs/driver" \
  python /home/kuangh/tools/CRISPR_array_related/scripts/step2c.py \
  --config "$CONFIG" --phases dedup review \
  > "$P/logs/step2c_prepare.$SLURM_JOB_ID.log" 2>&1
