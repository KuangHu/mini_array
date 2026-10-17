#!/bin/bash
#SBATCH --job-name=release
#SBATCH --output=/home/kuangh/logs/release_%j.out
#SBATCH --error=/home/kuangh/logs/release_%j.err
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --qos=standard
#SBATCH --partition=standard
set -euo pipefail
P=/groups/rubin/projects/kuang/out/CRISPR_array_related
conda run --no-capture-output -p "$P/envs/driver" \
  python /home/kuangh/tools/CRISPR_array_related/scripts/build_release.py \
  > "$P/logs/release.$SLURM_JOB_ID.log" 2>&1
