#!/bin/bash
#SBATCH --job-name=step2c_an
#SBATCH --output=/home/kuangh/logs/step2c_an_%j.out
#SBATCH --error=/home/kuangh/logs/step2c_an_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --qos=standard
#SBATCH --partition=standard
set -euo pipefail
P=/groups/rubin/projects/kuang/out/CRISPR_array_related
conda run --no-capture-output -p "$P/envs/driver" \
  python /home/kuangh/tools/CRISPR_array_related/scripts/step2c_analyze.py \
  > "$P/logs/step2c_analyze.$SLURM_JOB_ID.log" 2>&1
