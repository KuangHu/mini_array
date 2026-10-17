#!/bin/bash
#SBATCH --job-name=step2c_rp
#SBATCH --output=/home/kuangh/logs/step2c_rp_%j.out
#SBATCH --error=/home/kuangh/logs/step2c_rp_%j.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --qos=standard
#SBATCH --partition=standard
set -euo pipefail
P=/groups/rubin/projects/kuang/out/CRISPR_array_related
conda run --no-capture-output -p "$P/envs/driver" \
  python /home/kuangh/tools/CRISPR_array_related/scripts/step2c_cards_and_report.py \
  > "$P/logs/step2c_report.$SLURM_JOB_ID.log" 2>&1
