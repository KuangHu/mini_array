#!/bin/bash
#SBATCH --job-name=coord_fin
#SBATCH --output=/home/kuangh/logs/coord_fin_%j.out
#SBATCH --error=/home/kuangh/logs/coord_fin_%j.err
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --qos=standard
#SBATCH --partition=standard
set -euo pipefail
P=/groups/rubin/projects/kuang/out/CRISPR_array_related
conda run --no-capture-output -p "$P/envs/driver" \
  python /home/kuangh/tools/CRISPR_array_related/scripts/finalize_coordinate_catalog.py \
  > "$P/logs/coord_finalize.$SLURM_JOB_ID.log" 2>&1
