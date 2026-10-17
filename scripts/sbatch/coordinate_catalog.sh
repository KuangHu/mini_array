#!/bin/bash
#SBATCH --job-name=coord_cat
#SBATCH --output=/home/kuangh/logs/coord_cat_%j.out
#SBATCH --error=/home/kuangh/logs/coord_cat_%j.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --qos=standard
#SBATCH --partition=standard
set -euo pipefail
P=/groups/rubin/projects/kuang/out/CRISPR_array_related
conda run --no-capture-output -p "$P/envs/driver" \
  python /home/kuangh/tools/CRISPR_array_related/scripts/build_coordinate_catalog.py \
  > "$P/logs/coord_catalog.$SLURM_JOB_ID.log" 2>&1
