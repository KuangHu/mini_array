#!/bin/bash
#SBATCH --job-name=mini_catalog_v2
#SBATCH --output=/home/kuangh/logs/mini_catalog_v2_%j.out
#SBATCH --error=/home/kuangh/logs/mini_catalog_v2_%j.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --qos=standard
#SBATCH --partition=standard
set -euo pipefail
P="${PROJECT_ROOT:-/groups/rubin/projects/kuang/out/CRISPR_array_related}"
CONFIG="${CONFIG:-/home/kuangh/tools/CRISPR_array_related/config/step2.yaml}"
conda run --no-capture-output -p "$P/envs/driver" \
  python /home/kuangh/tools/CRISPR_array_related/scripts/build_candidate_catalog.py \
  --config "$CONFIG"
