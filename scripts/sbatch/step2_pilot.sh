#!/bin/bash
#SBATCH --job-name=step2_pilot
#SBATCH --output=/home/kuangh/logs/step2_pilot_%j.out
#SBATCH --error=/home/kuangh/logs/step2_pilot_%j.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --qos=standard
#SBATCH --partition=standard
#
# Step 2 pilot — runs the complete search + pairing + dedup + context on the
# frozen 1000-genome batch and the 213 primary HQ queries. Executes under
# envs/identify (which has blastn + makeblastdb) but the Python glue uses
# envs/driver imports — we call step2.py through driver and point BLAST calls
# at identify's binaries via config.
#
# Explicit caps per user instruction:
#   - cpus_per_job=8
#   - mem_per_job_gb=32
#   - max_concurrent_jobs=4 (this script is a single job; the pilot driver is
#     single-process; the "4" cap applies if we fan-out subsidiary jobs later)
#   - wall_time_h=6

set -euo pipefail
P="${PROJECT_ROOT:-/groups/rubin/projects/kuang/out/CRISPR_array_related}"
CONFIG="${CONFIG:-/home/kuangh/tools/CRISPR_array_related/config/step2.yaml}"

mkdir -p "$P/results/step2_pilot" "$P/logs"

conda run --no-capture-output -p "$P/envs/driver" \
  python /home/kuangh/tools/CRISPR_array_related/scripts/step2.py \
  run --config "$CONFIG" \
  > "$P/logs/step2_pilot.$SLURM_JOB_ID.log" 2>&1
