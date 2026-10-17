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
P=/groups/rubin/projects/kuang/out/CRISPR_array_related

# Step 2C phases: dedup + review + proteomes (prodigal × 1000 genomes)
# downstream phases (families/homologs/associations/report) in step2c_analyze.sh
conda run --no-capture-output -p "$P/envs/driver" \
  python /home/kuangh/tools/CRISPR_array_related/scripts/step2c.py \
  --phases dedup review proteomes \
  > "$P/logs/step2c_prepare.$SLURM_JOB_ID.log" 2>&1
