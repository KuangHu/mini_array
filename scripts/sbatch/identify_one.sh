#!/bin/bash
#SBATCH --job-name=cri_identify
#SBATCH --output=/home/kuangh/logs/cri_identify_%j.out
#SBATCH --error=/home/kuangh/logs/cri_identify_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --qos=standard
#SBATCH --partition=standard
#
# Run CRISPRidentify on one staged FASTA.
#
# Usage (submitted from a login node):
#   sbatch scripts/sbatch/identify_one.sh <GENOME_UID>
#
# Expects:
#   - project root at /groups/rubin/projects/kuang/out/CRISPR_array_related
#   - staged FASTA at $P/input/staged/<GENOME_UID>.fna
#   - envs/identify built per CRISPRidentify/environment.yml
#
# Parameters mirror step1.md §3 and config/step1.yaml (identify section).

set -euo pipefail

G="${1:?Usage: sbatch identify_one.sh <GENOME_UID>}"
P="${PROJECT_ROOT:-/groups/rubin/projects/kuang/out/CRISPR_array_related}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

RUN="$P/work/identify/$G"
mkdir -p "$RUN/cwd" "$RUN/json" "$P/raw/identify" "$P/logs"

# Keep TF/BLAS single-threaded; parallelism is handled by --cpu.
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$RUN/cwd"

conda run --no-capture-output -p "$P/envs/identify" \
  python "$P/software/CRISPRidentify/CRISPRidentify.py" \
  --file "$P/input/staged/$G.fna" \
  --result_folder "$P/raw/identify/$G" \
  --json_report "$RUN/json" \
  --fasta_report True \
  --model ALL \
  --cas False \
  --is_element False \
  --strand True \
  --cpu "$THREADS" \
  --parallel True \
  --fast_run False \
  --min_repeats 3 \
  --min_len_rep 21 \
  --max_len_rep 55 \
  --min_len_spacer 18 \
  --max_len_spacer 78 \
  > "$P/logs/identify.$G.log" 2>&1
