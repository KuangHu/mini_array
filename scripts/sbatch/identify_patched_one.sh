#!/bin/bash
#SBATCH --job-name=cri_id_patch
#SBATCH --output=/home/kuangh/logs/cri_id_patch_%j.out
#SBATCH --error=/home/kuangh/logs/cri_id_patch_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --qos=standard
#SBATCH --partition=standard
#
# Patched identify runner for the main 100-genome pilot. Writes outputs under
# <P>/raw_patched/identify/<G>/. Does not touch the original raw/identify tree.
# Reuses the existing cctyper outputs (input FASTAs + cctyper version unchanged).

set -euo pipefail
G="${1:?Usage: sbatch identify_patched_one.sh <GENOME_UID>}"
P="${PROJECT_ROOT:-/groups/rubin/projects/kuang/out/CRISPR_array_related}"
T="${SLURM_CPUS_PER_TASK:-4}"

RUN="$P/work_patched/identify/$G"
mkdir -p "$RUN/cwd" "$RUN/json" "$P/raw_patched/identify" "$P/logs"

export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS="$T"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$RUN/cwd"

conda run --no-capture-output -p "$P/envs/identify" \
  python "$P/software/CRISPRidentify.patched/CRISPRidentify.py" \
  --file "$P/input/staged/$G.fna" \
  --result_folder "$P/raw_patched/identify/$G" \
  --json_report "$RUN/json" \
  --fasta_report True \
  --model ALL \
  --cas False \
  --is_element False \
  --strand True \
  --cpu "$T" \
  --parallel True \
  --fast_run False \
  --min_repeats 3 \
  --min_len_rep 21 \
  --max_len_rep 55 \
  --min_len_spacer 18 \
  --max_len_spacer 78 \
  > "$P/logs/identify_patched.$G.log" 2>&1
