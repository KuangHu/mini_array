#!/bin/bash
#SBATCH --job-name=smk_id_patch
#SBATCH --output=/home/kuangh/logs/smk_id_patch_%j.out
#SBATCH --error=/home/kuangh/logs/smk_id_patch_%j.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --qos=standard
#SBATCH --partition=standard
#
# Same as smoke_identify_one.sh but points at software/CRISPRidentify.patched/
# and writes outputs under smoke_test/raw_patched/ (does not touch the
# previously saved unpatched outputs under smoke_test/raw/).
set -euo pipefail
G="${1:?Usage: sbatch smoke_identify_patched_one.sh <GENOME_UID>}"
P="${PROJECT_ROOT:-/groups/rubin/projects/kuang/out/CRISPR_array_related}"
SMK="$P/smoke_test"
T="${SLURM_CPUS_PER_TASK:-4}"
RUN="$SMK/work_patched/identify/$G"
mkdir -p "$RUN/cwd" "$RUN/json" "$SMK/raw_patched/identify" "$SMK/logs"
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS="$T"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$RUN/cwd"
conda run --no-capture-output -p "$P/envs/identify" \
  python "$P/software/CRISPRidentify.patched/CRISPRidentify.py" \
  --file "$SMK/input/staged/$G.fna" \
  --result_folder "$SMK/raw_patched/identify/$G" \
  --json_report "$RUN/json" \
  --fasta_report True --model ALL --cas False --is_element False \
  --strand True --cpu "$T" --parallel True --fast_run False \
  --min_repeats 3 --min_len_rep 21 --max_len_rep 55 \
  --min_len_spacer 18 --max_len_spacer 78 \
  > "$SMK/logs/identify_patched.$G.log" 2>&1
