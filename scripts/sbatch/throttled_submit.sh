#!/bin/bash
# Throttled submitter: reads a file of sbatch commands and dispatches them
# while the user's queue depth stays below MAX_QUEUED. Prints each submitted
# jobid. Written for the batch-1000 scale-up where the SLURM QOS caps user
# submissions at 200.
#
# Usage:
#   bash throttled_submit.sh <commands_file> [MAX_QUEUED] [SLEEP_SEC]

set -eo pipefail   # (no -u) commands may reference $IDENTIFY / $CCTYPER set below
CMDS_FILE="${1:?Usage: throttled_submit.sh commands_file [MAX_QUEUED] [SLEEP_SEC]}"
MAX_QUEUED="${2:-180}"      # keep headroom under the 200 QOS limit
SLEEP_SEC="${3:-20}"
TOTAL=$(grep -c '^sbatch ' "$CMDS_FILE" || echo 0)
export IDENTIFY=/home/kuangh/tools/CRISPR_array_related/scripts/sbatch/identify_patched_one.sh
export CCTYPER=/home/kuangh/tools/CRISPR_array_related/scripts/sbatch/cctyper_one.sh

echo "throttled_submit: dispatching $TOTAL sbatch commands, keeping queue <= $MAX_QUEUED"
submitted=0
grep '^sbatch ' "$CMDS_FILE" | while read -r line; do
  while :; do
    q=$(squeue -u "$USER" -h 2>/dev/null | wc -l)
    if [ "$q" -lt "$MAX_QUEUED" ]; then
      break
    fi
    sleep "$SLEEP_SEC"
  done
  # Resolve shell variables ($IDENTIFY, $CCTYPER) that the generator left in-place
  eval "jid=\$($line)"
  submitted=$((submitted + 1))
  echo "$jid"
done
echo "done: $submitted/$TOTAL submitted"
