#!/bin/bash
# One-genome smoke submission. Submits identify + cctyper for the first manifest row.
# Run from a Rubin login node.
set -euo pipefail
cd "$(dirname "$0")/../.."
IDJOB=$(sbatch --parsable scripts/sbatch/identify_one.sh "GCA_001803565_1")
echo "identify jobid=$IDJOB"
CCJOB=$(sbatch --parsable scripts/sbatch/cctyper_one.sh "GCA_001803565_1" "meta")
echo "cctyper  jobid=$CCJOB"
echo "Watch with: squeue -j $IDJOB,$CCJOB"
echo "Expected outputs:"
echo "  /groups/rubin/projects/kuang/out/CRISPR_array_related/raw/identify/GCA_001803565_1/"
echo "  /groups/rubin/projects/kuang/out/CRISPR_array_related/raw/cctyper/GCA_001803565_1/"
echo "  /groups/rubin/projects/kuang/out/CRISPR_array_related/logs/identify.GCA_001803565_1.log"
echo "  /groups/rubin/projects/kuang/out/CRISPR_array_related/logs/cctyper.GCA_001803565_1.log"
