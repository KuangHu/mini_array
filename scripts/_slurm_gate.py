"""Shared SLURM gate — import `require_slurm()` at the top of any analysis
entry point. Running outside sbatch should hard-fail immediately."""
from __future__ import annotations

import os
import sys


def require_slurm() -> str:
    jid = os.environ.get("SLURM_JOB_ID")
    if not jid:
        sys.stderr.write(
            "ERROR: SLURM_JOB_ID is not set.\n"
            "This script is a workload and MUST run inside an sbatch job.\n"
            "Submit via `sbatch scripts/sbatch/<wrapper>.sh` instead.\n"
        )
        sys.exit(78)   # EX_CONFIG
    return jid
