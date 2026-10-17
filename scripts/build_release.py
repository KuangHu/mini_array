#!/usr/bin/env python3
"""Build the release package under /home/kuangh/tools/CRISPR_array_related/release/.

Does NOT push to GitHub. Produces a ready-to-push git tree. The user runs
`git remote add origin <url> && git push -u origin main` after review.
sbatch gate enforced.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _slurm_gate import require_slurm  # noqa: E402

CODE_ROOT = Path("/home/kuangh/tools/CRISPR_array_related")
DATA_ROOT = Path("/groups/rubin/projects/kuang/out/CRISPR_array_related")
RELEASE = CODE_ROOT / "release"

# Items from the code tree we WANT in the release (relative to CODE_ROOT)
CODE_WANTED = [
    "scripts",
    "config",
    "tests",
    "input_template",
    "step1.md",
    ".gitignore",
]

# Items from the data tree we WANT in the release (relative to DATA_ROOT)
DATA_WANTED = [
    ("results/coordinate_catalog", "results/coordinate_catalog"),
    ("results/sequence_catalog", "results/sequence_catalog"),
    ("results/batch1000_delivery_metrics", "results/batch1000_delivery_metrics"),
    ("results/REPORT_PATCHED.md", "results/REPORT_PATCHED.md"),
    ("results/REPORT_BATCH1000.md", "results/REPORT_BATCH1000.md"),
    ("results/REPORT_STEP2_PILOT.md", "results/REPORT_STEP2_PILOT.md"),
    ("results/REPORT_STEP2B.md", "results/REPORT_STEP2B.md"),
    ("provenance", "provenance"),
]

# Explicit excludes inside code dirs
CODE_EXCLUDES = (
    "__pycache__",
    ".pytest_cache",
    "data_root",            # symlink into /groups — don't follow
    "tmpsq31pz4b",          # one-off scratch dir if present
    "release",              # don't recurse into our own output
    ".git",                 # old .git if any
)


def _log(msg: str) -> None:
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def _copytree_filtered(src: Path, dst: Path) -> None:
    """Copy a directory, dropping cache / symlink / scratch artifacts."""
    def _ignore(_dir, names):
        return [n for n in names if n in CODE_EXCLUDES
                 or n.endswith(".pyc") or n.endswith(".pyo")]
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst, symlinks=False, ignore=_ignore)


def _run(cmd, cwd=None, check=True):
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if check and r.returncode != 0:
        raise RuntimeError(f"{' '.join(map(str, cmd))} failed:\n{r.stderr}")
    return r


def write_readme() -> None:
    txt = f"""# CRISPR array / mini-array coordinate + sequence catalog

Reproducible catalog of CRISPR arrays and mini-array candidates recovered
from the first 1000 bacterial representative genomes of GTDB r220, using
a patched CRISPRidentify fork and CRISPRCasTyper. This repository holds
code, reports, and the frozen coordinate + sequence deliverables.

## Scope

- **Source genomes**: 1000 of the 107,235 bac120 representative genomes
  in GTDB release r220 (public data).
- **Tools**: CRISPRidentify (patched; see `provenance/LOCKED.md`) and
  CRISPRCasTyper 1.9.0.
- **Primary standard**: Step-1 `seed_qc.min_score=0.75`; strict comparator
  at 0.90.
- **Not included here**: no full-database runs, no CDS/protein analysis,
  no RNA-protein association, no biological inference beyond coordinate
  recovery.
- **`allow_full_run: false`** — this is a pilot catalog, not a full-DB
  deliverable.

## Repository layout

```
scripts/                 Pipeline code (step1.py, step2.py, step2b*, step2c*,
                         build_coordinate_catalog.py, build_sequence_catalog.py,
                         export_pairs.py, …) + sbatch wrappers
config/step1.yaml        Frozen Step 1 configuration
config/step2.yaml        Frozen Step 2 configuration
tests/                   Unit + regression tests (46 passing)
step1.md                 Original Step 1 design spec (Chinese)

provenance/              Pinned software + resource fingerprints
  LOCKED.md              Frozen implementation spec (upstream commit, patch
                         hashes, resource SHA256s, env lockfiles)
  CRISPRidentify.commit  Upstream CRISPRidentify git commit used
  CRISPRCasTyper.commit  Upstream CRISPRCasTyper git commit used
  CRISPRloci.commit      CRISPRloci commit that supplied models_tandem.hmm
  patches/               Unified diffs applied to CRISPRidentify fork
  identify.conda.explicit.txt   envs/identify lockfile
  driver.conda.explicit.txt     envs/driver lockfile
  models_tandem.hmm.*    HMM resource provenance
  fingerprint.json       Cache key over patches + resources + envs
  upstream_issue_draft.md Draft issue for CRISPRidentify (not submitted)
  interface_notes.md     Upstream tool quirks observed
  pilot_manifest.build.json Deterministic 1000-genome manifest build record

results/
  coordinate_catalog/            Phase L1+L2 — array + mini-array coordinates
  sequence_catalog/              Phase L3+L4 — genomic-forward DNA for every
                                 call, mini-array, pair unit
  batch1000_delivery_metrics/    Phase H6 — per-run resource accounting (v1-v3)
  REPORT_PATCHED.md              100-genome pilot report after upstream patch
  REPORT_BATCH1000.md            1000-genome scale-up report
  REPORT_STEP2_PILOT.md          Step 2 anchor-pair pilot
  REPORT_STEP2B.md               Step 2B structural review
```

## Reproducing

See `REPRODUCING.md`. All compute runs under SLURM via `scripts/sbatch/`
wrappers; a `SLURM_JOB_ID` gate in every analysis entry point refuses
login-node execution.

## Non-negotiables honoured

- Only source of sequence: the staged FASTAs derived deterministically
  from GTDB r220 reps. No consensus / query substitution.
- All sequences in the catalog are genomic-forward DNA; no automatic
  reverse-complement.
- Ambiguous bases preserved as-is; `has_ambiguous_bases` columns flag
  them without rewriting.
- Mini-array strand is `.` on BED — repeat alignment direction is a
  BLAST artefact, not a biological orientation.
- Pair-level correspondence preserved in `mini_array_units.*`; identical
  component sequences are NOT deduplicated.
- Thresholds (identify `min_score`, Step 2 anchor/secondary cutoffs,
  spacer window, flank requirement) are the frozen L2-finalised values.
- No new biology, no CDS/protein/HMM, no RNA-protein association.
- `allow_full_run: false`.

## License

MIT — see `LICENSE`.

## Caveats

- Query library is small (213 NR100 HQ repeats) — repeat families not
  represented in these 213 queries are invisible to Step 2.
- `repeat_2_completeness='degenerate' = 0` is a consequence of the
  frozen Step 2 secondary thresholds, not a biological claim.
- Mini-array count of 40 reflects the post-filter cut (no identify
  repeat units in ±2 kb window); the raw pair catalogue is 15,472 pairs
  under `results/sequence_catalog/pairs/`.
- CRISPRCasTyper does not publish per-repeat genomic coordinates;
  `array_units_cctyper_unavailable.tsv` labels this explicitly.

## Related

- CRISPRidentify: https://github.com/BackofenLab/CRISPRidentify
- CRISPRCasTyper: https://github.com/Russel88/CRISPRCasTyper
- CRISPRloci (source of models_tandem.hmm):
  https://github.com/BackofenLab/CRISPRloci
- GTDB r220: https://gtdb.ecogenomic.org/

"""
    (RELEASE / "README.md").write_text(txt)


def write_license() -> None:
    year = "2026"
    txt = f"""MIT License

Copyright (c) {year} CRISPR array pilot contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
    (RELEASE / "LICENSE").write_text(txt)


def write_gitignore() -> None:
    (RELEASE / ".gitignore").write_text("""__pycache__/
*.pyc
*.pyo
.pytest_cache/
.mypy_cache/
.ipynb_checkpoints/
*.log
# Big data not shipped in the repo:
input/staged/
raw/
raw_patched/
work/
work_patched/
envs/
software/
logs/
.conda_cache/
# Scratch
*.swp
.DS_Store
""")


def write_reproducing() -> None:
    txt = """# Reproducing the catalog

This repository ships the frozen coordinate + sequence catalog and the
code that produced it. The underlying raw outputs (staged FASTAs, conda
envs, software clones, cctyper per-genome output, patched identify
per-genome output) are NOT shipped — they are large and derivable.

## Prerequisites

- SLURM cluster with `sbatch` reachable from the submit node.
- Conda / mamba.
- GTDB r220 bac120 representative genomes available locally.

## Steps

1. Install envs from the locked explicit files in `provenance/`:
   ```
   mamba create -p envs/identify --file provenance/identify.conda.explicit.txt
   mamba create -p envs/driver   --file provenance/driver.conda.explicit.txt
   ```

2. Pin the three upstream tools to the recorded commits:
   ```
   git clone --recurse-submodules https://github.com/BackofenLab/CRISPRidentify.git software/CRISPRidentify
   git -C software/CRISPRidentify checkout $(cat provenance/CRISPRidentify.commit)

   git clone https://github.com/Russel88/CRISPRCasTyper.git software/CRISPRCasTyper
   git -C software/CRISPRCasTyper checkout $(cat provenance/CRISPRCasTyper.commit)

   git clone --recurse-submodules https://github.com/BackofenLab/CRISPRloci.git software/CRISPRloci
   git -C software/CRISPRloci checkout $(cat provenance/CRISPRloci.commit)
   ```

3. Create the patched CRISPRidentify fork:
   ```
   rsync -a --exclude='.git' software/CRISPRidentify/ software/CRISPRidentify.patched/
   mkdir -p software/CRISPRidentify.patched/tools/hmm_search
   cp software/CRISPRloci/CRISPRidentify/tools/hmm_search/models_tandem.hmm \\
      software/CRISPRidentify.patched/tools/hmm_search/
   patch -p0 -d software/CRISPRidentify.patched/ < provenance/patches/CRISPRidentify.py.patch
   patch -p0 -d software/CRISPRidentify.patched/components/ < provenance/patches/components_evaluation.py.patch
   cp provenance/patches/_resources.py.new software/CRISPRidentify.patched/components/_resources.py
   ```

4. Build the pilot manifest + stage FASTAs:
   ```
   python scripts/build_pilot_manifest.py --config config/step1.yaml --n 1000
   python scripts/step1.py prepare --config config/step1.yaml
   ```

5. Submit Step 1 + Step 2 via `scripts/sbatch/*.sh` wrappers (see
   `step1.md` for the full procedure).

6. Rebuild the coordinate + sequence catalogs:
   ```
   sbatch scripts/sbatch/coordinate_catalog.sh
   sbatch scripts/sbatch/finalize_coordinate_catalog.sh
   sbatch scripts/sbatch/sequence_catalog.sh
   sbatch scripts/sbatch/export_pairs.sh
   ```

Every entry point enforces a `SLURM_JOB_ID` gate; direct invocation on a
login node will exit with code 78.
"""
    (RELEASE / "REPRODUCING.md").write_text(txt)


def main() -> int:
    require_slurm()
    t0 = time.time()

    RELEASE.mkdir(parents=True, exist_ok=True)
    _log(f"RELEASE = {RELEASE}")

    # ---- stage code ----
    _log("staging code tree …")
    for item in CODE_WANTED:
        src = CODE_ROOT / item
        if not src.exists():
            _log(f"  skip (missing): {item}")
            continue
        dst = RELEASE / item
        if src.is_dir():
            _copytree_filtered(src, dst)
        else:
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
        _log(f"  staged {item}")

    # ---- stage data deliverables ----
    _log("staging data deliverables from /groups …")
    for src_rel, dst_rel in DATA_WANTED:
        src = DATA_ROOT / src_rel
        if not src.exists():
            _log(f"  skip (missing): {src_rel}")
            continue
        dst = RELEASE / dst_rel
        if src.is_dir():
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(src, dst, symlinks=False)
        else:
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
        _log(f"  staged {dst_rel}")

    # ---- add README / LICENSE / REPRODUCING / .gitignore ----
    _log("writing README / LICENSE / REPRODUCING / .gitignore …")
    write_readme()
    write_license()
    write_reproducing()
    write_gitignore()

    # ---- git init + commit ----
    _log("git init + initial commit …")
    git_dir = RELEASE / ".git"
    if git_dir.exists():
        shutil.rmtree(git_dir)
    _run(["git", "init", "-b", "main"], cwd=RELEASE)
    _run(["git", "config", "user.email", "noreply@example.com"], cwd=RELEASE)
    _run(["git", "config", "user.name", "CRISPR Pilot Release"], cwd=RELEASE)
    _run(["git", "add", "-A"], cwd=RELEASE)
    _run(["git", "commit", "-m",
          "Initial release — frozen 1000-genome CRISPR array / mini-array "
          "coordinate + sequence catalog"], cwd=RELEASE)
    sha = _run(["git", "rev-parse", "HEAD"], cwd=RELEASE).stdout.strip()
    size_out = _run(["du", "-sh", str(RELEASE)]).stdout.strip()
    log_summary = {
        "release_root": str(RELEASE),
        "initial_commit": sha,
        "total_size_du": size_out,
        "elapsed_s": round(time.time() - t0, 1),
    }
    (RELEASE / "RELEASE_INFO.json").write_text(
        json.dumps(log_summary, indent=2) + "\n")
    _log(f"[done] commit {sha[:12]} ; {size_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
