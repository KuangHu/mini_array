# Interface notes — Step 1 upstream software

Observations from the actually installed commits. Update as behavior differs from upstream docs.

## CRISPRidentify @ 474896ce3c9b89a609e69ac9d2fcdf78e4e2aad4

- CLI help captured in `CRISPRidentify.help.txt` (77 lines).
- `--json_report` expects a directory path, not a boolean. Confirmed in help: `--json_report JSON_REPORT`.
- Argparse runs before TensorFlow import, so `--help` works on nodes without AVX. Full runs require AVX.
- Default behaviour caps that step1.yaml does not override and we should be aware of:
  - `--max_identical_spacers` default `4`
  - `--max_identical_cluster_spacers` default `3`
  - `--enhancement_max_min` and `--enhancement_start_end` default `True`. Per step1.md §3, passing the string `False` through argparse `type=bool` does not reliably disable these — do not try to turn them off in the config.
- `--cas` must be `False` (step1.md §2) — we handle Cas annotation via CRISPRCasTyper separately.

## CRISPRCasTyper @ a2d51831d4802533bb113bac23cf9df2853063b0

- `cctyper 1.9.0`. CLI help captured in `cctyper.help.txt`.
- Positional args: `cctyper <input.fasta> <output_prefix>`.
- `--prodigal` takes `single|meta` (default `single`). Project default per step1.yaml section 8: `single` for isolates, `meta` for MAG/phage/plasmid/unknown.
- `--ccd 10000` for operon detection window (step1.md §8).
- `--skip_blast` disables the known-repeat BLAST search; step1 does not use that path.
- `--no_plot` to skip the plot step on headless/batch runs.
- `--repeat_search_backend` auto-picks between myers/edlib/blast — upstream behavior; leave at auto.

## Node-local quirks

- Current Claude-side execution node reports `__archspec=westmere` (pre-AVX). TensorFlow 2.3.0 wheels in `envs/identify` crash with `Illegal instruction (core dumped)` on import. CRISPRidentify runs must therefore be submitted via SLURM to a compute node with AVX2 (standard Rubin cluster nodes).
- `sbatch`/`sinfo` fail from this node with `resolve_ctls_from_dns_srv: res_nsearch error: Unknown host`. The user submits sbatch wrappers from a login node.
- Conda shared pkg cache at `/groups/rubin/projects/kuang/trash/conda-pkgs` has intermittent NFS stale-handle errors. The identify env was built with `CONDA_PKGS_DIRS=$PROJECT_ROOT/.conda_cache/pkgs` to avoid this. Keep that override for future rebuilds.
