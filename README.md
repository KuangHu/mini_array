# CRISPR array / mini-array coordinate + sequence catalog

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

