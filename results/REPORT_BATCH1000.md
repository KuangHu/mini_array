# Step 1 — batch-1000 delivery report

Date: 2026-10-16
Approval: User review 2026-10-15 "按当前报告，受控 patched Step 1 pilot 通过",
scale-up approved up to 1000 inputs with explicit batch manifest.

**Scope:** 1000 GTDB r220 bacterial representatives (100 original pilot + 900
deterministic additions, same seed=42 round-robin across 196 phyla). Primary
threshold `seed_qc.min_score=0.75`; strict comparator `0.90`. All frozen
spec (software, resources, patches, models, env) per
`provenance/LOCKED.md`.

## tl;dr

- **1000/1000 identify + 1000/1000 cctyper success**, 0 execution failures.
- **239 HQ seeds** at primary 0.75 (213 NR100 repeats, 155 families).
- **75 HQ seeds** at strict 0.90 (65 NR100 repeats, 39 families, subset of
  0.75 → every family_id in 0.90 is identical to its 0.75 family_id).
- **652** coordinate-reliable bona-fide arrays in the known-array mask
  (preserves ALL coord-reliable bona-fide, not shrunk to HQ).
- 100 arrays failed the coordinate self-check (~0.14% of all candidates) —
  isolated upstream edge cases, flagged `invalid`, retained in audit tables.
- Patched classifier scored exactly as expected: isolate-heavy phyla produce
  higher HQ yield than fragmented MAGs; this is captured in stratification
  table, not treated as a bug.

---

## (a) Library sizes at both thresholds

| level | primary (0.75) | strict (0.90) |
|---|---:|---:|
| candidates_total | 69278 | 69278 |
| bona_fide_coord_pass | 652 | 652 |
| **hq_arrays** | **239** | **75** |
| reserve_arrays | 1318 | 1482 |
| context_only_arrays | 67621 | 67621 |
| invalid_arrays (coord-fail) | 100 | 100 |
| **nr100_repeats_hq** (search library) | **213** | **65** |
| **families_hq** | **155** | **39** |

Family IDs are stable across thresholds (0.90 ⊂ 0.75). The main search
library is `repeats.hq.nr100.fna` — **all** NR100 HQ consensuses are kept,
not just one per family. `known_array_regions.bed` has 652 rows in both
exports (defined as bona-fide + coord-check-pass, independent of HQ tier).

## (b) Filter loss bona-fide → HQ (primary 0.75)

- Bona-fide candidates: **658**
- Coord-check pass: **652** (6 fail = 0.9% upstream edge case)
- HQ at primary: **239**
- Reserve at primary: **413**

| gate that disqualified HQ | count |
|---|---:|
| `flank<100bp` | 401 |
| `ambiguous_bases_in_context` | 22 |
| `n_full_repeats<3` | 1 |

The dominant filter is **flank**: 97% of reserve bona-fides sit ≤100 bp from
a contig edge. This matches the earlier observation that bacterial MAGs
often fragment just outside array boundaries. Loosening the flank requirement
is out-of-scope for this round (user explicitly instructed "100 bp flank
不变").

## (c) Family cross-genome support (primary, 155 families)

The top 10 families appear in 2-9 genomes and 1-9 phyla each. These are
operational similarity families (90% norm. edit sim + 90% length ratio),
NOT evolutionary families.

Top 10 by #genomes supporting:

| family_id | n_repeats | n_arrays | n_genomes | n_phyla |
|---|---:|---:|---:|---:|
| fam_7a805f541ac1b6ee72cf | 7 | 9 | 9 | 9 |
| fam_6f393e61ef58ac27c10c | 7 | 11 | 7 | 5 |
| fam_32bff8536f7e4b0e82bd | 5 | 5 | 5 | 5 |
| fam_18747d6a513562247094 | 5 | 8 | 4 | 3 |
| fam_c746fdff9d0f000b08a8 | 3 | 4 | 4 | 3 |
| fam_88b45a39a9cd4a86b359 | 3 | 4 | 3 | 3 |
| fam_0cf6ecb62bb5aafb0abe | 2 | 2 | 2 | 2 |
| fam_12de37b63356e8a5bf00 | 2 | 2 | 2 | 2 |
| fam_1b23be31875941438c5a | 2 | 2 | 2 | 2 |
| fam_26dc06856c7f90ece5a6 | 2 | 3 | 2 | 1 |

Full per-family support table: `results/batch1000_delivery_metrics/c_family_support.tsv`.

## (d) Stratification

### by `source_type`

| source | n_genomes | candidates | bona-fide | hq (0.75) |
|---|---:|---:|---:|---:|
| MAG | 933 | 64681 | 570 | 165 |
| isolate | 67 | 4597 | **82** | **74** |

Isolate yield is dramatically higher per-genome:
`74 HQ / 67 isolates = 1.10 per genome` vs
`165 HQ / 933 MAGs = 0.18 per genome` — a **6× gap**, driven almost entirely
by MAG assembly fragmentation (next table).

### by assembly contiguity (contig count)

| bucket | n_genomes | candidates | bona-fide | hq (0.75) |
|---|---:|---:|---:|---:|
| 1-10 contigs | 63 | 6001 | 64 | **57** |
| 11-100 contigs | 285 | 17053 | 203 | 82 |
| 101-500 contigs | 567 | 37925 | 337 | 88 |
| >500 contigs | 85 | 8299 | 48 | 12 |

HQ-per-bona-fide drops from **89%** (1-10 contigs) to **25%** (>500 contigs)
as expected — fragmented assemblies truncate flank context.

### by phylum (top 25 by #genomes)

Full table in `results/batch1000_delivery_metrics/SUMMARY.md`. Notable
patterns:

- `Bacillota_B` yielded **15 HQ** from 7 genomes (highest density).
- `Bacillota_F` yielded 8 HQ from 7 genomes.
- Many phyla produced 0 HQ despite having bona-fide candidates (all failed
  flank or coord check) — e.g. `Bacteroidota` 0 HQ from 7 genomes with 0
  bona-fide after coord check.

## (e) Runtime, failure rate, disk

Measured over 1126 identify + 1310 cctyper SLURM jobs (pilot + batch1000 +
retries).

| tool | n_run | n_true_completed | n_duplicate_rejected | median_sec | mean_sec | max_sec | median_disk_bytes |
|---|---:|---:|---:|---:|---:|---:|---:|
| identify | 1126 | 1000 (dedup) | 0 | 96 | 130 | 1209 | 194,494 |
| cctyper | 1310 | 1000 (dedup) | 310 | 23 | 25 | 84 | 138,821 |

- `n_duplicate_rejected` are the cctyper safety-exit codes (OUT dir exists)
  — the sbatch wrapper refuses to overwrite existing output. Not real
  failures.
- 7 cctyper jobs exited `FAILED 1:0` on first attempt; they were
  automatically re-dispatched by the throttled submitter and succeeded.
- Zero `feature-extraction` failures (preflight validator catches
  resource-missing conditions before Pipeline starts).

Per-genome disk cost (median): identify ~190 KB, cctyper ~135 KB, staged
input ~2.9 MB. Combined ~3.3 MB per genome.

## (f) Full-DB (107,235 bac120 reps) extrapolation

Scaling per-genome medians:

| resource | per-genome | extrapolated to 107,235 |
|---|---:|---:|
| CPU time (identify + cctyper) | 476 s × 4 cores = 0.53 cpu-h | **14,180 cpu-h** |
| Wall time @ 30 concurrent jobs | — | **118 h ≈ 5 days** |
| identify output | 190 KB | 20.9 GB |
| cctyper output | 135 KB | 14.9 GB |
| staged input | 2.9 MB | 311 GB |
| **total disk** | — | **~346 GB** |

Assumptions:
- 4 cpu-per-job × 1 wall-second × job = 4 cpu-seconds.
- 30 jobs in flight is the sustained concurrency we observed under
  SLURM's current QOS limits (QOSMaxSubmitJobPerUserLimit=200, plus the
  partition scheduler effectively running ~10-30 jobs concurrently per
  user depending on load).
- No account for job queuing overhead — the 5-day estimate is wall-time
  assuming saturated flow.

## Non-negotiables honored this round

- Native primary `Bona-fide` required for HQ (Possible stays in reserve).
- 100-bp flank kept.
- Cas absence NOT disqualifying.
- CRISPRCasTyper disagreement NOT disqualifying.
- Query-seeds (NR100 HQ) and known-array mask (ALL bona-fide coord-reliable)
  kept as separate sets.
- Alternative / Possible Discarded / Low-score kept as soft context
  (`context_only`).
- `allow_full_run: false` preserved.
- Step 2 NOT started.
- Upstream issue NOT auto-submitted (draft at
  `provenance/upstream_issue_draft.md`).

## Delivery inventory

```
results/
  batch1000_delivery_metrics/            (all section-wise JSONs + SUMMARY.md)
    a_libraries.json
    b_filter_loss.json
    c_family_support.json / .tsv
    d_stratification.json
    e_runtime_disk.json
    f_full_db_extrapolation.json
    SUMMARY.md
  REPORT_BATCH1000.md                    (this file)
  REPORT_PATCHED.md                      (prior 100-genome report)
  step1_patched_min_score_0.75/          primary standard — recommended seed set
    arrays.{all,hq,reserve}.tsv
    array_units.tsv
    array_cas_links.tsv
    identify_vs_cctyper.tsv
    known_array_regions.bed              652 rows
    run_status.tsv
    repeats.hq.nr100.fna                 213 sequences
    repeat_sources.tsv
    repeat_family_members.tsv
    repeat_families.tsv                  155 families
    report.md
  step1_patched_min_score_0.90/          strict comparator subset
    (same files; 65 NR100, 39 families; 75 HQ arrays)
  step1_broken_feature_extraction_min_score_{0.75,0.90}/
    (preserved as evidence with BROKEN.md annotations)

provenance/
  LOCKED.md                              frozen spec
  CRISPRidentify.commit / CRISPRCasTyper.commit / CRISPRloci.commit
  models_tandem.hmm.{resource.json, hmmstat.txt, profiles.tsv}
  patches/                               patch diffs + SHA256SUMS
  fingerprint.json                       cache key
  env lockfiles + help text + interface notes
  upstream_issue_draft.md                NOT submitted
  pilot_manifest.build.json              1000-genome selection record
  pilot_manifest.build.pilot100.json     original 100-genome snapshot
```

## What this round does NOT authorize

- Full-DB (107k+) runs. The 118 h / 346 GB extrapolation is for planning
  only; sign-off required before any full run.
- Threshold or feature-definition adjustments.
- Loosening the 100-bp flank rule.
- Starting Step 2.
- Submitting the upstream issue draft.

## Open items for the next review

- Decide whether the 401 bona-fide calls lost to the flank rule should
  be revisited — options include a 50-bp "secondary HQ" tier (requires
  explicit approval), or waiting until Step 2's full-genome spacer walk
  provides independent cross-contig support.
- The isolate / MAG recovery gap (6×) is empirical, not a bug. Pilot 2
  could sample isolates more heavily if query diversity for the HQ seed
  library is a priority.
- CRISPRloci-sourced `models_tandem.hmm` has 79 profiles; no profile was
  observed to drive any array into Bona-fide on its own in this batch
  (`hmmr_score` contributes but doesn't dominate score for most hits).
  A per-profile hit-rate audit could inform whether a lighter HMM would
  speed up at scale.
