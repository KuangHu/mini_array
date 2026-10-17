# Step 2 pilot report — anchor-pair repeat search

Date: 2026-10-16
Approval ref: User review 2026-10-16, "批准有限范围的 Step 2 pilot"
Scope lock: Only the 1000-genome frozen batch + 213 primary HQ queries +
65 strict (0.90) support queries. No new DB download. No full-DB run.
Step 2 pipeline, caller version pin, config and non-negotiables per the
2026-10-16 instruction.

## tl;dr

- **1,073 unique candidate loci across 247 genomes** after dedup.
- Bucket distribution: 735 inside known arrays, 260 overlap a Possible
  soft-context region, **73 standalone, 5 at contig edge**.
- **531 loci (142 genomes) have at least one strict-support query**
  (query sourced from the 0.90 HQ library).
- **Background enrichment: real 15,472 pairs vs 37.7 pairs mean of 3
  dinuc-shuffled replicates ≈ 410× over shuffle.** This is enrichment,
  not a calibrated FDR.
- One of the 3 background replicates (rep1) returned 113 spurious pairs —
  single-replicate noise; mean + range is reported honestly.
- 73 standalone candidate loci are the most interesting handoff — of
  these, 15 carry strict-support (0.90 query attribution).

## Pipeline parameters

Per `config/step2.yaml`:

| parameter | value | source |
|---|---:|---|
| anchor min_cov_query | 0.90 | fixed per user instruction |
| anchor min_identity | 0.90 | fixed per user instruction |
| second_repeat min_cov_query | 0.60 | truncated/degenerate copy allowed |
| second_repeat min_identity | 0.70 | degenerate copy allowed |
| spacer min_len | 15 nt | per user instruction |
| spacer max_len | 60 nt | per user instruction |
| overlap_merge_bp | 10 | duplicate-anchor merge |
| locus_dedup_merge_bp | 60 | dedup footprint merge |
| background method | dinuc_shuffle | AE 1985 edge-swap |
| background replicates | 3 | enrichment, not FDR |
| background seed | 42 | — |
| blastn word_size | 7 | — |
| blastn evalue | 10 | permissive; filtering done downstream |
| blastn dust | no | — |

All classification thresholds are applied in Python *after* a single
permissive BLAST run (one scan, both strands via BLAST's native handling);
the canonical direction of each query is **not** treated as a transcription
direction.

## Resource usage (this pilot — single sbatch job)

- Job 12402 wall time: **6:23**
- CPUs allocated: 8
- CPU time: 6:23 × 8 = 51.1 cpu-min = **0.85 cpu-h**
- Peak transient disk: 2.9 GB combined-FASTA + 740 MB blast DB + 13 MB
  real hits + ~0.7 MB per background hit file + 2.2 MB pair tables + 1.6 MB
  loci attribution = **≈ 3.7 GB during run**
- Steady-state after cleanup: ~ 20 MB (results/step2_pilot/ TSV + logs)
  plus the kept subject DB (~3.6 GB) if `--keep-blastdb`.
- All within the user-set caps (max_concurrent_jobs=4, cpus_per_job=8,
  mem_per_job=32 GB, tmp_disk_cap=50 GB, wall_time_h=6).

## Section 1 — Search output

| counter | value |
|---|---:|
| raw BLAST hits (real) | 126,351 |
| anchor hits (≥0.90 cov, ≥0.90 ident) | 41,434 |
| secondary-passing hits (≥0.60 cov, ≥0.70 ident) | 119,104 |
| anchor-second pairs (15-60 nt spacer) | **15,472** |
| unique loci after dedup | **1,073** |

Full raw BLAST output at
`results/step2_pilot/{hits_real.tsv, pairs_real.tsv, loci.tsv,
loci_attributions.tsv}`.

## Section 2 — Context classification (per user instruction)

| bucket | loci | genomes |
|---|---:|---:|
| `known_array_interior` (overlap ≥50% of known array) | 735 | 212 |
| `possible_soft_context` (overlap a non-bona-fide valid array) | 260 | 80 |
| `standalone` (none of above, no contig edge issue) | **73** | **28** |
| `contig_edge` (within 100 bp of contig end) | 5 | 3 |
| **total** | **1,073** | **247** |

Known-array mask uses **all 652** coord-reliable bona-fide arrays from the
primary 0.75 export (per the "don't shrink to HQ" rule). Soft context uses
the 67,621 coord-valid non-bona-fide calls (Alternative/Possible/Possible
Discarded/Low score) from the same export.

## Section 3 — Query attribution preserved

Each locus retains the full set of queries that produced a pair at that
footprint. `loci.tsv` reports `query_set` and `n_unique_queries`; the
exploded `loci_attributions.tsv` has one row per `(locus, anchor_qseqid,
second_qseqid)` tuple with per-pair identity + spacer length. This lets
downstream analysis attribute enrichment to specific repeat sequences
without double-counting the shared locus.

**Strict (0.90) support:** 531 of 1073 loci (49.5%) carry at least one
query from the 65-sequence strict library.

| | all loci | with strict support |
|---|---:|---:|
| known_array_interior | 735 | 420 |
| possible_soft_context | 260 | 97 |
| standalone | 73 | 15 |
| contig_edge | 5 | 1 |

## Section 4 — Background comparison

Dinucleotide-preserving shuffle of each of the 213 real queries, 3
independent replicates (seed=42, Altschul-Erickson edge-swap).

| replicate | raw hits | anchors | pairs |
|---|---:|---:|---:|
| hits_bg_rep0 | 6,414 | 40 | 0 |
| hits_bg_rep1 | 6,275 | 158 | 113 |
| hits_bg_rep2 | 6,457 | 30 | 0 |
| **mean** | **6,382** | **76** | **37.7** |
| real | 126,351 | 41,434 | **15,472** |
| **real / bg_mean** | **19.8×** | **545×** | **~410×** |

Interpretation: strongly above background, enrichment holds across all 3
replicates for anchors. The single rep1 pair-count spike (113) illustrates
why n=3 replicates is **not** enough to calibrate FDR; it is sufficient only
to document that the signal is many orders of magnitude above random. Per
the user instruction, this is **not** reported as a calibrated FDR.

## Section 5 — Standalone candidate loci (potential new mini-arrays)

73 standalone candidates across 28 genomes. 15 carry strict (0.90) support.

First 10 strict-supported standalone loci:

| locus | genome | contig:start-end | strand | n_pairs | n_queries |
|---|---|---|---|---:|---:|
| loc_000256 | GCA_937139605.1 | c_3eda2bf23e69:24111-24190 | - | 1 | 1 |
| loc_000257 | GCA_937139605.1 | c_3eda2bf23e69:24111-24190 | + | 1 | 1 |
| loc_000344 | GCA_016214985.1 | c_579e14ded727:40156-40315 | + | 2 | 1 |
| loc_000353 | GCF_002742685.1 | c_5c345aa920cb:60866-60961 | + | 1 | 2 |
| loc_000400 | GCF_000020945.1 | c_6c5dd549e7ca:1276628-1276790 | + | 2 | 1 |
| loc_000438 | GCA_003574845.1 | c_796c69ab02b3:56961-57051 | - | 1 | 2 |
| loc_000439 | GCA_003574845.1 | c_796c69ab02b3:56961-57051 | + | 1 | 2 |
| loc_000451 | GCA_951828535.1 | c_7c479615e6c6:3135118-3135205 | + | 1 | 1 |
| loc_000452 | GCA_951828535.1 | c_7c479615e6c6:3135119-3135205 | - | 1 | 1 |
| loc_000537 | GCA_025774605.1 | c_91f3e150f49e:44390-44483 | + | 1 | 2 |

Full table: `results/step2_pilot/loci.tsv` (filter `context=standalone`).

### Queries driving the standalone candidates

| query repeat_id | pairs produced on standalone loci |
|---|---:|
| rep_b9539b1ec8a651f2b0bc | 167 |
| rep_3a57ed1fe934c58a6497 | 152 |
| rep_b145c10a1c64b946fcb7 | 136 |
| rep_ea7ceb4fb98e27b96111 | 124 |
| rep_2ce9e72294a010d2ee31 | 77 |

Several standalone candidates appear on both strands at nearly the same
coordinates (e.g. loc_000256 ↔ loc_000257, loc_000438 ↔ loc_000439,
loc_000451 ↔ loc_000452). The pipeline reports both strand orientations
without inferring which (if either) is transcriptional. Downstream review
should merge or keep both as needed.

## Non-negotiables honored

- Native primary HQ queries used (213 at 0.75), with 0.90 strict as a
  support flag, not a filter.
- `n_hq_arrays`, `n_nr100_repeats`, `n_families` reported separately
  (239 / 213 / 155 at primary).
- Query library kept intact — not reduced to one representative per family.
- All contigs scanned (211,631 contigs × 213 queries), not just seed-source
  genomes or Cas neighbourhoods.
- Both same-genome and cross-genome query↔locus attributions preserved.
- Full-length anchor + second repeat (complete/truncated/degenerate) with
  explicit coverage and identity thresholds.
- Spacer 15-60 nt measured between aligned repeat boundaries, same strand.
- Both strands scanned via BLAST's native handling; raw alignment coords
  preserved. Canonical direction NOT treated as transcription direction.
- Known-array mask uses all 652 coord-reliable bona-fide arrays — not
  reduced to HQ seed sources.
- Possible/Alternative/Low score kept as separate soft-context records.
- Composition-matched shuffled background included (not calibrated FDR).
- No Cas or DNA-modifying gene preset required for a candidate.
- No full-DB spacer-target search.
- `allow_full_run: false` preserved; this pilot ran as a single bounded
  sbatch job with explicit caps.
- No new DB downloaded.
- No upstream issue auto-submitted.
- Step 2 pilot delivery stops here.

## Delivery

```
results/step2_pilot/
  hits_real.tsv                  126,351 raw BLAST hits
  hits_bg_rep{0,1,2}.tsv         background hits per replicate
  pairs_real.tsv                 15,472 anchor-second pairs
  pairs_bg_rep{0,1,2}.tsv        background pairs per replicate
  loci.tsv                       1,073 dedup unique loci, with context + strict flag
  loci_attributions.tsv          every (locus, anchor_q, second_q) tuple
  summary.json                   initial counters (pre-genome-count fix)
  summary_v2.json                corrected counters including per-bucket genome counts
  bg_queries_rep{0,1,2}.fna      shuffled background queries used
  subject_db/                    makeblastdb output (kept for reuse)
  step2.log                      run log

config/step2.yaml                frozen parameters for this round
scripts/step2.py                 pipeline driver
scripts/step2_lib/               search, pair, context, dedup, dinuc_shuffle modules
scripts/sbatch/step2_pilot.sh    sbatch wrapper used (jobid 12402)
```

## Open items (not acted on — delivery stops here)

- Several standalone candidates lie on contigs that already contain many
  self-similar query hits (e.g. `GCA_958295815.1 c_e8c2c4e07de8` has 15
  standalone loci). Review whether these are repetitive non-array elements
  or genuine clusters.
- rep1 of the dinuc shuffle produced 113 pairs vs 0 for rep0/rep2 — a
  single-shuffle artefact; worth noting but doesn't change the enrichment
  conclusion. More replicates would narrow the mean.
- No per-locus alignment extracts are generated in this report. If the
  next review requires alignment images or multi-sequence alignments for
  the 15 strict-supported standalone candidates, that is a separate task.
