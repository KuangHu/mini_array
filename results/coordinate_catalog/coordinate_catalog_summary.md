# Coordinate catalog — final export semantics

Pure coordinate catalog. No biology. No new DB. No parameter change. 0-based, half-open on BED. All analysis ran via sbatch; a SLURM_JOB_ID gate in every entry point refuses login-node execution.

## 1. CRISPR arrays — ALL tool calls (not independent trusted loci)

`crispr_arrays.{tsv,bed}` is the union of every call produced by either CRISPRidentify or CRISPRCasTyper on the frozen 1000-genome batch. **This is NOT a set of independent trusted arrays** — it includes Possible/Low score/Alternative/Possible Discarded identify calls AND Untrusted cctyper calls. Each row carries `source_tool` + `source_id` so attributions are preserved. Rows from the two tools MAY overlap on the same genomic interval; **adding counts across tools would double-count overlapping loci**.

- Total rows: **70207**
  - `crispridentify`: 69178
  - `crisprcastyper`: 1029
- Identify / cctyper call-class distribution:
  - `Low score`: 66657
  - `Possible`: 905
  - `Bona-fide`: 652
  - `Possible Discarded`: 600
  - `Trusted`: 553
  - `Untrusted`: 476
  - `Alternative`: 364

## 2. Trusted-call view

`crispr_arrays_trusted.{tsv,bed}` filters the ALL-tool-call table to rows that are **coordinate-valid AND** one of:
- `source_tool=crispridentify` AND `category_raw=Bona-fide`, OR
- `source_tool=crisprcastyper` AND `score_or_trusted=True`.

This is a per-call view. We DO NOT merge CRISPRidentify Bona-fide + cctyper Trusted calls that fall on the same genomic interval, and we DO NOT report an independent locus count by summing tool rows. Downstream work should dedup per (contig_uid, start0, end0) if a per-locus count is needed.

- Trusted-call rows: **1205**
  - `crispridentify`: 652
  - `crisprcastyper`: 553
- coordinate_in_bounds=true (by filter definition): 1205

## 3. Coordinate validation (split)

Validation is now reported in two independent columns:
- `coordinate_in_bounds` (true/false): did [start0, end0) fit inside the staged contig?
- `has_ambiguous_bases` (true/false): does the slice contain any non-ACGT base? Ambiguity is NOT reported as 'out of bounds'.

### CRISPR arrays (all tool calls)
- coordinate_in_bounds=true: 70207
- coordinate_in_bounds=false: 0
- has_ambiguous_bases=true: 1155
- has_ambiguous_bases=false: 69052

### Mini-array + edge candidates
- coordinate_in_bounds=true: 51
- has_ambiguous_bases=true: 0

## 4. Mini-array semantics

- `candidate_mini_array` physical loci: **40**
- `edge_unresolved` physical loci (separate file): **11**
- The mini_arrays.bed `strand` column is forced to `.` — the underlying hypothesis comes from a BLAST query alignment and is NOT an independent transcription-direction call. The TSV preserves `strand_hypotheses` AND a `strand_source` column that identifies the raw source (step2_pair_alignment). Treat the strand only as a hint for aligning queries, never as a biological orientation.

### mini + edge strand-hypothesis distribution
| strand_hypotheses | count |
|---|---:|
| `plus` | 22 |
| `minus` | 20 |
| `minus|plus` | 9 |

### The '78 both-strand loci' statement — scope

- `78` refers to the full set of **259 physical loci** (all structural classes) in `step2c_pilot/physical_candidate_evidence.tsv` where the `strand_hypotheses` column is `minus|plus` (both strand-pair reps collapsed under one physical_locus_id by Step 2C dedup).
- Within the 40 `candidate_mini_array` rows kept in `mini_arrays.tsv`: **8** have `strand_hypotheses=minus|plus`.
- Within the 11 `edge_unresolved` rows: **1** have both strands.
- Within the remaining 209 physical loci (ordinary_array_fragment / terminal_extension / repeat_like_ambiguous — all classes excluded from the mini-array file): **69**.

## 5. mini_array_units.tsv — one row per pair

- Rows: **70** (= pairs), covering **51** physical loci.
- New columns: `pair_id`, `physical_locus_id`, `locus_class` (= the mini/edge candidate's structural_class).
- Pair-level breakdown by locus_class:
  - `candidate_mini_array`: 57
  - `edge_unresolved`: 13
- Pair-level repeat_2 completeness distribution:
  - `full`: 43
  - `partial`: 27
- Pair-level strand (from Step 2 BLAST pair, NOT a biological call):
  - `plus`: 40
  - `minus`: 30

### Caveat: degenerate=0 is a secondary-threshold artefact

The Step 2 pair finder requires each secondary (second-repeat) hit to pass `identity ≥ 0.70 AND coverage ≥ 0.60`. Hits that would have been labelled `degenerate` under a looser rule were filtered out BEFORE pair formation. Therefore `repeat_2_completeness='degenerate' = 0` in this catalog is a **consequence of the frozen secondary thresholds**, not a statement that biologically degenerate second repeats do not exist. Changing the thresholds is explicitly out of scope for this round.

## 6. Non-negotiables honoured

- No new search, no new biology, no CDS, no protein, no HMM.
- Candidate boundaries, classification thresholds, search parameters: unchanged.
- `allow_full_run: false` preserved. No new DB. No upstream issue submission.
- All analysis sbatch-only (SLURM gate enforced at import).

