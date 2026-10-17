# Step 2 pair-level sequence export (all 15,472 pairs)

Pure extraction. Every pair in `step2_pilot/pairs_real.tsv` has its full span + three component sequences dumped from the staged FASTA (genomic-forward DNA). Pair-level correspondence is preserved: identical component sequences across pairs are NOT deduped; each pair keeps its own records.

## Query library — origin confirmation

All anchors / seconds come from the 213 NR100 HQ repeat sequences in `results/step1_patched_min_score_0.75/repeats.hq.nr100.fna`. That library was built from the 652 Bona-fide + coordinate-pass identify arrays on the same frozen 1000-genome batch, deduped (exact + RC). No external repeat library is used.

## Record counts

- Pairs exported: **15472** (of 15,472 raw pairs; 0 dropped because contig/sequence lookup failed)
- Full-span FASTA: `all_pairs.fna`  (15472 records)
- Unit FASTA: `all_pair_units.fna` (46416 records = pairs × 3 unit types)
- Index: `all_pairs_index.tsv` (15472 rows × 28 columns)

## Step 2 context bucket distribution

| bucket | n_pairs |
|---|---:|
| `known_array_interior` | 14311 |
| `possible_soft_context` | 969 |
| `standalone` | 139 |
| `suspected_terminal_extension` | 47 |
| `contig_edge` | 6 |

## FASTA ID conventions

- `all_pairs.fna`: `step2_pair_<NNNNNN>`
- `all_pair_units.fna`: `step2_pair_<NNNNNN>|{repeat_1,spacer,repeat_2}`
- `repeat_1` is always the UPSTREAM repeat of the pair (lower start0). The `repeat_*_role_in_pair` column in the index records which was the original BLAST anchor vs secondary.

## Non-negotiables

- Only source of sequence: `input/staged/<genome_uid>.fna`.
- All sequences genomic-forward; no RC applied.
- Ambiguous bases preserved as-is.
- No deduplication (pair-level correspondence intact).
- No new search, no boundary / threshold / parameter change.
- `allow_full_run: false` preserved. sbatch-only.
