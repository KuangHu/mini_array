# Step 2 v2: corrected pairs and separate candidate records

Status: implemented and tested locally; the real 1000-genome dataset has
**not** been rerun. Staged genomes and raw BLAST HSP files are not included
in this repository. Historical results remain unchanged. No corrected
biological candidate count is claimed by this change.

## Classification and evidence

Each pair must use the **same query ID**, same contig, same BLAST alignment
strand, and an aligned-boundary spacer of 15–60 bp. We deliberately use a
conservative common-model requirement instead of inventing a new threshold
for cross-model compatibility. Distinct query attributions to the same
genomic locus are retained separately; incompatible left-A/right-B hits
alone cannot support a pair. This does not establish biological function.

At least one repeat must have query coverage ≥0.90 and identity ≥0.90.
It may occur on either genomic side, on either alignment strand. Each
secondary repeat must have coverage ≥0.60 and identity ≥0.70.

Coverage is `(abs(qend - qstart) + 1) / qlen`, not BLAST alignment-column
length divided by query length. Raw query endpoints, query length and gap
open count are retained, with unrounded coverage. This removes inflation
from query-gap columns; it does not claim that a gapped alignment is an
exact contiguous copy of the query.

| Pair kind | Length completeness | Sequence identity |
| --- | --- | --- |
| `mini_array` | Both repeats cover ≥0.90 of the common query | At least one ≥0.90; other ≥0.70 |
| `half_repeat` | One full anchor; other coverage ≥0.60 and <0.90 | Short copy ≥0.70 |

`half_repeat` is an operational name for a shortened copy, **not exactly
50% length**. Copies below 60% coverage remain outside this search.
Each copy has independent `completeness=full/truncated` and
`sequence_quality=conserved/degenerate` fields. A full-coverage copy with
75% identity is full and degenerate, not truncated.

In exports, `repeat_1` always means genomic-left and `repeat_2` means
genomic-right. Either can be the shortened copy. All exported DNA is sliced
from the staged forward genomic sequence, including minus-strand hits.
BED strand is `.`; BLAST alignment orientation is not transcription direction.

## Contig edges and counting

- The initial pair screen excludes either flank <100 bp before locus
  deduplication and accepted pair/background counts. Such HSP pairs go to
  `pairs_*.edge_excluded.tsv` only.
- Final structural review requires the existing full **2,000 bp on each
  side**. Missing flank takes priority even when a known array overlaps.
  These loci are `edge_excluded`, absent from accepted candidate views,
  candidate FASTAs, protein neighborhoods and candidate cards.
- Initial `all_pairs.tsv` is the 100 bp screened **search evidence** table,
  including ordinary arrays and loci subsequently excluded by structural
  review. It is not the final candidate list. Use `candidate_loci.tsv` and
  the mini-array/half-repeat views for final candidates.
- One physical locus can have both kinds of pair evidence. Its pair rows
  occur in the corresponding separate views. Sum of the two view counts
  can exceed the physical-locus union; `summary.json` reports the overlap.

Structural unit checks now count repeat units, not spacer units. A lone
repeat in a low-confidence call is no longer sufficient by itself to label
an ordinary array fragment; the existing bona-fide unit evidence rule and
the ≥2 repeat rule remain. CRISPRCasTyper evidence is still array-level;
this revision does not provide missing per-copy CCTyper coordinates.

## Run on the original staged inputs

Use the locked Step 1 inputs and environments in `REPRODUCING.md`.
Deploy this checkout to the existing cluster code root. Configure all v2
paths in `config/step2.yaml`; default directories are `step2_v2`,
`step2b_v2`, `step2c_v2`, and `candidate_catalog_v2` under `results/`.
Use fresh Step 2 and catalog directories for each run; the program refuses
nonempty output directories to prevent mixing old controls or exports.
`PROJECT_ROOT` can override the data root; `CONFIG` can select another
config in the wrappers. The wrappers retain the existing cluster code path.

Submit these four jobs in order, waiting for successful completion before
each next job, or use SLURM `afterok` dependencies:

```bash
sbatch scripts/sbatch/step2_pilot.sh
sbatch scripts/sbatch/step2b_explainer.sh
sbatch scripts/sbatch/step2c_prepare.sh
sbatch scripts/sbatch/candidate_catalog_v2.sh
```

The retained `step2_pilot.sh` filename now launches v2 via the config.
Step 2C preparation performs physical deduplication and structural review;
protein prediction is not required for these catalogs. Optional protein
analysis is separate and must use the same v2 config.

To reprocess available raw BLAST HSPs without running BLAST again:

```bash
HITS_DIR=/path/to/original/step2_pilot sbatch scripts/sbatch/step2_pilot.sh
```

The raw directory must contain `hits_real.tsv` and, with the default
config, `hits_bg_rep0.tsv` through `hits_bg_rep2.tsv`. The 15-column outfmt
is specified in `config/step2.yaml`. Copy background-query FASTAs alongside
them when available so the shuffle audit can also run. Missing requested
HSP controls fail rather than silently becoming zero counts. A deliberately
real-only run can use `step2.py run --no-background` inside a SLURM job.
Raw-hit replay inherits the original BLAST search's recall and reporting
limits; it corrects interpretation of reported HSPs, not missing HSPs.
Do not replay old `pairs_real.tsv`: it cannot recover omitted left-partial
matches or reliably repair gap-inflated coverage.

Run manifests record method version, resolved configuration, input/code
hashes and raw HSP hashes. Review and catalog outputs verify upstream hashes.
An old rounded pair export without v2 provenance is rejected.

## Current outputs

`results/step2_v2/` contains accepted `pairs_real.tsv`, separate
`pairs_real.mini_array.tsv` and `pairs_real.half_repeat.tsv`, edge audit
records, locus attribution tables, and equivalent background outputs.

`results/candidate_catalog_v2/` contains:

- `candidate_loci.tsv`: unique accepted physical loci, with evidence kinds
  and counts for each kind.
- `candidate_pairs.tsv`: pair-level support linked to physical locus IDs.
- `mini_arrays.tsv` and `half_repeats.tsv`: separate physical-locus views.
- `mini_arrays.pairs.tsv` and `half_repeats.pairs.tsv`: separate detailed
  pair records, including both repeat completeness/identity attributes.
- Corresponding `.bed`, `.fna`, and `.index.tsv` files. Each kind's FASTA
  includes physical-locus, full-pair and three component records;
  `record_type` in the index distinguishes these. IDs and sequence hashes
  are checked in tests. No sequence deduplication removes correspondence.
- `all_pairs.tsv/.fna/.index.tsv`: initial search evidence, not final calls.
- `edge_excluded.pairs.tsv` and `edge_excluded.review.tsv`: audit only;
  no accepted candidate sequences for these records.
- `summary.json`: dynamic counts, physical-locus overlap, file hashes.

The frozen all-array coordinate/sequence catalogs still describe Step 1;
they are not rewritten by this candidate exporter. The old catalog-builder
scripts remain historical utilities and require explicit `--legacy`.
To reproduce the original release exactly, use its pinned commit.

## Validation

Local Windows/Python 3.13 validation: **55 passed, 12 skipped**. The skips
are 11 CRISPRidentify integration tests requiring the external tool/model
environment and one cluster-layout check requiring mounted cluster paths.
This is not a rerun of BLAST or the 1000-genome scientific results.

New tests exercise both anchor sides and both strands, reverse-complement
invariance, cross-query rejection, overlapping HSP alternatives, gap
coverage, identity/completeness separation, edge boundaries and precedence,
legacy-input rejection, real/background symmetry, missing controls, empty
results, mixed-kind loci, and genomic sequence/hash round trips through
synthetic replay → Step 2B → Step 2C → catalog export.

The optional no-edlib fallback now uses exact Levenshtein distance rather
than SequenceMatcher, fixing a pre-existing local family-clustering test
failure. Locked cluster environments with edlib keep the existing edlib
path; published Step 1 data are untouched.

```bash
python -m pytest tests -q
```
