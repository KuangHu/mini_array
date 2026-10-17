# Sequence catalog (frozen)

Pure genomic-forward DNA export from the frozen coordinate catalog. No new biology, no new search, no re-identification. All sequences are literal slices of `input/staged/<genome_uid>.fna` using the 0-based half-open coordinates in `results/coordinate_catalog/`. Ambiguous bases (N and other non-ACGT) are preserved as-is.

## Scope

- This catalog covers ONLY the frozen 1000-genome batch and the calls that were already recovered under the current coordinate catalog (allow_full_run: false).
- `crispr_arrays.all_calls.fna` is the union of every call from CRISPRidentify + CRISPRCasTyper. **It is NOT a set of independent trusted arrays**; rows from the two tools may overlap on the same genomic interval.
- `crispr_arrays.trusted_calls.fna` applies the coordinate_catalog's trusted filter (identify Bona-fide OR cctyper Trusted) and is still per-call. Downstream work must dedup if an independent-locus count is needed.
- 210 unit records = 70 pairs × 3 unit_types (repeat_1, spacer, repeat_2), spanning both `candidate_mini_array` and `edge_unresolved` loci. The `locus_class` column in `sequence_index.tsv` partitions them.
- The existing degenerate-repeat coverage limitation from Step 2's secondary thresholds remains unchanged (see the coordinate_catalog summary for the note).

## Expected vs actual record counts (post-write validation)

| fasta_file | expected | actual |
|---|---:|---:|
| `crispr_arrays.all_calls.fna` | 70207 | 70207 |
| `crispr_arrays.trusted_calls.fna` | 1205 | 1205 |
| `mini_arrays.fna` | 40 | 40 |
| `edge_unresolved_candidates.fna` | 11 | 11 |
| `mini_array_units.fna` | 210 | 210 |

### Validation results

- missing contig lookups: 0
- length mismatches (len != end0 - start0): 0
- sequence sha256 mismatches: 0
- all per-file counts match expected values

## FASTA ID conventions

| file | ID format |
|---|---|
| crispr_arrays.all_calls.fna | `<source_tool>|<source_id>` |
| crispr_arrays.trusted_calls.fna | `<source_tool>|<source_id>` |
| mini_arrays.fna | `<physical_locus_id>` |
| edge_unresolved_candidates.fna | `<physical_locus_id>` |
| mini_array_units.fna | `<pair_id>|<unit_type>` where unit_type ∈ {repeat_1, spacer, repeat_2} |

Every ID is unique within its own FASTA file. `sequence_index.tsv` maps every ID back to its row in the frozen coordinate catalog.

## Deliverables

```
results/sequence_catalog/
  crispr_arrays.all_calls.fna
  crispr_arrays.trusted_calls.fna
  mini_arrays.fna
  edge_unresolved_candidates.fna
  mini_array_units.fna
  array_units_patched.tsv            # identify repeat/spacer units (copy of patched output)
  array_units_cctyper_unavailable.tsv # cctyper calls: 'unit annotation unavailable'
  sequence_index.tsv
  input_fingerprints.json            # coordinate catalog + staged FASTA fingerprints
  README.md
```

## Non-negotiables honoured

- Sole source of sequence: `input/staged/<genome_uid>.fna`. No consensus/query substitution. No re-identification of units.
- All sequence output is genomic-forward DNA; no automatic reverse-complement based on strand hypothesis. Strand hypotheses stay in `coordinate_catalog/mini_arrays.tsv`.
- Ambiguous bases preserved as-is; `has_ambiguous_bases` column in `sequence_index.tsv` flags any non-ACGT content.
- Count levels preserved. Identical unit sequences are NOT deduped — pair-level correspondence must be kept.
- No parameter / threshold / boundary change. No new search. No CDS / protein / HMM analysis.
- All execution via sbatch; SLURM_JOB_ID gate enforced at import.
- `allow_full_run: false` preserved.

## Input fingerprints

```json
{
  "input_coordinate_catalog": {
    "crispr_arrays.tsv": {
      "path": "/groups/rubin/projects/kuang/out/CRISPR_array_related/results/coordinate_catalog/crispr_arrays.tsv",
      "sha256": "c8abf7fdd01325bc6fa7b3eb3624e0d1fe246dfd1037dfd030c9e56af419028c",
      "size": 22047322
    },
    "crispr_arrays_trusted.tsv": {
      "path": "/groups/rubin/projects/kuang/out/CRISPR_array_related/results/coordinate_catalog/crispr_arrays_trusted.tsv",
      "sha256": "8de38b133fe834b381b61192f897efb76b4ebb363bd80535d3306dd993a43d33",
      "size": 304556
    },
    "mini_arrays.tsv": {
      "path": "/groups/rubin/projects/kuang/out/CRISPR_array_related/results/coordinate_catalog/mini_arrays.tsv",
      "sha256": "caa77008f5c3215bac9debf43416da24ceec8a9775c64f4f557ca3a8d9757948",
      "size": 13244
    },
    "edge_unresolved_candidates.tsv": {
      "path": "/groups/rubin/projects/kuang/out/CRISPR_array_related/results/coordinate_catalog/edge_unresolved_candidates.tsv",
      "sha256": "12463054ddf8a159b481d80202785f1ee49ad5546235fe03fceee3c53c2991c9",
      "size": 3633
    },
    "mini_array_units.tsv": {
      "path": "/groups/rubin/projects/kuang/out/CRISPR_array_related/results/coordinate_catalog/mini_array_units.tsv",
      "sha256": "7fa5ba65b94df290d1f510c11a65f3bf87a89c9f22013b2099c7acd4820262b6",
      "size": 15531
    },
    "coordinate_catalog_summary.md": {
      "path": "/groups/rubin/projects/kuang/out/CRISPR_array_related/results/coordinate_catalog/coordinate_catalog_summary.md",
      "sha256": "018ad6afd29908d5e8ea95b027529ce2046e9a197c01978694860dafe8fa0921",
      "size": 5073
    }
  },
  "input_staged_fastas": {
    "n_files": 1000,
    "listing_sha256": "f754e2b9ef4b65a35d8249447bb879f1e1901228ae52b782e2357cdd7a1b9822",
    "note": "Listing-level fingerprint (filename:size) over all 1000 staged FASTAs. Individual per-file hashes omitted for speed; each staged FASTA is deterministic from Step 1's staging code (see provenance/LOCKED.md)."
  }
}
```
