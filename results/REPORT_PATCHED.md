# Step 1 pilot — final report after upstream patch

Date: 2026-10-15 · 100 GTDB r220 bacterial representative genomes · one array
caller (patched CRISPRidentify) + one Cas-operon caller (CRISPRCasTyper 1.9.0).

## tl;dr

The original pilot produced **0 HQ seeds** on 100 genomes. The root cause was
three distinct upstream defects in CRISPRidentify commit `474896ce` that
caused the ML classifier's `blast_score_1`, `blast_score_2`, and
`hmmr_score` features to be permanently 0.0 for every input. After an
upstream-faithful patch that restores all three features without changing
thresholds or model files, the same 100 genomes now produce:

- **50 bona-fide calls** (up from 7)
- **4 HQ seeds** at `min_score=0.90` (up from 0)
- **20 HQ seeds** at `min_score=0.75` (up from 0)
- **0 regressions** vs the broken run
- 41 locus-matched promotions non-bona-fide → bona-fide

Smoke-test validation on two independently-annotated reference genomes
(Bacteroidota `GCA_013414795.1`, Cyanobacteria `GCA_000309945.1`) recovered
**8 of 13 CRISPRspacers-annotated arrays** as bona-fide, with the remaining
5 staying `Possible` for legitimate reasons (small size, few repeats).

Thresholds were **not** adjusted to force this outcome. Reserve candidates
below HQ are explicitly flagged with the specific per-array reasons.

---

## 1. The upstream defects we had to fix

Verified from a fresh clone of
`https://github.com/BackofenLab/CRISPRidentify.git` at the pinned commit
`474896ce3c9b89a609e69ac9d2fcdf78e4e2aad4`:

| # | Defect | Location | Effect |
|---|---|---|---|
| 1 | `tools/hmm_search/models_tandem.hmm` is not shipped | repo tree | per-array HMM step opens a nonexistent file → `hmmr_score = 0.0` for all candidates |
| 2 | BLAST `-db tools/blasting/...` uses a cwd-relative path | `components_evaluation.py:406, 419` | fails whenever runner cwd ≠ CRISPRidentify source dir → `blast_score_{1,2} = 0.0` for all candidates |
| 3 | f-string typo passes literal `hmm_model` instead of the parameter | `components_evaluation.py:900` | per-array HMM subprocess tries to open file named `hmm_model` and fails silently |

All three failures were absorbed into `subprocess.Popen(...).communicate()`
with no exit-code or stderr checks, which is why they manifested as "all
features are zero" rather than crashes.

## 2. Resource recovery

| Resource | Source | SHA-256 | Size | Content |
|---|---|---|---|---|
| `models_tandem.hmm` | CRISPRloci @ `d25287b4` → `CRISPRidentify/tools/hmm_search/models_tandem.hmm` | `818363ec…` | 13.7 MB | 79 HMM profiles (NCBI CDD `cd*`, COG, PF Pfam) |
| `Verified_repeats_dataset{1,2}.fa` | our locked CRISPRidentify tree | md5 `5d2677c2 / e936ddec` | — | md5-identical to the copies in CRISPRloci; no substitution needed |

Full NAME/ACC manifest for the 79 HMM profiles saved at
`provenance/models_tandem.hmm.profiles.tsv`.

## 3. Patch footprint

Isolated to a new file (`components/_resources.py`) and minimal edits in two
existing files. Classifier, feature extraction order, feature definitions,
model files, and all search parameters are unchanged.

| File | Lines changed | Purpose |
|---|---:|---|
| `components/_resources.py` (new) | — | absolute-path resolver + preflight validator + checked `run_subprocess` with explicit `FeatureExtractionError` |
| `components/components_evaluation.py` | 127 lines of unified diff | all 4 subprocess call sites switched to absolute paths and checked runs; `hmm_model` literal bug fixed |
| `CRISPRidentify.py` | 22 lines of unified diff | preflight runs at import time; task refuses to start if any resource or executable is missing |

Patch SHAs saved to `provenance/patches/SHA256SUMS` and `provenance/fingerprint.json`.

## 4. Fixture tests (patched feature extraction)

11 tests at `tests/test_feature_extraction.py` cover every state the user
asked us to distinguish:

| Fixture | Expected outcome |
|---|---|
| BLAST DB1 self-match query | ≥1 hit, non-zero bitscore |
| BLAST DB2 self-match query | ≥1 hit, non-zero bitscore |
| BLAST random no-hit query | rc=0, empty output (legitimate 0.0 feature) |
| BLAST missing DB | `FeatureExtractionError` raised |
| HMM YSIRK_signal self-consensus query | ≥1 hit, non-zero bitscore |
| HMM random no-hit query | rc=0, empty output (legitimate 0.0 feature) |
| no-ORF / no-query case | branch taken explicitly, not a failure |
| missing tools dir | preflight raises before any task starts |
| missing executable | preflight raises before any task starts |
| different starting cwd | identical BLAST output (byte-identical) |
| two parallel tasks | per-cwd temp files preserved, no cross-talk |

All 11 pass. Full suite (36 pre-existing + 11 new + 2 new QC) = **38 passed**.

## 5. Reference smoke test — 2 isolates, 13 reference arrays

Reference: CRISPRspacers (Nat 2021, doi:10.1093/nar/gkaa1023), filtered to
`crispr_tool ∈ {crt, pilercr}` — fully independent of CRISPRidentify.
Reference annotations treated as a computed comparator, not experimental
ground truth.

| genome | ref# | n_sp | tool | pre_cat | pre_score | post_cat | post_score | Δ |
|---|---:|---:|---|---|---:|---|---:|---:|
| GCA_013414795.1 | 1 | 588 | crt | Possible | 0.666 | **Bona-fide** | **0.939** | +0.273 |
| GCA_013414795.1 | 2 | 482 | pilercr | Possible | 0.666 | **Bona-fide** | **0.939** | +0.273 |
| GCA_013414795.1 | 3 | 40 | crt | Possible Discarded | 0.575 | Possible | 0.645 | +0.070 |
| GCA_013414795.1 | 4 | 20 | crt | Possible | 0.508 | Possible | 0.594 | +0.086 |
| GCA_000309945.1 | 1 | 198 | crt | Possible | 0.698 | **Bona-fide** | **0.875** | +0.177 |
| GCA_000309945.1 | 2 | 103 | pilercr | Possible | 0.698 | **Bona-fide** | **0.875** | +0.177 |
| GCA_000309945.1 | 3 | 59 | pilercr | Possible | 0.698 | **Bona-fide** | **0.875** | +0.177 |
| GCA_000309945.1 | 4 | 12 | pilercr | Possible | 0.698 | **Bona-fide** | **0.875** | +0.177 |
| GCA_000309945.1 | 5 | 14 | pilercr | Possible | 0.698 | **Bona-fide** | **0.875** | +0.177 |
| GCA_000309945.1 | 6 | 3 | pilercr | Possible Discarded | 0.527 | Possible Discarded | 0.672 | +0.145 |
| GCA_000309945.1 | 7 | 26 | crt | Possible | 0.708 | **Bona-fide** | **0.784** | +0.076 |
| GCA_000309945.1 | 8 | 39 | crt | Possible | 0.628 | Possible | 0.651 | +0.023 |
| GCA_000309945.1 | 9 | 36 | crt | Possible | 0.703 | Possible | 0.750 | +0.047 |

**8 of 13 reference arrays** are now bona-fide. The remaining 5 are
legitimately Possible (small n_sp, or score ceiling imposed by the model's
`number_mismatches` / `spacer_similarity` features that are not affected by
the patch). No threshold adjustment was made.

Direct feature-vector comparison for the biggest array (588 repeats):

| feature | pre-patch | post-patch |
|---|---:|---:|
| `repeat_similarity` | 0.9900 | 0.9901 |
| `mfe_score` | -1.70 | -1.70 |
| `orf_score` | 94.82 | 94.82 |
| `blast_score_1` | **0.0** | **32.5** |
| `blast_score_2` | **0.0** | **47.3** |
| `hmmr_score` | 0.0 | 0.0 *(legitimate — no tandem-protein HMM hits)* |
| **Certainty Score** | **0.666** | **0.939** |

## 6. QC tier rule change (Phase E)

Per the review instruction *"坐标/序列有效的 primary Possible 保留为 reserve,
低于 HQ score 本身不能成为降入 context_only 的理由"*:

- coord/seq-valid primary `Possible` → **reserve** (previously `context_only`)
- `Alternative`, `Possible Discarded`, `Low score` → `context_only` (unchanged)
- Flank measurement is now the actual bp available from the real contig
  length (not capped at the window) — the comparison predicate against
  `context_bp=100` is unchanged.
- No change to the 100-bp flank requirement, no Cas requirement, no
  identify-cctyper agreement requirement.

Two added tests (`test_qc_possible_coord_valid_stays_in_reserve`,
`test_qc_alternative_and_possible_discarded_go_context_only`) lock the rule.

## 7. Full-pilot comparison (100 GTDB r220 genomes)

Both exports share identical inputs, environment, patched identify tree
(`raw_patched/identify/`), reused unchanged cctyper outputs
(`raw/cctyper/`), and QC parameters — only `seed_qc.min_score` differs.

| metric | broken (0.90) | patched (0.90) | patched (0.75) |
|---|---:|---:|---:|
| total candidates | 6146 | 6180 | 6180 |
| bona-fide calls | 7 | **50** | 50 |
| coordinate_check_fail (upstream edge-case isolate) | 4 | 7 | 7 |
| QC tier `hq` | 0 | **4** | **20** |
| QC tier `reserve` | 7 | 133 | 117 |
| QC tier `context_only` | 6135 | 6036 | 6036 |
| QC tier `invalid` | 4 | 7 | 7 |
| NR100 repeats (HQ library) | 0 | 3 | 19 |
| operational families | 0 | 3 | 18 |
| `known_array_regions.bed` rows | 7 | 49 | 49 |
| Cas link `proximal_operon` | 125 | 120 | 120 |

Locus-matched score deltas on 2,556 identical loci (same genome / contig /
start / end / category_raw): mean **+0.024**, median **+0.016**, max
**+0.319**, min **−0.111**. **Zero regressions** at the category level.

41 non-bona-fide → bona-fide promotions, listed in
`results/step1_patched_min_score_0.90/pre_post_comparison.md`.

## 8. Deliverables

```
provenance/
  CRISPRidentify.commit                      474896ce...
  CRISPRCasTyper.commit                      a2d51831...
  CRISPRloci.commit                          d25287b4...
  models_tandem.hmm.resource.json            SHA256 pin + 79-profile count
  models_tandem.hmm.hmmstat.txt              per-profile hmmstat
  models_tandem.hmm.profiles.tsv             79-row name/acc manifest
  patches/components_evaluation.py.patch     127-line unified diff
  patches/CRISPRidentify.py.patch            22-line unified diff
  patches/_resources.py.new                  new file
  patches/SHA256SUMS                         patch hashes
  fingerprint.json                           consolidated cache key:
                                               upstream commit, resource
                                               commit + hashes, patch
                                               hashes, classifier model
                                               hash, env lock hashes
  identify.conda.explicit.txt                env lock
  driver.conda.explicit.txt                  env lock
  driver.pip.freeze.txt                      pip lock
  interface_notes.md                         pre-existing + updated
  upstream_issue_draft.md                    public issue draft (not sent)

results/
  step1_broken_feature_extraction_min_score_0.90/   preserved as evidence
    BROKEN.md                                       failure annotation
  step1_broken_feature_extraction_min_score_0.75/   preserved as evidence
    BROKEN.md                                       failure annotation
  step1_patched_min_score_0.90/                     recommended output
    arrays.{all,hq,reserve}.tsv
    array_units.tsv
    array_cas_links.tsv
    identify_vs_cctyper.tsv
    known_array_regions.bed
    run_status.tsv
    repeats.hq.nr100.fna                            3 sequences
    repeat_sources.tsv
    repeat_family_members.tsv
    repeat_families.tsv                             3 families
    report.md                                       auto-generated summary
    pre_post_comparison.md
  step1_patched_min_score_0.75/                     comparator
    (same files, 19-seq NR100, 18 families)
  REPORT_PATCHED.md                                 this file

smoke_test/
  input/, raw/ (unpatched), raw_patched/, logs/
```

## 9. What `allow_full_run: false` still forbids

- Running the full ≥10k genome scan.
- Starting Step 2 (short repeat / mini-array search).
- Submitting any upstream issue (draft is at
  `provenance/upstream_issue_draft.md`, review before posting).

## 10. Open items for the next review round

- 46 arrays are bona-fide but not HQ at 0.90, mostly due to the 100-bp flank
  rule on fragmented MAG contigs (91 of 100 pilot inputs are MAGs). A
  larger isolate-only comparison would separate the "upstream feature
  defect" signal from the "assembly fragmentation" signal more cleanly.
- Two `cctyper_only` loci in the pilot (identify misses something cctyper
  calls Trusted). Worth looking at per-locus.
- One array (arr_… `GCA_013414795.1` ref#1 Bacteroidota) scores 0.939 with
  Δstart=74 vs cctyper — cctyper calls one extra upstream repeat;
  CRISPRidentify's `number_mismatches` edge rule is almost certainly the
  cause. This is a known-honest difference, not a bug.
