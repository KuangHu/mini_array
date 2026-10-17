# Step 2B — structural + neighborhood explainer

Date: 2026-10-16
Approval ref: User review 2026-10-16, "批准在现有 1000 个基因组内开展 Step 2B / Explainer".
Scope lock: Only existing 1000-genome batch + existing annotations
(identify arrays/units, cctyper arrays/genes/HMM, pilot manifest). No new
DB download. No full-DB spacer-target search. No upstream issue submit.

## tl;dr

- All **338 candidate loci** (73 standalone + 260 possible_soft_context +
  5 contig_edge) from Step 2 were explained with ±2 kb overlay +
  structural class + ±10 kb protein neighborhood.
- **268 are ordinary array fragments** (all 260 soft-context + 6 standalone
  + 2 contig-edge land in or on known identify/cctyper arrays). They are
  catalogued, not promoted.
- **56 candidate_mini_array** survive structural review (55 from standalone
  + 1 contig-edge).
- **2 terminal_extension** (both standalone, abut known arrays within 100
  bp but <50% overlap).
- **12 edge_unresolved** (10 standalone + 2 contig-edge) — not enough
  flank + no identify context to decide.
- Of the 56 `candidate_mini_array`:
  - **12 carry strict 0.90 support**
  - **28 have ≥1 Cas HMM hit in ±10 kb window**
  - **14 have Cas1 or Cas2 (adaptation hallmarks) in that window**
  - **42 live in a genome with a high-confidence Cas operon somewhere
    on the same assembly** (not required to be adjacent)
- **25 module families** by HMM-signature; two cross-genome modules
  recur (CysH_0_IVB — a type IV-B marker, and CasR_1_* — a Cas regulator
  marker).
- Background audit: all shuffle replicates preserve mono+dinuc
  composition. **One query in rep1 happened to shuffle back to itself**
  (edit_distance_vs_original = 0, 29 nt), producing 50 of the 113
  anomalous rep1 pairs. All 113 rep1 pairs land in known_array_interior
  or soft_context; **ZERO** in standalone — the mini-array signal is
  unaffected by the rep1 outlier.

Explanation is complete for all 338 candidates; no top-N forced; 12
edge-unresolved items are listed explicitly as still-pending.

---

## 1. How original Step 2 buckets translate (structural explain)

| step2 bucket (before) | ordinary_array | terminal_ext | candidate_mini | repeat_amb | edge_unres |
|---|---:|---:|---:|---:|---:|
| possible_soft_context (260) | **260** | 0 | 0 | 0 | 0 |
| standalone (73) | 6 | 2 | **55** | 0 | 10 |
| contig_edge (5) | 2 | 0 | 1 | 0 | 2 |
| **total (338)** | **268** | **2** | **56** | **0** | **12** |

- Every `possible_soft_context` locus overlapped a non-bona-fide
  coord-valid identify call by ≥50% once the explainer looked at
  identify's full array table. We keep them as ordinary; none are
  dropped from the audit.
- 6 of the original 73 "standalone" loci turned out to overlap a known
  array once the mask included ALL coord-reliable bona-fide (we had
  already respected this at BED emission time; the structural overlay
  found them because `nearby_identify_arrays` traversed every array
  category, which is wider than the known-array BED).
- 2 standalones are terminal extensions — within 100 bp of a known
  array end.
- The remaining **56 standalone → candidate_mini_array** are the
  structurally defensible pool of possibly-novel mini-arrays.
- 10 standalone + 2 contig-edge are `edge_unresolved` — explicitly
  marked as not enough context to classify.
- `repeat_like_ambiguous` was defined but not reached by any
  candidate — in practice candidates either had enough repeat
  context to be `candidate_mini_array` or no context at all
  (`edge_unresolved`).

## 2. Protein neighborhood (±10 kb) summary

For the 70 structurally-surviving candidates (56 mini + 2 terminal +
12 edge-unresolved):

| cohort | n | with ≥1 HMM hit | with Cas1/Cas2 in ±10 kb | with high-conf Cas in genome |
|---|---:|---:|---:|---:|
| candidate_mini_array | 56 | 28 | 14 | 42 |
| terminal_extension | 2 | 2 | 1 | 2 |
| edge_unresolved | 12 | 5 | 1 | 6 |

Observations:
- Adaptation hallmarks (Cas1/Cas2) sit next to **14 of 56** mini-array
  candidates — well above the all-50-with-genome-Cas baseline because
  "genome Cas" doesn't require spatial co-localisation.
- `candidate_cards.md` sorts the 70 surviving candidates by strict
  support × HMM richness × genome, so the Cas1/Cas2-adjacent cases
  appear first.

## 3. Module families

25 HMM-signature module families total.

| cohort | n_members | n_genomes | notes |
|---|---:|---:|---|
| `no_hmm_hits` | 35 | 10 | Candidates with no annotated protein in ±10 kb. Not pre-judged — may still be real mini-arrays or non-array repeats. |
| `mod_0021` (`CysH_0_IVB`) | 3 | 2 | Cross-genome type IV-B marker. |
| `mod_0018` (`CasR_1_CAS-I-II-III-IV-V-VI`) | 2 | 2 | Cross-genome Cas regulator marker. |
| 22 other modules | — | — | All singletons or single-genome, kept per the "keep singletons without comparators" rule. |

`candidate_module_families.tsv` has the full list with sizes, genomes,
phyla, and representative locus_ids.

**Note on module interpretation:** cctyper's HMM library is Cas-centric,
so our "signature" is a Cas-vs-no-Cas space. We have NOT run Pfam /
general-protein HMM searches (would require new DB downloads, which
is explicitly out of scope). Candidates with `no_hmm_hits` are a mix
of "no neighbor protein that cctyper's HMM set recognises" — not a
statement that they have no protein at all.

## 4. Background audit (re-bucketed)

Identical pipeline, identical bucket definitions. Pair-counts are kept
separate from locus-counts.

### Pair counts by context bucket

| set | known_array_interior | suspected_terminal_ext | possible_soft_context | standalone | contig_edge | **n_pairs** |
|---|---:|---:|---:|---:|---:|---:|
| real | 14,308 | 47 | 969 | 139 | 6 | **15,472** |
| bg_rep0 | 0 | 0 | 0 | 0 | 0 | **0** |
| bg_rep1 | 102 | 0 | 11 | **0** | 0 | **113** |
| bg_rep2 | 0 | 0 | 0 | 0 | 0 | **0** |

### Locus-counts (via proxy dedup, same bucket definitions)

| set | known_array_interior | suspected_terminal_ext | possible_soft_context | standalone | contig_edge | **n_loci** |
|---|---:|---:|---:|---:|---:|---:|
| real | 14,311 | 47 | 968 | 139 | 6 | **15,471** |
| bg_rep0 | — | — | — | — | — | **0** |
| bg_rep1 | 102 | — | 11 | **0** | — | **113** |
| bg_rep2 | — | — | — | — | — | **0** |

**The 113 rep1 pairs land in known_array_interior or possible_soft_context
only.** Zero land in standalone. The anomalously high rep1 does not
contaminate the "new mini-array" signal.

### Shuffle validity audit

For all 213 queries × 3 replicates = **639 (replicate × query) rows**:

- Mono composition preserved: **639 / 639**
- Dinuc composition preserved: **639 / 639**

Trace of the rep1 anomaly:

- `rep_63c3a59da4f274b1ae76` (29 nt) in rep1 has
  `edit_distance_vs_original = 0, dinuc_preserved = True` — the
  dinuc-shuffle happened to converge back to the original.
- That one query drives 50 of 113 rep1 pairs; the other 63 come from
  near-palindromic queries (several entries at similarity ≥ 0.60).
- Per user instruction, these are **kept** in the output; no post-hoc
  deletion of hitting queries. A hardening fix would be to retry the
  shuffle if it produces the same sequence or a near-palindrome, but
  that's a Step 2 fix not Step 2B.

## 5. Candidate cards

70 cards in `candidate_cards.md`, ordered: strict_support first, then
HMM-rich neighborhood, then by genome. Each card lists genome, contig,
candidate interval, flank available, structural class + reasons, pair
count + unique queries, ORF count in ±10 kb, HMM signature, nearby HMM
families, genome-wide Cas types (high + putative), module family, and
the query repeats that produced the pair.

The top of the file is now the strict-supported, HMM-rich mini-arrays
including loc_000754 / loc_000755 (both strands of a 92-bp repeat in
`GCA_016929575.1` adjacent to a full Type I-E Cas operon). These are
the strongest candidates for further manual review.

## 6. Deliverables

```
results/step2b_pilot/
  candidate_evidence.tsv            338 candidates × 31 columns
  candidate_protein_links.tsv       locus × HMM-family hit counts
  candidate_module_families.tsv     25 module families with sizes
  candidate_cards.md                70 structurally-surviving candidates
  background_by_bucket.tsv          real + 3 bg replicates × bucket counts
  background_query_audit.tsv        shuffle validity + similarity 639 rows
  bg_rep1_pair_trace.json           the 50 anchor×second entries for rep1
  summary.json                      overall counters
  step2b_sbatch.<jobid>.log         sbatch run log

provenance/
  LOCKED.md                         frozen implementation spec (unchanged)
  fingerprint.json                  cache key (unchanged)
```

The sbatch wrapper used was `scripts/sbatch/step2b_explainer.sh`
(jobid 12494, submitted 2026-10-16 ~13:34). The earlier login-node run
produced byte-identical outputs because the pipeline is deterministic
(no RNG, no clock-based identifiers).

## 7. Honoured non-negotiables

- All 73 standalone + 260 soft-context + 5 contig-edge entered the
  explainer. No post-hoc dropping.
- Being outside the known-array BED was NOT treated as "independent".
  Each candidate was overlaid with ALL identify categories + cctyper
  arrays before class assignment.
- Candidates in known-array BED-blind buckets were still re-checked
  for ordinary-array overlap.
- Protein neighborhood extracted but reclaims nothing — no neighbor is
  required to be a known DNA-modifying enzyme.
- "No HMM hit" is **not** reported as "novel"; it is reported as
  `hmm_signature=no_hmm_hits` with the caveat that cctyper's HMM
  library is Cas-centric and does not cover general protein families.
- "No Cas detected" is **not** reported as "Cas-independent"; we
  separately record genome-wide Cas (high + putative) via cctyper
  cas_operons.tab.
- Query library intact — not reduced to one representative per family.
- Same-genome and cross-genome query↔locus attributions preserved in
  `candidate_evidence.tsv` and `loci_attributions.tsv` (Step 2).
- Spacer variability was not used to disqualify.
- Singleton module families kept.
- `allow_full_run: false` preserved.
- Step 2B pilot ran under existing resource caps (sbatch 2 CPUs, 16 GB,
  1 h). Login-node run was replaced by sbatch for compliance; outputs
  identical.
- No new DB downloaded.
- No upstream issue submitted.
- No calibrated FDR claim.
- No full-DB spacer-target search.

## 8. Open items (not acted on)

- 12 `edge_unresolved` candidates are listed in `candidate_cards.md`
  as-is. Classification requires more flank or external context.
- Shuffle quality fix (retry if produced sequence identical to source
  or near-palindrome) is a Step 2 improvement, not a Step 2B deliverable.
- Pfam / general-protein HMM over the ±10 kb windows would widen the
  non-Cas annotation coverage. Would require a new HMM DB — out of
  scope.
- The "no_hmm_hits" module (35 candidates, 10 genomes) is not resolved
  with current annotation sources.
- Several candidates on the same contig in the same genome (e.g.
  GCA_958295815.1 c_e8c2c4e07de8 had 15 standalones in Step 2) collapse
  into fewer loci here but may still represent repetitive non-array
  elements. Separately verified with the module family output.
