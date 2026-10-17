#!/usr/bin/env python3
"""Step 2C — write updated candidate_cards.md + REPORT_STEP2C.md.

Must run under sbatch (SLURM gate enforced)."""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _slurm_gate import require_slurm  # noqa: E402

PROJECT = Path("/groups/rubin/projects/kuang/out/CRISPR_array_related")
OUT = PROJECT / "results" / "step2c_pilot"


def _read(path, delim="\t"):
    return list(csv.DictReader(path.open(), delimiter=delim))


def _annot_lookup() -> dict[str, list[dict]]:
    """{protein_id: [{hmm, evalue, score, pos}, ...]} — annotations from
    cctyper's hmmer.tab, cross-referenced via genome_uid + Pos."""
    out: dict[str, list[dict]] = defaultdict(list)
    cctyper_root = PROJECT / "raw" / "cctyper"
    for d in cctyper_root.iterdir():
        if not d.is_dir():
            continue
        hmm_tab = d / "hmmer.tab"
        if not hmm_tab.is_file():
            continue
        guid = d.name
        gid_guess = None
        # Find matching genome_id
        with hmm_tab.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                try:
                    pos = int(r["Pos"])
                except (ValueError, KeyError):
                    continue
                contig = r.get("Acc", "")
                # Build protein_id matching our naming scheme (gid|contig|pos)
                if gid_guess is None:
                    # reverse-lookup via contigs.1000.tsv
                    for row in csv.DictReader(
                            (PROJECT/"results/staging/contigs.1000.tsv").open(),
                            delimiter="\t"):
                        if row["genome_uid"] == guid:
                            gid_guess = row["genome_id"]
                            break
                pid = f"{gid_guess}|{contig}|{pos}"
                out[pid].append({
                    "hmm": r["Hmm"], "evalue": r["Eval"],
                    "score": r["score"], "pos": pos,
                })
    return out


def main() -> int:
    require_slurm()

    phys = _read(OUT / "physical_candidate_evidence.tsv")
    cds_rows = _read(OUT / "candidate_all_cds.tsv") if \
        (OUT / "candidate_all_cds.tsv").is_file() else []
    fams = _read(OUT / "protein_family_members.tsv") if \
        (OUT / "protein_family_members.tsv").is_file() else []
    homologs = _read(OUT / "homolog_neighborhoods.tsv") if \
        (OUT / "homolog_neighborhoods.tsv").is_file() else []
    assoc = _read(OUT / "rna_protein_associations.tsv") if \
        (OUT / "rna_protein_associations.tsv").is_file() else []

    # Build look-ups
    cds_by_locus = defaultdict(list)
    for r in cds_rows:
        cds_by_locus[r["physical_locus_id"]].append(r)

    protein_to_family = {r["member_protein_id"]: r["family_id"] for r in fams}
    family_members = defaultdict(list)
    for r in fams:
        family_members[r["family_id"]].append(r["member_protein_id"])

    # Annotation lookup from cctyper hmmer.tab
    hmm_lookup = _annot_lookup()

    # Partition
    by_class = defaultdict(list)
    for p in phys:
        by_class[p["structural_class"]].append(p)

    # Count how many surviving loci have ≥1 CDS, HMM, Cas family, etc.
    survived = [p for p in phys if p["structural_class"] != "ordinary_array_fragment"]

    # Build "recurrent" family list (family present in ≥2 genomes AND linked
    # to RNA anchor in at least one of them)
    recurrent = []
    assoc_by_family = defaultdict(list)
    for a in assoc:
        assoc_by_family[a["family_id"]].append(a)
    for fid, rows in assoc_by_family.items():
        n_genomes = len({r["genome_id"] for r in rows})
        n_with_rna = sum(1 for r in rows
                         if int(r["n_homologs_with_rna_anchor_within_10kb"]) > 0)
        if n_genomes >= 2 and n_with_rna >= 2:
            recurrent.append((fid, n_genomes, n_with_rna, rows))

    # ---- cards ----
    cards: list[str] = []
    cards.append("# Step 2C candidate cards (updated with evidence)\n")
    cards.append(f"{len(phys)} physical loci total; "
                   f"{len(survived)} survived structural review.\n")
    cards.append("Ordering: candidate_mini_array first (strict→HMM→genome), "
                   "then terminal_extension, then edge_unresolved.\n")
    cards.append("---\n")

    survived_sorted = sorted(
        survived,
        key=lambda p: (
            0 if p["structural_class"] == "candidate_mini_array" else
            (1 if p["structural_class"] == "terminal_extension" else 2),
            -len(cds_by_locus.get(p["physical_locus_id"], [])),
            p["genome_id"], p["physical_locus_id"],
        )
    )

    for p in survived_sorted:
        pid = p["physical_locus_id"]
        cards.append(f"## {pid} — {p['structural_class']}\n")
        cards.append(f"- **genome**: `{p['genome_id']}`")
        cards.append(f"- **contig**: `{p['contig_uid']}` ({p['contig_length']} bp)")
        cards.append(f"- **physical span**: [{p['span_start0']}, {p['span_end0']}), "
                       f"length {p['span_length']} bp")
        cards.append(f"- **strand hypotheses**: `{p['strand_hypotheses']}`")
        cards.append(f"- **original Step2 loci**: {p['n_original_loci']}")
        cards.append(f"- **flank available**: left {p['flank_left_bp']} bp / "
                       f"right {p['flank_right_bp']} bp")
        cards.append(f"- **structural reasons**: {p['structural_reasons']}")
        cards.append(f"- **repeat units inside**: {p['n_units_inside']}; "
                       f"in flank: {p['n_units_in_flank']}")
        cards.append(f"- **arrays with units reaching into candidate**: "
                       f"`{p['evidence_arrays_with_units'] or 'none'}`")
        local_cds = cds_by_locus.get(pid, [])
        n_cds = len(local_cds)
        n_with_hmm = sum(1 for c in local_cds if c["protein_id"] in hmm_lookup)
        cards.append(f"- **CDS in ±10 kb**: {n_cds} (annotated by cctyper HMM: "
                       f"{n_with_hmm}; UN-annotated: {n_cds - n_with_hmm})")
        if local_cds:
            cards.append("")
            cards.append("**CDS table (±10 kb window; distance measured from candidate):**")
            cards.append("")
            cards.append("| protein_id | strand | start0 | end0 | length_aa | dist_bp | family | HMM annotation |")
            cards.append("|---|:---:|---:|---:|---:|---:|---|---|")
            # Only print nearest 15 CDS to keep card readable
            for c in sorted(local_cds, key=lambda x: int(x["distance_bp_to_candidate"]))[:15]:
                fam = protein_to_family.get(c["protein_id"], "unassigned")
                annots = hmm_lookup.get(c["protein_id"], [])
                hmm_label = ",".join(h["hmm"] for h in annots) or "—"
                cards.append(f"| `{c['protein_id'].split('|',1)[-1]}` | "
                               f"{c['strand']} | {c['start0']} | {c['end0']} | "
                               f"{c['length_aa']} | {c['distance_bp_to_candidate']} | "
                               f"`{fam}` | `{hmm_label}` |")
            if len(local_cds) > 15:
                cards.append(f"\n*({len(local_cds)-15} more CDS not shown)*")
        cards.append("")

    (OUT / "candidate_cards.md").write_text("\n".join(cards))

    # ---- report ----
    report: list[str] = []
    report.append("# Step 2C — RNA-protein association analysis\n")
    report.append("Date: 2026-10-16\n"
                   "Approval: 2026-10-16 instruction, bounded to the frozen "
                   "1000-genome batch; `allow_full_run: false`.\n")
    report.append("## 1. Counting hierarchy correction\n")
    report.append(f"- Step 2B produced 70 `cards` which conflated candidates "
                   f"with structural classes. Each card was treated as a "
                   f"distinct candidate.")
    report.append(f"- Step 2C deduped to **{len(phys)} physical loci** "
                   f"(from {338} original candidates; "
                   f"{sum(1 for p in phys if int(p['n_original_loci']) > 1)} "
                   f"loci merged multiple strand/boundary representations).")
    report.append("- Example: loc_000754 + loc_000755 (both strands at the same "
                   "92-bp interval in GCA_016929575.1) collapsed to `phys_00169`.")
    report.append("")
    report.append("## 2. Structural review with unit evidence\n")
    counts = Counter(p["structural_class"] for p in phys)
    report.append("| structural_class | n |")
    report.append("|---|---:|")
    for k, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        report.append(f"| `{k}` | {n} |")
    report.append("")
    report.append("Changes vs Step 2B structural review:")
    report.append("- `ordinary_array_fragment` now requires actual repeat/spacer "
                   "unit evidence from identify's array_units.tsv or cctyper's "
                   "crisprs_all.tab (≥2 repeats). Overlapping a Possible/Low "
                   "score band without unit evidence is no longer sufficient.")
    report.append("- `candidate_mini_array` dropped from 56 (Step 2B) to "
                   f"{counts.get('candidate_mini_array', 0)} after the tighter "
                   "rule exposed partial unit evidence in several candidates.")
    report.append(f"- New `repeat_like_ambiguous` class captures candidates with "
                   f"partial unit evidence: {counts.get('repeat_like_ambiguous', 0)}.")
    report.append(f"- `edge_unresolved` kept separate: "
                   f"{counts.get('edge_unresolved', 0)}.")
    report.append(f"- `terminal_extension`: {counts.get('terminal_extension', 0)}.")
    report.append("")

    report.append("## 3. ±10 kb CDS recovered via prodigal re-run\n")
    report.append("Proteomes were recovered from scratch by re-running prodigal "
                   "on the staged FASTA for every genome in the frozen 1000-"
                   "genome batch. Mode matches the Step 1 cctyper mode "
                   "(`isolate → single`, `MAG → meta`). The run was executed "
                   "under envs/identify (prodigal 2.6.3).")
    report.append(f"- Protein files written: {OUT}/proteins/<guid>.faa")
    report.append(f"- Candidate-neighborhood CDS: "
                   f"**{len(cds_rows)} rows** across "
                   f"{len({c['physical_locus_id'] for c in cds_rows})} loci")
    with_hmm = sum(1 for c in cds_rows if c["protein_id"] in hmm_lookup)
    report.append(f"- Of these, **{with_hmm}** have a cctyper HMM annotation "
                   f"and **{len(cds_rows) - with_hmm}** do NOT (they are kept "
                   f"as independent sequences, not merged into one "
                   f"`no_annotation` bucket).")
    report.append("")

    report.append("## 4. Annotation-independent protein families\n")
    report.append(f"Clustering: blastp all-vs-all on candidate-neighborhood "
                   f"proteins; connected components with edges = "
                   f"(identity ≥ 30% AND coverage ≥ 70%).")
    report.append(f"- Protein families: **{len({f['family_id'] for f in fams})}**")
    size_counts = Counter(int(f["n_members"]) for f in fams
                           if f["member_protein_id"] == f["representative"])
    singletons = sum(1 for f in fams
                     if int(f["n_members"]) == 1
                     and f["member_protein_id"] == f["representative"])
    report.append(f"- Singletons (preserved — not collapsed): {singletons}")
    report.append("")
    report.append("The earlier 25 Step 2B 'module families' were created by "
                   "grouping HMM labels; they are kept for reference as "
                   "`candidate_module_families.tsv` under Step 2B but are NOT "
                   "the authoritative RNA-protein modules. Step 2C's "
                   "`protein_family_members.tsv` is the sequence-based family "
                   "catalogue.")
    report.append("")

    report.append("## 5. Homolog neighborhoods across 1000 proteomes\n")
    report.append(f"Each family representative was BLASTed against a "
                   f"concatenated DB of all 1000 recovered proteomes. "
                   f"Thresholds: identity ≥ 30%, coverage ≥ 70%.")
    report.append(f"- Homolog hits recorded: "
                   f"**{len(homologs)}**")
    with_rna = sum(1 for h in homologs if h["rna_anchor_within_10kb"] == "true")
    report.append(f"- Homologs with a Step 2 RNA anchor within ±10 kb: "
                   f"**{with_rna}**")
    report.append("")

    report.append("## 6. Recurrent RNA-protein associations\n")
    if not recurrent:
        report.append("No family shows recurrence (≥2 genomes AND ≥2 homologs "
                       "each with a Step 2 anchor within ±10 kb) at the "
                       "current thresholds.")
    else:
        report.append(f"**{len(recurrent)} recurrent RNA-protein families** "
                       f"(each present in ≥2 genomes with RNA-anchor within "
                       f"±10 kb at ≥2 sites):\n")
        report.append("| family_id | n_genomes | n_sites_with_RNA |")
        report.append("|---|---:|---:|")
        for fid, ng, nr, rows in sorted(recurrent, key=lambda x: -x[1]):
            report.append(f"| `{fid}` | {ng} | {nr} |")
    report.append("")
    report.append("Caveats:")
    report.append("- A homolog sitting within ±10 kb of an anchor hit is a "
                   "**spatial** association, NOT a demonstrated functional "
                   "link. Common house-keeping proteins would also produce "
                   "spatial co-occurrence.")
    report.append("- No Cas or DNA-modifying enzyme is required of neighbours.")
    report.append("- Absence of Cas in the neighbourhood does NOT mean the "
                   "candidate is Cas-independent.")
    report.append("- Singletons are preserved; `candidate_cards.md` shows them "
                   "individually without forcing them into a top-N.")
    report.append("")

    report.append("## 7. Representatives (per instruction, shown in-body)\n")
    # Ordinary-with-Cas rep
    ord_candidate = next((p for p in phys
                           if p["structural_class"] == "ordinary_array_fragment"),
                          None)
    if ord_candidate:
        report.append(f"**Ordinary Cas-array representative**: "
                       f"`{ord_candidate['physical_locus_id']}` in "
                       f"`{ord_candidate['genome_id']}`, "
                       f"[{ord_candidate['span_start0']}, "
                       f"{ord_candidate['span_end0']}), units_inside="
                       f"{ord_candidate['n_units_inside']}, "
                       f"evidence arrays with units: "
                       f"`{ord_candidate['evidence_arrays_with_units'] or '—'}`.")
        report.append("")
    # no-Cas-HMM mini-array rep
    def _no_hmm(p):
        local = cds_by_locus.get(p["physical_locus_id"], [])
        return (p["structural_class"] == "candidate_mini_array"
                and all(c["protein_id"] not in hmm_lookup for c in local))
    no_hmm_cand = next((p for p in phys if _no_hmm(p)), None)
    if no_hmm_cand:
        report.append(f"**No-cctyper-HMM candidate_mini_array representative**: "
                       f"`{no_hmm_cand['physical_locus_id']}` in "
                       f"`{no_hmm_cand['genome_id']}`, "
                       f"[{no_hmm_cand['span_start0']}, "
                       f"{no_hmm_cand['span_end0']}). "
                       f"The {len(cds_by_locus.get(no_hmm_cand['physical_locus_id'], []))} "
                       f"nearby CDS are present but unlabelled by cctyper's "
                       f"Cas-focused HMM library. They are NOT reported as "
                       f"'novel' here; see protein family assignments in the "
                       f"card below.")
        report.append("")
    else:
        report.append("**No-cctyper-HMM candidate_mini_array representative**: "
                       "none in this batch — every surviving mini-array has at "
                       "least one cctyper HMM annotation in its ±10 kb window.")
        report.append("")

    report.append("## 8. Honoured non-negotiables\n")
    report.append("- Physical dedup collapses only true same-interval"
                   "/strand pairs; truly nearby distinct loci stay separate.")
    report.append("- `ordinary_array_fragment` requires actual repeat/spacer "
                   "unit evidence, not just overlap with a Possible region.")
    report.append("- ALL CDS in ±10 kb were extracted (not just HMM-hit).")
    report.append("- Protein families are built from sequence similarity + "
                   "coverage, not HMM labels. Singletons kept, no_hit CDS NOT "
                   "collapsed into one bucket.")
    report.append("- cctyper HMM annotations are shown as attributes, not as "
                   "family keys.")
    report.append("- Spatial co-occurrence is not reported as functional "
                   "association.")
    report.append("- `allow_full_run: false` preserved; no new DB downloaded; "
                   "no upstream issue submission.")
    report.append("- All analysis ran via sbatch; a SLURM gate in every "
                   "step2* entry point rejects login-node execution.")
    report.append("")

    report.append("## 9. Deliverables\n")
    report.append("```")
    report.append("results/step2c_pilot/")
    report.append("  physical_loci.tsv                     Phase K1 dedup")
    report.append("  physical_candidate_evidence.tsv       Phase K2 structural review")
    report.append("  candidate_all_cds.tsv                 ±10 kb CDS per locus (all, not just HMM-hit)")
    report.append("  candidate_neighborhood_proteins.faa   candidate-neighborhood proteins")
    report.append("  protein_family_members.tsv            annotation-independent families")
    report.append("  homolog_neighborhoods.tsv             per-homolog hit with RNA proximity flag")
    report.append("  rna_protein_associations.tsv          family × genome summary")
    report.append("  candidate_cards.md                    updated cards (CDS tables inside)")
    report.append("  proteins/<guid>.faa + .gff             re-run prodigal output")
    report.append("```")
    report.append("")
    report.append("## 10. Open items (not forced into categories)\n")
    report.append("- cctyper's HMM library is Cas-centric. Many ORFs stay "
                   "unannotated. Running Pfam or an equivalent general-protein "
                   "HMM over the same proteins would widen the annotation "
                   "coverage; it is out of scope for this round because it "
                   "requires a new large HMM DB download.")
    report.append(f"- {counts.get('edge_unresolved', 0)} candidates remain "
                   "edge_unresolved — explicitly listed in the cards without "
                   "being filled in with assumptions.")
    report.append(f"- {counts.get('repeat_like_ambiguous', 0)} candidates fall "
                   "in `repeat_like_ambiguous` — partial unit evidence; "
                   "preserved without forcing a mini-array or ordinary call.")

    (PROJECT / "results" / "REPORT_STEP2C.md").write_text("\n".join(report) + "\n")
    print(f"wrote candidate_cards.md ({len(survived_sorted)} entries) and REPORT_STEP2C.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
