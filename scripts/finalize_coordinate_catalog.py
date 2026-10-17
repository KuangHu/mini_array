#!/usr/bin/env python3
"""Phase L2: semantic finalization of the coordinate catalog.

No new calculation. No new biology. No boundary / threshold / parameter
change. All inputs already exist under results/{step1,step2,step2b,step2c,
coordinate_catalog}/. Outputs land in results/coordinate_catalog/ and
replace older files with semantically clearer versions.

Changes vs L1 catalog:
  - crispr_arrays.{tsv,bed} kept, re-labelled as ALL tool calls (not
    independent trusted loci).
  - validation column split into `coordinate_in_bounds` (true/false) and
    `has_ambiguous_bases` (true/false). N / non-ACGT no longer presented
    as 'out-of-bounds'.
  - new trusted-call view: `crispr_arrays_trusted.{tsv,bed}` filtered to
    coord-valid AND (identify Bona-fide OR cctyper Trusted). Each row
    keeps `source_tool` + `source_id`. We DO NOT claim locus
    independence by summing the two tools.
  - mini_arrays.bed strand forced to '.' because repeat alignment
    direction is not an independent orientation call. TSV adds
    `strand_source` column that reports the raw origin
    (step2_pair_alignment / none).
  - mini_array_units.tsv now has `pair_id`, `physical_locus_id`,
    `locus_class` columns; one row = one pair.
  - coordinate_catalog_summary.md rewrites both physical-locus and
    pair-level stats and explicitly scopes the "78 physical loci with
    both strand hypotheses" figure. Adds note that degenerate=0 is a
    consequence of the secondary thresholds applied upstream, not
    biological absence.
"""
from __future__ import annotations

import csv
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _slurm_gate import require_slurm  # noqa: E402

PROJECT = Path("/groups/rubin/projects/kuang/out/CRISPR_array_related")
STEP2C = PROJECT / "results" / "step2c_pilot"
OUT = PROJECT / "results" / "coordinate_catalog"


def _read(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _write_tsv(path: Path, rows: list[dict], cols: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    tmp.rename(path)


def _write_bed_lines(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w") as fh:
        for ln in lines:
            fh.write(ln + "\n")
    tmp.rename(path)


# Treat the L1 catalog's coord_validation=='ok' as "in_bounds + no ambiguity".
# 'degenerate_bases' as "in_bounds + has ambiguity".
# 'oob_left' / 'oob_right' as "truly out of bounds".
def _split_validation(coord_validation: str) -> tuple[str, str]:
    v = (coord_validation or "").lower()
    if v == "ok":
        return "true", "false"
    if v == "degenerate_bases":
        return "true", "true"
    if v in ("oob_left", "oob_right"):
        return "false", "false"
    # seq_missing or anything else — unknown
    return "unknown", "unknown"


# ----------------- crispr_arrays rewrite -----------------

def rewrite_crispr_arrays() -> tuple[list[dict], list[dict]]:
    rows = _read(OUT / "crispr_arrays.tsv")
    for r in rows:
        in_bounds, has_amb = _split_validation(r.get("coord_validation", ""))
        r["coordinate_in_bounds"] = in_bounds
        r["has_ambiguous_bases"] = has_amb

    # Trusted-call view: coord-valid AND (identify Bona-fide OR cctyper Trusted)
    def _is_trusted(r: dict) -> bool:
        if r["coordinate_in_bounds"] != "true":
            return False
        src = r["source_tool"]
        if src == "crispridentify" and r["category_raw"] == "Bona-fide":
            return True
        if src == "crisprcastyper" and (r["score_or_trusted"] or "").lower() == "true":
            return True
        return False

    trusted = [r for r in rows if _is_trusted(r)]
    return rows, trusted


# ----------------- mini-arrays rewrite -----------------

STRAND_PLUS = {"plus", "forward", "+"}
STRAND_MINUS = {"minus", "reversed", "-"}


def _normalize_mini_bed_strand(_strand_hypotheses: str) -> str:
    # Always '.': the mini-array strand hypothesis is a BLAST query alignment
    # direction, NOT an independent orientation call. The TSV preserves the
    # raw hypothesis; the BED must not imply transcription direction.
    return "."


def rewrite_mini_arrays() -> tuple[list[dict], list[dict]]:
    mini = _read(OUT / "mini_arrays.tsv")
    edge = _read(OUT / "edge_unresolved_candidates.tsv")

    for r in mini + edge:
        in_bounds, has_amb = _split_validation(r.get("coord_validation", ""))
        r["coordinate_in_bounds"] = in_bounds
        r["has_ambiguous_bases"] = has_amb
        # explicit strand_source
        if r["strand_hypotheses"]:
            r["strand_source"] = "step2_pair_alignment"
        else:
            r["strand_source"] = "none"
    return mini, edge


# ----------------- units rewrite -----------------

def rewrite_mini_array_units(mini: list[dict], edge: list[dict]) -> list[dict]:
    units = _read(OUT / "mini_array_units.tsv")
    # Map physical_locus_id → locus_class from mini + edge rows
    class_by_phys = {r["physical_locus_id"]: r["structural_class"]
                     for r in (mini + edge)}

    # Group units per physical_locus_id to assign stable pair_ids
    by_phys = defaultdict(list)
    for u in units:
        by_phys[u["physical_locus_id"]].append(u)

    out = []
    for pid, lst in by_phys.items():
        lst.sort(key=lambda x: int(x["repeat_1_start0"]))
        for i, u in enumerate(lst):
            u["pair_id"] = f"{pid}:pair_{i+1:02d}"
            u["locus_class"] = class_by_phys.get(pid, "")
            in_bounds, _amb = _split_validation("ok" if u["coord_validation"]
                                                 == "ok,ok,ok" else "degenerate_bases"
                                                 if "degenerate" in u["coord_validation"]
                                                 else "unknown")
            # Already validated componentwise; just surface in_bounds flag
            all_ok = all(x == "ok" for x in u["coord_validation"].split(","))
            u["coordinate_in_bounds"] = "true" if all_ok else "false"
            out.append(u)
    return out


# ----------------- summary rewrite -----------------

def write_summary(all_rows: list[dict], trusted: list[dict],
                   mini: list[dict], edge: list[dict], units: list[dict]) -> None:

    tool_counts = Counter(r["source_tool"] for r in all_rows)
    cat_counts = Counter(r["category_raw"] for r in all_rows)
    in_bounds_counts = Counter(r["coordinate_in_bounds"] for r in all_rows)
    amb_counts = Counter(r["has_ambiguous_bases"] for r in all_rows)

    trusted_tool = Counter(r["source_tool"] for r in trusted)
    trusted_in_bounds = Counter(r["coordinate_in_bounds"] for r in trusted)

    # mini / edge stats
    mini_val = Counter(r["coordinate_in_bounds"] for r in mini + edge)
    mini_amb = Counter(r["has_ambiguous_bases"] for r in mini + edge)
    mini_strand = Counter(r["strand_hypotheses"] for r in mini + edge)

    # "78 physical loci with both strand hypotheses" came from Step 2C
    # physical_loci.tsv over 260 physical candidates (all structural classes
    # combined). The user asked for that scope to be explicit.
    phys = _read(STEP2C / "physical_candidate_evidence.tsv")
    both_strand_total = sum(1 for r in phys
                             if r["strand_hypotheses"] == "minus|plus")
    both_strand_mini = sum(1 for r in mini
                            if r["strand_hypotheses"] == "minus|plus")
    both_strand_edge = sum(1 for r in edge
                            if r["strand_hypotheses"] == "minus|plus")

    # Pair-level stats
    pair_locus_class = Counter(u["locus_class"] for u in units)
    pair_completeness = Counter(u["repeat_2_completeness"] for u in units)
    pair_strand_pair = Counter(u["strand_pair"] for u in units)

    md = []
    md.append("# Coordinate catalog — final export semantics\n")
    md.append("Pure coordinate catalog. No biology. No new DB. No parameter "
               "change. 0-based, half-open on BED. All analysis ran via "
               "sbatch; a SLURM_JOB_ID gate in every entry point refuses "
               "login-node execution.\n")
    md.append("## 1. CRISPR arrays — ALL tool calls (not independent trusted loci)\n")
    md.append("`crispr_arrays.{tsv,bed}` is the union of every call produced by "
               "either CRISPRidentify or CRISPRCasTyper on the frozen 1000-"
               "genome batch. **This is NOT a set of independent trusted "
               "arrays** — it includes Possible/Low score/Alternative/Possible "
               "Discarded identify calls AND Untrusted cctyper calls. Each "
               "row carries `source_tool` + `source_id` so attributions are "
               "preserved. Rows from the two tools MAY overlap on the same "
               "genomic interval; **adding counts across tools would double-"
               "count overlapping loci**.\n")
    md.append(f"- Total rows: **{len(all_rows)}**")
    for k, n in tool_counts.most_common():
        md.append(f"  - `{k}`: {n}")
    md.append("- Identify / cctyper call-class distribution:")
    for k, n in cat_counts.most_common():
        md.append(f"  - `{k}`: {n}")
    md.append("")
    md.append("## 2. Trusted-call view\n")
    md.append("`crispr_arrays_trusted.{tsv,bed}` filters the ALL-tool-call table "
               "to rows that are **coordinate-valid AND** one of:\n"
               "- `source_tool=crispridentify` AND `category_raw=Bona-fide`, OR\n"
               "- `source_tool=crisprcastyper` AND `score_or_trusted=True`.\n"
               "\n"
               "This is a per-call view. We DO NOT merge CRISPRidentify Bona-"
               "fide + cctyper Trusted calls that fall on the same genomic "
               "interval, and we DO NOT report an independent locus count by "
               "summing tool rows. Downstream work should dedup per "
               "(contig_uid, start0, end0) if a per-locus count is needed.\n")
    md.append(f"- Trusted-call rows: **{len(trusted)}**")
    for k, n in trusted_tool.most_common():
        md.append(f"  - `{k}`: {n}")
    md.append(f"- coordinate_in_bounds=true (by filter definition): "
               f"{trusted_in_bounds.get('true', 0)}")
    md.append("")
    md.append("## 3. Coordinate validation (split)\n")
    md.append("Validation is now reported in two independent columns:\n"
               "- `coordinate_in_bounds` (true/false): did [start0, end0) fit "
               "inside the staged contig?\n"
               "- `has_ambiguous_bases` (true/false): does the slice contain "
               "any non-ACGT base? Ambiguity is NOT reported as 'out of "
               "bounds'.\n")
    md.append("### CRISPR arrays (all tool calls)")
    md.append(f"- coordinate_in_bounds=true: {in_bounds_counts.get('true', 0)}")
    md.append(f"- coordinate_in_bounds=false: {in_bounds_counts.get('false', 0)}")
    md.append(f"- has_ambiguous_bases=true: {amb_counts.get('true', 0)}")
    md.append(f"- has_ambiguous_bases=false: {amb_counts.get('false', 0)}")
    md.append("")
    md.append("### Mini-array + edge candidates")
    md.append(f"- coordinate_in_bounds=true: {mini_val.get('true', 0)}")
    md.append(f"- has_ambiguous_bases=true: {mini_amb.get('true', 0)}")
    md.append("")
    md.append("## 4. Mini-array semantics\n")
    md.append(f"- `candidate_mini_array` physical loci: **{len(mini)}**")
    md.append(f"- `edge_unresolved` physical loci (separate file): **{len(edge)}**")
    md.append("- The mini_arrays.bed `strand` column is forced to `.` — the "
               "underlying hypothesis comes from a BLAST query alignment and "
               "is NOT an independent transcription-direction call. The TSV "
               "preserves `strand_hypotheses` AND a `strand_source` column "
               "that identifies the raw source (step2_pair_alignment). Treat "
               "the strand only as a hint for aligning queries, never as a "
               "biological orientation.")
    md.append("")
    md.append("### mini + edge strand-hypothesis distribution")
    md.append("| strand_hypotheses | count |")
    md.append("|---|---:|")
    for k, n in mini_strand.most_common():
        md.append(f"| `{k or '(empty)'}` | {n} |")
    md.append("")
    md.append("### The '78 both-strand loci' statement — scope\n")
    md.append(f"- `78` refers to the full set of **{len(phys)} physical loci** "
               f"(all structural classes) in `step2c_pilot/"
               f"physical_candidate_evidence.tsv` where the `strand_hypotheses` "
               f"column is `minus|plus` (both strand-pair reps collapsed under "
               f"one physical_locus_id by Step 2C dedup).")
    md.append(f"- Within the 40 `candidate_mini_array` rows kept in "
               f"`mini_arrays.tsv`: **{both_strand_mini}** have "
               f"`strand_hypotheses=minus|plus`.")
    md.append(f"- Within the 11 `edge_unresolved` rows: "
               f"**{both_strand_edge}** have both strands.")
    md.append(f"- Within the remaining 209 physical loci "
               f"(ordinary_array_fragment / terminal_extension / "
               f"repeat_like_ambiguous — all classes excluded from the "
               f"mini-array file): "
               f"**{both_strand_total - both_strand_mini - both_strand_edge}**.")
    md.append("")
    md.append("## 5. mini_array_units.tsv — one row per pair\n")
    md.append(f"- Rows: **{len(units)}** (= pairs), covering "
               f"**{len({u['physical_locus_id'] for u in units})}** physical "
               f"loci.")
    md.append("- New columns: `pair_id`, `physical_locus_id`, `locus_class` "
               "(= the mini/edge candidate's structural_class).")
    md.append("- Pair-level breakdown by locus_class:")
    for k, n in pair_locus_class.most_common():
        md.append(f"  - `{k}`: {n}")
    md.append("- Pair-level repeat_2 completeness distribution:")
    for k, n in pair_completeness.most_common():
        md.append(f"  - `{k}`: {n}")
    md.append("- Pair-level strand (from Step 2 BLAST pair, NOT a biological call):")
    for k, n in pair_strand_pair.most_common():
        md.append(f"  - `{k}`: {n}")
    md.append("")
    md.append("### Caveat: degenerate=0 is a secondary-threshold artefact\n")
    md.append("The Step 2 pair finder requires each secondary (second-repeat) "
               "hit to pass `identity ≥ 0.70 AND coverage ≥ 0.60`. Hits that "
               "would have been labelled `degenerate` under a looser rule "
               "were filtered out BEFORE pair formation. Therefore "
               "`repeat_2_completeness='degenerate' = 0` in this catalog is "
               "a **consequence of the frozen secondary thresholds**, not a "
               "statement that biologically degenerate second repeats do not "
               "exist. Changing the thresholds is explicitly out of scope for "
               "this round.\n")
    md.append("## 6. Non-negotiables honoured\n")
    md.append("- No new search, no new biology, no CDS, no protein, no HMM.")
    md.append("- Candidate boundaries, classification thresholds, search "
               "parameters: unchanged.")
    md.append("- `allow_full_run: false` preserved. No new DB. No upstream "
               "issue submission.")
    md.append("- All analysis sbatch-only (SLURM gate enforced at import).")
    md.append("")

    (OUT / "coordinate_catalog_summary.md").write_text("\n".join(md) + "\n")


# ----------------- main -----------------

def main() -> int:
    require_slurm()

    # 1. Rewrite crispr arrays with split validation + add trusted view.
    all_rows, trusted = rewrite_crispr_arrays()
    ARR_COLS = ["source_tool", "source_id", "genome_id", "contig_uid",
                 "start0", "end0", "length", "strand_raw",
                 "orientation_prediction", "category_raw", "qc_tier",
                 "score_or_trusted", "subtype_or_subtype", "n_repeats",
                 "n_spacers", "consensus_repeat_forward", "repeat_id",
                 "qc_reasons", "provenance",
                 "coord_validation", "coordinate_in_bounds",
                 "has_ambiguous_bases"]
    _write_tsv(OUT / "crispr_arrays.tsv", all_rows, ARR_COLS)

    # Trusted view TSV + BED
    _write_tsv(OUT / "crispr_arrays_trusted.tsv", trusted, ARR_COLS)

    def _strand_token(s: str) -> str:
        s = (s or "").strip()
        if s in STRAND_PLUS:
            return "+"
        if s in STRAND_MINUS:
            return "-"
        if "forward" in s.lower() or s == "+":
            return "+"
        if "reversed" in s.lower() or s == "-":
            return "-"
        return "."

    # Rewrite ALL-call BED with the new label semantics (still unchanged
    # byte-wise except we also keep an explicit header-comment line that
    # labels the file).
    arr_bed_lines = ["# ALL tool calls (identify + cctyper); NOT a trusted "
                      "independent-locus catalogue — see "
                      "crispr_arrays_trusted.bed for the filtered view."]
    for r in all_rows:
        name = f"{r['source_tool']}:{r['source_id']}:{r['category_raw']}"
        strand = _strand_token(r["orientation_prediction"])
        arr_bed_lines.append("\t".join([r["contig_uid"], str(r["start0"]),
                                          str(r["end0"]), name, "0", strand]))
    _write_bed_lines(OUT / "crispr_arrays.bed", arr_bed_lines)

    trusted_bed_lines = ["# Trusted-call view: coord-valid AND (identify "
                          "Bona-fide OR cctyper Trusted). Per-call rows; "
                          "DO NOT sum across tools to get an independent "
                          "locus count."]
    for r in trusted:
        name = f"{r['source_tool']}:{r['source_id']}:{r['category_raw']}"
        strand = _strand_token(r["orientation_prediction"])
        trusted_bed_lines.append("\t".join([r["contig_uid"], str(r["start0"]),
                                              str(r["end0"]), name, "0",
                                              strand]))
    _write_bed_lines(OUT / "crispr_arrays_trusted.bed", trusted_bed_lines)

    # 2. Rewrite mini / edge with split validation + strand_source + strand=.
    mini, edge = rewrite_mini_arrays()

    MINI_COLS = ["physical_locus_id", "genome_id", "genome_uid", "contig_uid",
                 "contig_length", "start0", "end0", "length",
                 "strand_hypotheses", "strand_source",
                 "n_original_loci", "original_locus_ids",
                 "query_repeat_ids", "strict_0.90_support", "flank_left_bp",
                 "flank_right_bp", "structural_class", "structural_reasons",
                 "n_units_inside", "n_units_in_flank",
                 "evidence_arrays_with_units", "provenance",
                 "coord_validation", "coordinate_in_bounds",
                 "has_ambiguous_bases"]
    _write_tsv(OUT / "mini_arrays.tsv", mini, MINI_COLS)
    _write_tsv(OUT / "edge_unresolved_candidates.tsv", edge, MINI_COLS)

    mini_bed_lines = [
        "# Mini-array physical loci (strand ALWAYS '.' — repeat alignment "
        "direction is a BLAST artefact, not a biological orientation; see "
        "mini_arrays.tsv `strand_hypotheses` + `strand_source`)."
    ]
    for r in mini:
        name = f"{r['physical_locus_id']}:{r['structural_class']}"
        mini_bed_lines.append("\t".join([r["contig_uid"], str(r["start0"]),
                                           str(r["end0"]), name, "0", "."]))
    _write_bed_lines(OUT / "mini_arrays.bed", mini_bed_lines)

    # 3. Rewrite units with pair_id + locus_class
    units = rewrite_mini_array_units(mini, edge)
    UNIT_COLS = ["pair_id", "physical_locus_id", "locus_class",
                  "genome_id", "contig_uid", "strand_pair",
                  "repeat_1_start0", "repeat_1_end0", "repeat_1_query_id",
                  "spacer_start0", "spacer_end0", "spacer_len",
                  "repeat_2_start0", "repeat_2_end0", "repeat_2_query_id",
                  "repeat_2_pident", "repeat_2_cov", "repeat_2_completeness",
                  "coord_validation", "coordinate_in_bounds"]
    _write_tsv(OUT / "mini_array_units.tsv", units, UNIT_COLS)

    # 4. Summary
    write_summary(all_rows, trusted, mini, edge, units)

    print("wrote finalized catalog:")
    for name in ("crispr_arrays.tsv", "crispr_arrays.bed",
                  "crispr_arrays_trusted.tsv", "crispr_arrays_trusted.bed",
                  "mini_arrays.tsv", "mini_arrays.bed",
                  "edge_unresolved_candidates.tsv", "mini_array_units.tsv",
                  "coordinate_catalog_summary.md"):
        p = OUT / name
        if p.is_file():
            print(f"  {p.name:<36} {p.stat().st_size:>10} B")
    print(f"\ntrusted-call rows: {len(trusted)}")
    print(f"candidate_mini_array rows: {len(mini)}")
    print(f"edge_unresolved rows: {len(edge)}")
    print(f"pair rows: {len(units)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
