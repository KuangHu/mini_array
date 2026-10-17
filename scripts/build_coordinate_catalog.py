#!/usr/bin/env python3
"""Build the Step 1/2/2B/2C coordinate catalog with no biological inference.

Outputs (all 0-based, half-open on BED; TSV mirrors keep both raw and
derived columns):

  crispr_arrays.tsv / crispr_arrays.bed      — identify + cctyper arrays
  mini_arrays.tsv / mini_arrays.bed          — Step 2B/2C candidate_mini_array
  mini_array_units.tsv                       — repeat_1, spacer, repeat_2
  coordinate_catalog_summary.md              — provenance + counts + validations

Preserves category, provenance (which caller produced the record), physical
dedup, strand hypotheses, supporting query identifiers, repeat completeness
flag.
"""
from __future__ import annotations

import csv
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _slurm_gate import require_slurm  # noqa: E402

PROJECT = Path("/groups/rubin/projects/kuang/out/CRISPR_array_related")
STEP1 = PROJECT / "results" / "step1_patched_min_score_0.75"
STEP2 = PROJECT / "results" / "step2_pilot"
STEP2B = PROJECT / "results" / "step2b_pilot"
STEP2C = PROJECT / "results" / "step2c_pilot"
OUT = PROJECT / "results" / "coordinate_catalog"


def _read(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _load_staged(contig_uids_needed: set[str],
                 staged_dir: Path, contigs_tsv: Path) -> dict[str, str]:
    """Return {contig_uid: sequence} for only the contigs we need to validate.
    Keeps memory bounded even though we have 211k contigs on disk."""
    want_by_genome: dict[str, set[str]] = defaultdict(set)
    with contigs_tsv.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["contig_uid"] in contig_uids_needed:
                want_by_genome[r["genome_uid"]].add(r["contig_uid"])
    seqs: dict[str, str] = {}
    for guid, cids in want_by_genome.items():
        f = staged_dir / f"{guid}.fna"
        if not f.is_file():
            continue
        cur = None
        parts: list[str] = []
        with f.open() as fh:
            for ln in fh:
                if ln.startswith(">"):
                    if cur is not None and cur in cids:
                        seqs[cur] = "".join(parts)
                    cur = ln[1:].split()[0]
                    parts = []
                else:
                    parts.append(ln.strip())
            if cur is not None and cur in cids:
                seqs[cur] = "".join(parts)
    return seqs


def _validate(contig_seq: str, start0: int, end0: int) -> str:
    """Return 'ok' / 'oob_left' / 'oob_right' / 'degenerate_bases'."""
    L = len(contig_seq)
    if start0 < 0:
        return "oob_left"
    if end0 > L:
        return "oob_right"
    slice_ = contig_seq[start0:end0]
    if any(c not in "ACGT" for c in slice_):
        return "degenerate_bases"
    return "ok"


# ---------------- crispr_arrays ----------------

def build_crispr_arrays(contig_to_meta: dict) -> list[dict]:
    """Union of identify arrays (all categories, coord-pass only) and cctyper
    arrays, with provenance preserved. Does NOT merge records from the two
    callers — each caller's call is one row with distinct `source_tool`."""
    rows: list[dict] = []
    # identify
    for r in _read(STEP1 / "arrays.all.tsv"):
        if r["coordinate_check_pass"] != "true":
            continue
        try:
            s0, e0 = int(r["start0"]), int(r["end0"])
        except ValueError:
            continue
        rows.append({
            "source_tool": "crispridentify",
            "source_id": r["array_id"],
            "genome_id": r["genome_id"],
            "contig_uid": r["contig_uid"],
            "start0": s0,
            "end0": e0,
            "length": e0 - s0,
            "strand_raw": r["strand_raw"],
            "orientation_prediction": r["orientation_prediction"],
            "category_raw": r["category_raw"],
            "qc_tier": r["qc_tier"],
            "score_or_trusted": r["score"],
            "subtype_or_subtype": "",
            "n_repeats": r["n_repeats"],
            "n_spacers": r["n_spacers"],
            "consensus_repeat_forward": r["consensus_forward"],
            "repeat_id": r["repeat_id"],
            "qc_reasons": r["qc_reasons"],
            "provenance": "step1_patched_min_score_0.75/arrays.all.tsv",
        })
    # cctyper (reuse the raw crisprs_all.tab for each genome in raw/cctyper/)
    cctyper_root = PROJECT / "raw" / "cctyper"
    for d in sorted(cctyper_root.iterdir()):
        if not d.is_dir():
            continue
        tab = d / "crisprs_all.tab"
        if not tab.is_file():
            continue
        guid = d.name
        gid = None
        with tab.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                try:
                    s0 = int(r["Start"]) - 1
                    e0 = int(r["End"])
                except (ValueError, KeyError):
                    continue
                if gid is None:
                    # reverse-lookup via contig_to_meta (any contig in this genome)
                    for cu, m in contig_to_meta.items():
                        if m["genome_uid"] == guid:
                            gid = m["genome_id"]
                            break
                rows.append({
                    "source_tool": "crisprcastyper",
                    "source_id": r["CRISPR"],
                    "genome_id": gid or guid,
                    "contig_uid": r["Contig"],
                    "start0": s0,
                    "end0": e0,
                    "length": e0 - s0,
                    "strand_raw": "",
                    "orientation_prediction": "",
                    "category_raw": "Trusted" if r.get("Trusted") == "True" else "Untrusted",
                    "qc_tier": "",
                    "score_or_trusted": r.get("Trusted", ""),
                    "subtype_or_subtype": r.get("Subtype", ""),
                    "n_repeats": r.get("N_repeats", ""),
                    "n_spacers": "",
                    "consensus_repeat_forward": r.get("Consensus_repeat", ""),
                    "repeat_id": "",
                    "qc_reasons": "",
                    "provenance": f"raw/cctyper/{guid}/crisprs_all.tab",
                })
    rows.sort(key=lambda r: (r["genome_id"], r["contig_uid"], int(r["start0"])))
    return rows


# ---------------- mini_arrays (candidate_mini_array only + edge_unresolved) ----

def build_mini_arrays() -> tuple[list[dict], list[dict]]:
    """Return (mini_rows, edge_rows).

    `mini_rows`: physical loci whose structural_class == 'candidate_mini_array'.
    `edge_rows`: structural_class == 'edge_unresolved'.
    Both carry physical_locus_id, span, strand hypotheses, supporting queries.
    """
    phys = _read(STEP2C / "physical_candidate_evidence.tsv")
    phys_dedup = _read(STEP2C / "physical_loci.tsv")
    orig_evidence = _read(STEP2B / "candidate_evidence.tsv")
    # original locus_id → row
    by_locus_id = {r["locus_id"]: r for r in orig_evidence}
    # physical_locus_id → list[original locus_id]
    orig_by_phys = {r["physical_locus_id"]: r["original_locus_ids"].split(",")
                    for r in phys_dedup}

    mini_rows, edge_rows = [], []
    for r in phys:
        orig = orig_by_phys.get(r["physical_locus_id"], [])
        orig_rows = [by_locus_id[o] for o in orig if o in by_locus_id]
        queries = set()
        strict = False
        for o in orig_rows:
            for q in (o.get("query_set", "") or "").split(","):
                if q.strip():
                    queries.add(q.strip())
            if o.get("strict_support") == "true":
                strict = True
        out = {
            "physical_locus_id": r["physical_locus_id"],
            "genome_id": r["genome_id"],
            "genome_uid": r["genome_uid"],
            "contig_uid": r["contig_uid"],
            "contig_length": r["contig_length"],
            "start0": int(r["span_start0"]),
            "end0": int(r["span_end0"]),
            "length": int(r["span_end0"]) - int(r["span_start0"]),
            "strand_hypotheses": r["strand_hypotheses"],
            "n_original_loci": r["n_original_loci"],
            "original_locus_ids": orig_by_phys.get(r["physical_locus_id"], []),
            "query_repeat_ids": ",".join(sorted(queries)),
            "strict_0.90_support": "true" if strict else "false",
            "flank_left_bp": r["flank_left_bp"],
            "flank_right_bp": r["flank_right_bp"],
            "structural_class": r["structural_class"],
            "structural_reasons": r["structural_reasons"],
            "n_units_inside": r["n_units_inside"],
            "n_units_in_flank": r["n_units_in_flank"],
            "evidence_arrays_with_units": r["evidence_arrays_with_units"],
            "provenance": "step2c_pilot/physical_candidate_evidence.tsv",
        }
        out["original_locus_ids"] = ",".join(out["original_locus_ids"])
        if r["structural_class"] == "candidate_mini_array":
            mini_rows.append(out)
        elif r["structural_class"] == "edge_unresolved":
            edge_rows.append(out)
    return mini_rows, edge_rows


# ---------------- mini_array_units (repeat_1, spacer, repeat_2) --------------

def build_mini_array_units(mini_rows: list[dict], edge_rows: list[dict]
                            ) -> list[dict]:
    """For every physical locus in `mini_rows` or `edge_rows`, pull the Step 2
    pairs that formed each supporting original locus, and emit one row per
    (physical_locus_id, pair) with repeat_1 + spacer + repeat_2 coords.

    Completeness of repeat_2:
      - 'full'       : second_cov ≥ 0.90 AND second_pident ≥ 0.90
      - 'partial'    : second_cov ∈ [0.60, 0.90) with pident ≥ 0.70
      - 'degenerate' : otherwise (anything that still passed secondary thresholds)
    """
    # Load Step 2 pairs — index by contig + anchor_start0 since the same
    # contig can host many pairs.
    pairs_tsv = STEP2 / "pairs_real.tsv"
    pairs_rows = _read(pairs_tsv)
    pairs_by_contig = defaultdict(list)
    for p in pairs_rows:
        pairs_by_contig[p["contig"]].append(p)

    units: list[dict] = []
    for r in mini_rows + edge_rows:
        cuid = r["contig_uid"]
        s0 = int(r["start0"])
        e0 = int(r["end0"])
        # Find pairs whose footprint is inside this physical locus
        for p in pairs_by_contig.get(cuid, []):
            a0 = int(p["anchor_start0"])
            a1 = int(p["anchor_end0"])
            s2_0 = int(p["second_start0"])
            s2_1 = int(p["second_end0"])
            foot_start = min(a0, s2_0)
            foot_end = max(a1, s2_1)
            if foot_start < s0 or foot_end > e0:
                continue
            # Completeness of repeat_2
            cov = float(p["second_cov"])
            pid = float(p["second_pident"])
            if cov >= 0.90 and pid >= 90.0:
                completeness = "full"
            elif cov >= 0.60 and pid >= 70.0:
                completeness = "partial"
            else:
                completeness = "degenerate"
            # Normalise so repeat_1 is always upstream (lower start0) regardless
            # of which was anchor/second
            if a0 <= s2_0:
                r1_s, r1_e, r1_q = a0, a1, p["anchor_qseqid"]
                r2_s, r2_e, r2_q = s2_0, s2_1, p["second_qseqid"]
                r2_pident = float(p["second_pident"])
                r2_cov = float(p["second_cov"])
            else:
                r1_s, r1_e, r1_q = s2_0, s2_1, p["second_qseqid"]
                r2_s, r2_e, r2_q = a0, a1, p["anchor_qseqid"]
                r2_pident = float(p["anchor_pident"])
                r2_cov = float(p["anchor_cov"])
            sp_s = r1_e
            sp_e = r2_s
            units.append({
                "physical_locus_id": r["physical_locus_id"],
                "genome_id": r["genome_id"],
                "contig_uid": cuid,
                "strand_pair": p["strand"],
                "repeat_1_start0": r1_s,
                "repeat_1_end0": r1_e,
                "repeat_1_query_id": r1_q,
                "spacer_start0": sp_s,
                "spacer_end0": sp_e,
                "spacer_len": sp_e - sp_s,
                "repeat_2_start0": r2_s,
                "repeat_2_end0": r2_e,
                "repeat_2_query_id": r2_q,
                "repeat_2_pident": r2_pident,
                "repeat_2_cov": r2_cov,
                "repeat_2_completeness": completeness,
            })
    units.sort(key=lambda u: (u["physical_locus_id"],
                              u["repeat_1_start0"]))
    return units


# ---------------- BED emitters ----------------

def _write_bed(path: Path, rows: list[dict], *,
                contig_key: str, start_key: str, end_key: str,
                name_key: str, score_key: str | None, strand_key: str | None
                ) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for r in rows:
            name = r.get(name_key, ".")
            score = r.get(score_key, "0") if score_key else "0"
            # BED score must be integer 0-1000; clamp or 0 if non-numeric
            try:
                sc = max(0, min(1000, int(round(float(score) * 1000))
                                 if 0 <= float(score) <= 1 else int(float(score))))
            except (ValueError, TypeError):
                sc = 0
            strand = r.get(strand_key, ".") if strand_key else "."
            if "|" in strand:          # strand_hypotheses may be "minus|plus"
                strand = "."
            if strand in ("plus", "forward", "+"):
                strand = "+"
            elif strand in ("minus", "reversed", "-"):
                strand = "-"
            elif strand not in ("+", "-"):
                strand = "."
            fh.write("\t".join([r[contig_key],
                                 str(r[start_key]),
                                 str(r[end_key]),
                                 str(name),
                                 str(sc),
                                 strand]) + "\n")


def _write_tsv(path: Path, rows: list[dict], cols: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


# ---------------- main ----------------

def main() -> int:
    require_slurm()
    OUT.mkdir(parents=True, exist_ok=True)

    # ------ load contig lookup once ------
    contigs = {}
    with (PROJECT / "results" / "staging" / "contigs.1000.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            contigs[r["contig_uid"]] = r

    print(f"[load] contigs={len(contigs)}")

    # ------ crispr_arrays ------
    print("[K-arrays] assembling identify + cctyper arrays…")
    arr_rows = build_crispr_arrays(contigs)
    print(f"[K-arrays] rows={len(arr_rows)} "
          f"(crispridentify {sum(1 for r in arr_rows if r['source_tool']=='crispridentify')}, "
          f"crisprcastyper {sum(1 for r in arr_rows if r['source_tool']=='crisprcastyper')})")

    # Validate arrays against the staged FASTA (sample — would take minutes
    # for the full 70k-plus rows if sequences loaded in full. Use targeted
    # load: only contigs that host any array.)
    arr_contigs = {r["contig_uid"] for r in arr_rows}
    print(f"[K-arrays] contigs touched by arrays: {len(arr_contigs)}")
    arr_seqs = _load_staged(arr_contigs,
                             PROJECT / "input" / "staged",
                             PROJECT / "results" / "staging" / "contigs.1000.tsv")
    val_counts = Counter()
    for r in arr_rows:
        s = arr_seqs.get(r["contig_uid"])
        if s is None:
            val_counts["seq_missing"] += 1
            r["coord_validation"] = "seq_missing"
            continue
        v = _validate(s, int(r["start0"]), int(r["end0"]))
        val_counts[v] += 1
        r["coord_validation"] = v
    print(f"[K-arrays] validation: {dict(val_counts)}")
    _write_tsv(
        OUT / "crispr_arrays.tsv", arr_rows,
        cols=["source_tool", "source_id", "genome_id", "contig_uid",
              "start0", "end0", "length", "strand_raw",
              "orientation_prediction", "category_raw", "qc_tier",
              "score_or_trusted", "subtype_or_subtype", "n_repeats",
              "n_spacers", "consensus_repeat_forward", "repeat_id",
              "qc_reasons", "provenance", "coord_validation"]
    )

    # BED name = source_id + category; strand from orientation_prediction
    for r in arr_rows:
        r["_name"] = f"{r['source_tool']}:{r['source_id']}:{r['category_raw']}"
        r["_strand"] = r["orientation_prediction"] or "."
    _write_bed(OUT / "crispr_arrays.bed", arr_rows,
                contig_key="contig_uid", start_key="start0", end_key="end0",
                name_key="_name", score_key=None, strand_key="_strand")

    # ------ mini_arrays + edge ------
    print("[K-mini] building mini-arrays + edge_unresolved…")
    mini_rows, edge_rows = build_mini_arrays()
    print(f"[K-mini] candidate_mini_array={len(mini_rows)}, "
          f"edge_unresolved={len(edge_rows)}")

    # Validate against FASTA
    mini_contigs = {r["contig_uid"] for r in mini_rows + edge_rows}
    mini_seqs = _load_staged(mini_contigs,
                              PROJECT / "input" / "staged",
                              PROJECT / "results" / "staging" / "contigs.1000.tsv")
    for r in mini_rows + edge_rows:
        s = mini_seqs.get(r["contig_uid"])
        if s is None:
            r["coord_validation"] = "seq_missing"
            continue
        r["coord_validation"] = _validate(s, int(r["start0"]), int(r["end0"]))

    MINI_COLS = ["physical_locus_id", "genome_id", "genome_uid", "contig_uid",
                 "contig_length", "start0", "end0", "length",
                 "strand_hypotheses", "n_original_loci", "original_locus_ids",
                 "query_repeat_ids", "strict_0.90_support", "flank_left_bp",
                 "flank_right_bp", "structural_class", "structural_reasons",
                 "n_units_inside", "n_units_in_flank",
                 "evidence_arrays_with_units", "provenance",
                 "coord_validation"]
    _write_tsv(OUT / "mini_arrays.tsv", mini_rows, cols=MINI_COLS)
    _write_tsv(OUT / "edge_unresolved_candidates.tsv", edge_rows, cols=MINI_COLS)

    # BED emit for mini-arrays only
    for r in mini_rows:
        r["_name"] = f"{r['physical_locus_id']}:{r['structural_class']}"
        r["_strand"] = (r["strand_hypotheses"] if "|" not in r["strand_hypotheses"]
                        else ".")
    _write_bed(OUT / "mini_arrays.bed", mini_rows,
                contig_key="contig_uid", start_key="start0", end_key="end0",
                name_key="_name", score_key=None, strand_key="_strand")

    # ------ mini_array_units ------
    print("[K-units] extracting repeat_1 + spacer + repeat_2 per pair…")
    unit_rows = build_mini_array_units(mini_rows, edge_rows)
    print(f"[K-units] rows={len(unit_rows)}")

    # Validate repeat coordinates
    needed = {u["contig_uid"] for u in unit_rows}
    unit_seqs = _load_staged(needed,
                              PROJECT / "input" / "staged",
                              PROJECT / "results" / "staging" / "contigs.1000.tsv")
    for u in unit_rows:
        s = unit_seqs.get(u["contig_uid"])
        if s is None:
            u["coord_validation"] = "seq_missing"
            continue
        validations = []
        for pair in (("repeat_1_start0", "repeat_1_end0"),
                      ("spacer_start0", "spacer_end0"),
                      ("repeat_2_start0", "repeat_2_end0")):
            validations.append(_validate(s, int(u[pair[0]]), int(u[pair[1]])))
        u["coord_validation"] = ",".join(validations)
    _write_tsv(
        OUT / "mini_array_units.tsv", unit_rows,
        cols=["physical_locus_id", "genome_id", "contig_uid",
              "strand_pair", "repeat_1_start0", "repeat_1_end0",
              "repeat_1_query_id", "spacer_start0", "spacer_end0",
              "spacer_len", "repeat_2_start0", "repeat_2_end0",
              "repeat_2_query_id", "repeat_2_pident", "repeat_2_cov",
              "repeat_2_completeness", "coord_validation"]
    )

    # ------ summary ------
    cat_counts_arr = Counter(r["category_raw"] for r in arr_rows)
    tool_counts = Counter(r["source_tool"] for r in arr_rows)
    cat_counts_mini = Counter(r["structural_class"] for r in mini_rows + edge_rows)
    completeness_counts = Counter(u["repeat_2_completeness"] for u in unit_rows)
    coord_val_arr = Counter(r["coord_validation"] for r in arr_rows)
    coord_val_mini = Counter(r["coord_validation"] for r in mini_rows + edge_rows)
    coord_val_units = Counter(u["coord_validation"] for u in unit_rows)

    md = []
    md.append("# Coordinate catalog summary\n")
    md.append("Pure coordinate catalog assembled from existing Step 1, Step 2, "
               "Step 2B, Step 2C outputs. No biology. No new DB. No parameter "
               "change. 0-based, half-open on BED.\n")

    md.append("## 1. CRISPR arrays\n")
    md.append(f"- Rows total: **{len(arr_rows)}**")
    md.append("- By source tool:")
    for k, n in tool_counts.most_common():
        md.append(f"  - `{k}`: {n}")
    md.append("")
    md.append("- Identify category distribution:")
    for k, n in cat_counts_arr.most_common():
        md.append(f"  - `{k}`: {n}")
    md.append("")
    md.append("- Coordinate validation (against staged FASTA):")
    for k, n in coord_val_arr.most_common():
        md.append(f"  - `{k}`: {n}")
    md.append("")

    md.append("## 2. Mini-arrays (candidate_mini_array)\n")
    md.append(f"- Rows: **{len(mini_rows)}**")
    md.append(f"- Edge_unresolved (kept separate): {len(edge_rows)}")
    md.append("- Structural class distribution (mini + edge, from Step 2C):")
    for k, n in cat_counts_mini.most_common():
        md.append(f"  - `{k}`: {n}")
    md.append("")
    md.append("- Coordinate validation:")
    for k, n in coord_val_mini.most_common():
        md.append(f"  - `{k}`: {n}")
    md.append("")

    md.append("## 3. Mini-array units (repeat_1, spacer, repeat_2)\n")
    md.append(f"- Rows: **{len(unit_rows)}** across "
               f"{len({u['physical_locus_id'] for u in unit_rows})} physical loci")
    md.append("- repeat_2 completeness distribution:")
    for k, n in completeness_counts.most_common():
        md.append(f"  - `{k}`: {n}")
    md.append("")
    md.append("- Per-component coordinate validation "
               "(three slots: repeat_1 / spacer / repeat_2):")
    for k, n in coord_val_units.most_common():
        md.append(f"  - `{k}`: {n}")
    md.append("")

    md.append("## 4. Provenance\n")
    md.append("| output | source |")
    md.append("|---|---|")
    md.append("| crispr_arrays.{tsv,bed} | step1_patched_min_score_0.75/arrays.all.tsv + raw/cctyper/*/crisprs_all.tab |")
    md.append("| mini_arrays.{tsv,bed} | step2c_pilot/physical_candidate_evidence.tsv (structural_class=candidate_mini_array) |")
    md.append("| edge_unresolved_candidates.tsv | step2c_pilot/physical_candidate_evidence.tsv (structural_class=edge_unresolved) |")
    md.append("| mini_array_units.tsv | step2_pilot/pairs_real.tsv restricted to pairs whose footprint sits inside a mini or edge locus |")
    md.append("")
    md.append("## 5. BED file conventions\n")
    md.append("- 0-based, half-open (standard BED).")
    md.append("- `name` field: `<source_tool>:<source_id>:<category_raw>` for arrays; "
               "`<physical_locus_id>:<structural_class>` for mini-arrays.")
    md.append("- `strand` column: `+` / `-` where orientation is known; `.` when "
               "the record carries multiple strand hypotheses or no orientation.")
    md.append("")
    md.append("## 6. Non-negotiables honoured\n")
    md.append("- No biology. No HMM. No CDS extraction. No protein. No clustering.")
    md.append("- Dedup only merges records at the physical_locus_id level that "
               "represent the same genomic interval (strand hypotheses stay "
               "underneath).")
    md.append("- Ordinary arrays, terminal extensions, ambiguous regions are NOT "
               "folded into the mini-array list — they stay in `crispr_arrays.tsv` "
               "or `edge_unresolved_candidates.tsv` with their category.")
    md.append("- Component coordinates (repeat_1, spacer, repeat_2) retain their "
               "query attribution and the completeness flag of repeat_2.")
    md.append("- Every record carries a `coord_validation` column with the result "
               "of comparing the start/end interval against the staged FASTA.")
    md.append("")

    (OUT / "coordinate_catalog_summary.md").write_text("\n".join(md) + "\n")
    print(f"wrote {OUT}/{{crispr_arrays.tsv, crispr_arrays.bed, "
          f"mini_arrays.tsv, mini_arrays.bed, mini_array_units.tsv, "
          f"edge_unresolved_candidates.tsv, coordinate_catalog_summary.md}}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
