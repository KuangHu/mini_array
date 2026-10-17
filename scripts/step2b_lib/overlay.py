"""Overlay existing identify + cctyper annotations onto a candidate region.

For every candidate we gather:
  - all identify arrays (any category) on the same contig whose footprint
    intersects or sits within `expand_bp` of the candidate
  - all cctyper crisprs on same contig under the same distance rule
  - all cctyper ORFs within ±flank_bp of the candidate (reused to compute
    "additional repeats within 2 kb" via identify_units)

These overlays feed `classify.structural_class`.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


def load_identify_arrays(arrays_all_tsv: Path) -> dict[str, list[dict]]:
    """{contig_uid: [row…]} with int start0/end0 and category_raw."""
    out: dict[str, list[dict]] = defaultdict(list)
    with arrays_all_tsv.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["coordinate_check_pass"] != "true":
                continue
            try:
                r["start0"] = int(r["start0"])
                r["end0"] = int(r["end0"])
            except ValueError:
                continue
            out[r["contig_uid"]].append(r)
    for v in out.values():
        v.sort(key=lambda x: x["start0"])
    return out


def load_identify_units(units_tsv: Path) -> dict[str, list[dict]]:
    """{contig_uid: [unit…]} — every repeat/spacer with start0/end0."""
    out: dict[str, list[dict]] = defaultdict(list)
    with units_tsv.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                r["start0"] = int(r["start0"])
                r["end0"] = int(r["end0"])
            except ValueError:
                continue
            out[r["contig_uid"]].append(r)
    return out


def load_cctyper_crisprs(root_cctyper: Path) -> dict[tuple[str, str], list[dict]]:
    """{(genome_uid, contig_uid): [crispr…]} using 1-based inclusive from cctyper."""
    out: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for d in root_cctyper.iterdir():
        tab = d / "crisprs_all.tab"
        if not tab.is_file():
            continue
        with tab.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                try:
                    s0 = int(r["Start"]) - 1
                    e0 = int(r["End"])
                except ValueError:
                    continue
                out[(d.name, r["Contig"])].append({**r, "start0": s0, "end0": e0})
    return out


def load_cctyper_genes(genome_uid: str, raw_cctyper_root: Path
                        ) -> dict[str, list[dict]]:
    """{contig_uid: [ORF…]} per genome; call sites iterate per candidate."""
    out: dict[str, list[dict]] = defaultdict(list)
    tab = raw_cctyper_root / genome_uid / "genes.tab"
    if not tab.is_file():
        return out
    with tab.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                r["Start"] = int(r["Start"])
                r["End"] = int(r["End"])
                r["Strand"] = int(r["Strand"])
                r["Pos"] = int(r["Pos"])
            except (ValueError, KeyError):
                continue
            out[r["Contig"]].append(r)
    return out


def load_cctyper_hmmer(genome_uid: str, raw_cctyper_root: Path
                       ) -> dict[str, list[dict]]:
    """{contig_uid: [HMM-hit row…]}. hmmer.tab's `Acc` column is contig_uid."""
    out: dict[str, list[dict]] = defaultdict(list)
    tab = raw_cctyper_root / genome_uid / "hmmer.tab"
    if not tab.is_file():
        return out
    with tab.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                r["start"] = int(r["start"])
                r["end"] = int(r["end"])
                r["Pos"] = int(r["Pos"])
            except (ValueError, KeyError):
                continue
            out[r["Acc"]].append(r)
    return out


def _overlap(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def overlay_candidate(contig_uid: str, start0: int, end0: int, *,
                      identify_arrays: dict[str, list[dict]],
                      identify_units: dict[str, list[dict]],
                      cctyper_crisprs_contig: list[dict],
                      cctyper_genes_contig: list[dict],
                      cctyper_hmm_contig: list[dict],
                      flank_bp: int) -> dict:
    """Return structured overlay facts: nearby identify arrays, cctyper
    arrays, repeat units within flank, ORFs within flank with HMM hits."""
    window_s = start0 - flank_bp
    window_e = end0 + flank_bp

    nearby_arrays = []
    for a in identify_arrays.get(contig_uid, []):
        if _overlap(window_s, window_e, a["start0"], a["end0"]) > 0:
            nearby_arrays.append({
                "array_id": a["array_id"],
                "category_raw": a["category_raw"],
                "score": a["score"],
                "start0": a["start0"],
                "end0": a["end0"],
                "qc_tier": a["qc_tier"],
                "overlap_bp": _overlap(start0, end0, a["start0"], a["end0"]),
                "distance_bp": min(abs(a["start0"] - end0),
                                   abs(a["end0"] - start0)),
            })

    nearby_cctyper = []
    for cr in cctyper_crisprs_contig:
        if _overlap(window_s, window_e, cr["start0"], cr["end0"]) > 0:
            nearby_cctyper.append({
                "cctyper_crispr": cr["CRISPR"],
                "subtype": cr.get("Subtype"),
                "trusted": cr.get("Trusted"),
                "start0": cr["start0"],
                "end0": cr["end0"],
                "n_repeats": cr.get("N_repeats"),
                "overlap_bp": _overlap(start0, end0, cr["start0"], cr["end0"]),
            })

    nearby_units = []
    for u in identify_units.get(contig_uid, []):
        if _overlap(window_s, window_e, u["start0"], u["end0"]) > 0:
            nearby_units.append({
                "array_id": u["array_id"],
                "unit_kind": u["unit_kind"],
                "unit_index": u["unit_index"],
                "start0": u["start0"],
                "end0": u["end0"],
                "inside_candidate": _overlap(start0, end0,
                                              u["start0"], u["end0"]) > 0,
            })

    nearby_orfs = []
    for g in cctyper_genes_contig:
        g_s0 = g["Start"] - 1
        g_e0 = g["End"]
        if _overlap(window_s, window_e, g_s0, g_e0) > 0:
            nearby_orfs.append({
                "pos": g["Pos"],
                "strand": g["Strand"],
                "start0": g_s0,
                "end0": g_e0,
            })

    orf_hmm_by_pos: dict[int, list[dict]] = defaultdict(list)
    for h in cctyper_hmm_contig:
        if _overlap(window_s, window_e, h["start"], h["end"]) > 0:
            orf_hmm_by_pos[h["Pos"]].append({
                "hmm": h["Hmm"],
                "evalue": h["Eval"],
                "score": h["score"],
                "cov_hmm": h.get("Cov_hmm"),
            })
    for o in nearby_orfs:
        o["hmm_hits"] = orf_hmm_by_pos.get(o["pos"], [])

    return {
        "nearby_identify_arrays": nearby_arrays,
        "nearby_cctyper_crisprs": nearby_cctyper,
        "nearby_identify_units": nearby_units,
        "nearby_orfs": nearby_orfs,
        "window_start0": window_s,
        "window_end0": window_e,
    }
