"""identify × cctyper locus cross-table (classification only, no scoring).

step1.md spirit: cctyper is NOT ground truth. We emit a per-locus matched
table so you can see where the two callers agree, disagree, or exist alone.
Matching rule: same contig_uid AND non-zero genomic overlap. For a given
array we attach the single best-overlapping cctyper locus (by overlap bp),
and vice versa.
"""
from __future__ import annotations

import csv
from pathlib import Path

from . import stage as _stage


def _read_cctyper_crisprs(path: Path) -> list[dict]:
    """Return rows from cctyper's crisprs_all.tab (1-based inclusive coords)."""
    if not path.is_file():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _overlap_bp(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def build_cross_rows(arrays: list, project_root: Path,
                     manifest_rows: list[dict]) -> list[dict]:
    """Return cross-table rows per step1.md Q (identify vs cctyper)."""
    # Group identify-side arrays by genome_id → contig_uid.
    by_contig: dict[tuple[str, str], list] = {}
    for a in arrays:
        by_contig.setdefault((a.genome_id, a.contig_uid), []).append(a)

    rows: list[dict] = []
    seen_cctyper: set[tuple[str, str, str]] = set()

    for m in manifest_rows:
        gid = m["genome_id"]
        guid = _stage.sanitize(gid)
        crisprs_tab = (project_root / "raw" / "cctyper" / guid /
                       "crisprs_all.tab")
        ct_rows = _read_cctyper_crisprs(crisprs_tab)
        # Build index of cctyper loci on this genome by contig.
        ct_by_contig: dict[str, list[dict]] = {}
        for r in ct_rows:
            ct_by_contig.setdefault(r["Contig"], []).append(r)

        # Identify-side: emit a row per parsed array.
        for a in by_contig.get((gid, a.contig_uid if False else a.contig_uid), []):
            pass  # placeholder; real loop below
        # Correct loop over all identify arrays belonging to this genome:
        for a in [x for x in arrays if x.genome_id == gid]:
            ct_candidates = ct_by_contig.get(a.contig_uid, [])
            best: tuple[int, dict] | None = None
            for r in ct_candidates:
                b0 = int(r["Start"]) - 1  # cctyper 1-based inclusive
                b1 = int(r["End"])         # 0-based exclusive after +1 - 1
                ov = _overlap_bp(a.start0, a.end0, b0, b1)
                if ov == 0:
                    continue
                if best is None or ov > best[0]:
                    best = (ov, r)
            if best is None:
                rows.append({
                    "genome_id": gid, "contig_uid": a.contig_uid,
                    "classification": "identify_only",
                    "identify_array_id": a.array_id,
                    "identify_category": a.category,
                    "identify_score": f"{a.score:.6f}",
                    "identify_n_repeats": a.n_repeats,
                    "identify_start0": a.start0,
                    "identify_end0": a.end0,
                    "cctyper_crispr_id": "",
                    "cctyper_trusted": "",
                    "cctyper_subtype": "",
                    "cctyper_n_repeats": "",
                    "cctyper_start0": "",
                    "cctyper_end0": "",
                    "overlap_bp": 0,
                    "start_delta": "",
                    "end_delta": "",
                })
            else:
                ov, r = best
                b0 = int(r["Start"]) - 1
                b1 = int(r["End"])
                rows.append({
                    "genome_id": gid, "contig_uid": a.contig_uid,
                    "classification": "both",
                    "identify_array_id": a.array_id,
                    "identify_category": a.category,
                    "identify_score": f"{a.score:.6f}",
                    "identify_n_repeats": a.n_repeats,
                    "identify_start0": a.start0,
                    "identify_end0": a.end0,
                    "cctyper_crispr_id": r["CRISPR"],
                    "cctyper_trusted": r.get("Trusted", ""),
                    "cctyper_subtype": r.get("Subtype", ""),
                    "cctyper_n_repeats": r.get("N_repeats", ""),
                    "cctyper_start0": b0,
                    "cctyper_end0": b1,
                    "overlap_bp": ov,
                    "start_delta": a.start0 - b0,
                    "end_delta": a.end0 - b1,
                })
                seen_cctyper.add((gid, r["Contig"], r["CRISPR"]))

        # cctyper-only rows: any cctyper CRISPR with no identify-side match.
        for r in ct_rows:
            key = (gid, r["Contig"], r["CRISPR"])
            if key in seen_cctyper:
                continue
            b0 = int(r["Start"]) - 1
            b1 = int(r["End"])
            # Confirm there's no overlapping identify call on this contig
            overlap = False
            for a in by_contig.get((gid, r["Contig"]), []):
                if _overlap_bp(a.start0, a.end0, b0, b1) > 0:
                    overlap = True
                    break
            if overlap:
                continue
            rows.append({
                "genome_id": gid, "contig_uid": r["Contig"],
                "classification": "cctyper_only",
                "identify_array_id": "",
                "identify_category": "",
                "identify_score": "",
                "identify_n_repeats": "",
                "identify_start0": "",
                "identify_end0": "",
                "cctyper_crispr_id": r["CRISPR"],
                "cctyper_trusted": r.get("Trusted", ""),
                "cctyper_subtype": r.get("Subtype", ""),
                "cctyper_n_repeats": r.get("N_repeats", ""),
                "cctyper_start0": b0,
                "cctyper_end0": b1,
                "overlap_bp": 0,
                "start_delta": "",
                "end_delta": "",
            })
    return rows


CROSS_COLS = [
    "genome_id", "contig_uid", "classification",
    "identify_array_id", "identify_category", "identify_score",
    "identify_n_repeats", "identify_start0", "identify_end0",
    "cctyper_crispr_id", "cctyper_trusted", "cctyper_subtype",
    "cctyper_n_repeats", "cctyper_start0", "cctyper_end0",
    "overlap_bp", "start_delta", "end_delta",
]


def write_tsv(rows: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(".tsv.tmp")
    with tmp.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CROSS_COLS, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    tmp.rename(out_path)
