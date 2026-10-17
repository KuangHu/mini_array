"""Per-genome run_status.tsv (step1.md §9)."""
from __future__ import annotations

import csv
import json
from pathlib import Path

from . import stage as _stage

STATUS_COLS = [
    "genome_id", "genome_uid",
    "identify_status", "identify_arrays", "identify_bona_fide",
    "cctyper_status", "cctyper_operons", "cctyper_high_operons",
    "notes",
]


def _check_identify(raw_dir: Path) -> tuple[str, int, int, str]:
    sj = raw_dir / "Complete_json_summary.json"
    if not sj.is_file():
        return "missing", 0, 0, "Complete_json_summary.json not written"
    try:
        with sj.open() as fh:
            summary = json.load(fh)
    except Exception as exc:
        return "parse_failed", 0, 0, f"json load failed: {exc}"
    total = 0
    bona = 0
    for _, blob in summary.items():
        for cat, items in blob.get("Arrays", {}).items():
            total += len(items or {})
            if cat == "Bona-fide":
                bona += len(items or {})
    if total == 0:
        return "success_zero_arrays", 0, 0, "no arrays detected"
    return "success", total, bona, ""


def _check_cctyper(raw_dir: Path) -> tuple[str, int, int, str]:
    args = raw_dir / "arguments.tab"
    if not args.is_file():
        return "missing", 0, 0, "cctyper arguments.tab absent"
    ops = raw_dir / "cas_operons.tab"
    put = raw_dir / "cas_operons_putative.tab"
    orph = raw_dir / "cas_operons_orphan.tab"
    def _count(p: Path) -> int:
        if not p.is_file():
            return 0
        with p.open() as fh:
            return sum(1 for _ in fh) - 1
    total = max(0, _count(ops)) + max(0, _count(put)) + max(0, _count(orph))
    high = max(0, _count(ops))
    status = "success" if total > 0 else "success_zero_operons"
    return status, total, high, ""


def build_rows(project_root: Path, manifest_rows: list[dict]) -> list[dict]:
    rows: list[dict] = []
    for m in manifest_rows:
        gid = m["genome_id"]
        guid = _stage.sanitize(gid)
        id_status, id_total, id_bona, id_note = _check_identify(
            project_root / "raw" / "identify" / guid
        )
        cc_status, cc_total, cc_high, cc_note = _check_cctyper(
            project_root / "raw" / "cctyper" / guid
        )
        notes = "; ".join(x for x in [id_note, cc_note] if x)
        rows.append({
            "genome_id": gid, "genome_uid": guid,
            "identify_status": id_status, "identify_arrays": id_total,
            "identify_bona_fide": id_bona,
            "cctyper_status": cc_status, "cctyper_operons": cc_total,
            "cctyper_high_operons": cc_high,
            "notes": notes,
        })
    return rows


def write_tsv(rows: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(".tsv.tmp")
    with tmp.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=STATUS_COLS, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    tmp.rename(out_path)
