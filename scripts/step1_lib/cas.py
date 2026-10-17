"""Parse CRISPRCasTyper output and link to parsed arrays (step1.md §8).

Reads these cctyper tables per genome:
  cas_operons.tab           — Cas-based operons with Prediction/Best_type
  cas_operons_putative.tab  — same shape, lower-confidence predictions
  cas_operons_orphan.tab    — single-gene operons (we treat these as putative only)

Associates each parsed array with any operon on the SAME contig whose
`gap_bp <= max_local_gap_bp`. Does NOT cross contigs even if both live
in the same assembly (step1.md §8 rule).

Emits:
  array_cas_links.tsv — one row per (array, operon) pair + a synthetic
                        no_proximal_operon row for arrays with no local hit.
"""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Operon:
    contig_uid: str
    operon_id: str              # e.g. "c_<...>@1"
    start1: int                 # cctyper is 1-based inclusive on both ends
    end1: int
    prediction: str             # "I-C", "III-D", "Unknown", "False", ...
    best_type: str              # same scale, used for conflict detection
    best_score: float
    confidence: str             # "high" | "putative" | "orphan"
    genes: str                  # as reported (list repr), kept opaque
    source_table: str


def _read_cctyper_table(path: Path, confidence: str) -> list[Operon]:
    """cctyper writes TSVs with header + Python list reprs — parse opaquely."""
    if not path.is_file():
        return []
    out: list[Operon] = []
    with path.open(newline="") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        if "Contig" not in idx or "Operon" not in idx:
            return []
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < len(header):
                continue
            try:
                start1 = int(cols[idx["Start"]])
                end1 = int(cols[idx["End"]])
            except (ValueError, KeyError):
                continue
            try:
                best_score = float(cols[idx.get("Best_score", -1)])
            except (ValueError, IndexError):
                best_score = float("nan")
            out.append(Operon(
                contig_uid=cols[idx["Contig"]],
                operon_id=cols[idx["Operon"]],
                start1=start1, end1=end1,
                prediction=cols[idx.get("Prediction", -1)] if "Prediction" in idx else "",
                best_type=cols[idx.get("Best_type", -1)] if "Best_type" in idx else "",
                best_score=best_score,
                confidence=confidence,
                genes=cols[idx["Genes"]] if "Genes" in idx else "",
                source_table=path.name,
            ))
    return out


def parse_cctyper_operons(cctyper_dir: Path) -> list[Operon]:
    operons: list[Operon] = []
    operons += _read_cctyper_table(cctyper_dir / "cas_operons.tab", "high")
    operons += _read_cctyper_table(cctyper_dir / "cas_operons_putative.tab", "putative")
    operons += _read_cctyper_table(cctyper_dir / "cas_operons_orphan.tab", "orphan")
    # The three tables can overlap on the same operon_id; keep the first (high)
    # occurrence per operon_id.
    seen: set[str] = set()
    deduped: list[Operon] = []
    for o in operons:
        if o.operon_id in seen:
            continue
        seen.add(o.operon_id)
        deduped.append(o)
    return deduped


def _gap_bp(a_start0: int, a_end0: int, op_start1: int, op_end1: int) -> int:
    """Minimum gap between array [a_start0, a_end0) and operon [op_start1, op_end1]."""
    # Convert operon to 0-based half-open to make math uniform.
    op_s0 = op_start1 - 1
    op_e0 = op_end1  # 1-based inclusive end + 1 → 0-based exclusive
    return max(0, a_start0 - op_e0, op_s0 - a_end0)


def _normalize_prediction(pred: str) -> str:
    p = (pred or "").strip()
    if not p or p == "False":
        return "unknown"
    return p


def _is_high_confidence(op: Operon) -> bool:
    return op.confidence == "high" and _normalize_prediction(op.prediction) != "unknown"


def link_arrays_to_operons(arrays: list, operons_by_genome: dict[str, list[Operon]],
                           max_gap_bp: int) -> list[dict]:
    """Return per-array link rows per step1.md §8 (table at end of section)."""
    rows: list[dict] = []
    by_contig: dict[tuple[str, str], list[Operon]] = {}
    for genome_id, ops in operons_by_genome.items():
        for op in ops:
            by_contig.setdefault((genome_id, op.contig_uid), []).append(op)

    for a in arrays:
        key = (a.genome_id, a.contig_uid)
        candidates: list[tuple[Operon, int]] = []
        for op in by_contig.get(key, []):
            gap = _gap_bp(a.start0, a.end0, op.start1, op.end1)
            if gap <= max_gap_bp:
                candidates.append((op, gap))

        if not candidates:
            rows.append({
                "array_id": a.array_id,
                "genome_id": a.genome_id,
                "contig_uid": a.contig_uid,
                "cas_link_status": "no_proximal_operon",
                "operon_id": "",
                "operon_confidence": "",
                "subtype_prediction": "",
                "gap_bp": "",
                "candidate_operons": "[]",
            })
            continue

        high_subtypes = sorted({
            _normalize_prediction(o.prediction)
            for o, _ in candidates
            if _is_high_confidence(o)
        })
        high_subtypes = [s for s in high_subtypes if s != "unknown"]
        has_putative = any(o.confidence != "high" for o, _ in candidates)
        has_only_putative = not any(_is_high_confidence(o) for o, _ in candidates)

        if has_only_putative:
            status = "putative_only"
        elif len(high_subtypes) == 1 and len([o for o, _ in candidates if _is_high_confidence(o)]) == 1:
            status = "proximal_operon"
        elif len(high_subtypes) == 1:
            status = "multiple_operons_same_subtype"
        else:
            status = "ambiguous_local"

        # emit one row per candidate operon so step1.md §8's "保留多对多关系"
        # requirement is respected downstream.
        for op, gap in sorted(candidates, key=lambda t: t[1]):
            rows.append({
                "array_id": a.array_id,
                "genome_id": a.genome_id,
                "contig_uid": a.contig_uid,
                "cas_link_status": status,
                "operon_id": op.operon_id,
                "operon_confidence": op.confidence,
                "subtype_prediction": _normalize_prediction(op.prediction),
                "gap_bp": gap,
                "candidate_operons": json.dumps(
                    sorted({o.operon_id for o, _ in candidates}),
                    separators=(",", ":"),
                ),
            })
    return rows


CAS_LINK_COLS = [
    "array_id", "genome_id", "contig_uid",
    "cas_link_status", "operon_id", "operon_confidence", "subtype_prediction",
    "gap_bp", "candidate_operons",
]


def write_cas_links_tsv(rows: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(".tsv.tmp")
    with tmp.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CAS_LINK_COLS, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    tmp.rename(out_path)


def load_all_operons(project_root: Path, manifest_rows: list[dict]) -> dict:
    """Return (operons_by_genome, missing_list, status_counts)."""
    from . import stage as _stage
    by_genome: dict[str, list[Operon]] = {}
    missing: list[str] = []
    for row in manifest_rows:
        gid = row["genome_id"]
        guid = _stage.sanitize(gid)
        cct_dir = project_root / "raw" / "cctyper" / guid
        if not cct_dir.is_dir():
            missing.append(gid)
            continue
        by_genome[gid] = parse_cctyper_operons(cct_dir)
    return by_genome, missing
