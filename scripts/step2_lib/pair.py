"""Anchor + second-repeat pairing for Step-2.

Rules per user instruction:
  - Full-length repeat anchor (coverage, identity thresholds — "anchor" config).
  - Second repeat may be complete, truncated, or degenerate (secondary config).
  - Spacer length window measured between the aligned boundaries of the two
    repeats on the SAME contig, SAME strand. Window 15-60 nt by default.
  - Both strands scanned; alignment orientation retained raw. Canonical query
    direction is NOT treated as transcription direction.
"""
from __future__ import annotations

import csv
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from . import search


@dataclass
class RepeatPair:
    contig: str
    strand: str          # "plus" / "minus"
    anchor_qseqid: str
    anchor_pident: float
    anchor_cov: float
    anchor_start0: int
    anchor_end0: int
    second_qseqid: str
    second_pident: float
    second_cov: float
    second_start0: int
    second_end0: int
    spacer_len: int
    spacer_start0: int
    spacer_end0: int
    anchor_is_strict: bool   # whether anchor also passed anchor thresholds (vs secondary)
    second_is_full: bool     # whether second ALSO passed anchor thresholds


def classify_hits(hits: list[search.HitRow],
                   anchor_min_cov: float, anchor_min_id: float,
                   second_min_cov: float, second_min_id: float
                   ) -> tuple[list[search.HitRow], list[search.HitRow]]:
    """Return (anchors, secondary_candidates). A hit is a secondary if it
    passes the (looser) secondary thresholds — anchors are a subset."""
    anchors, secondary = [], []
    for h in hits:
        if search.hit_passes(h, second_min_id, second_min_cov):
            secondary.append(h)
        if search.hit_passes(h, anchor_min_id, anchor_min_cov):
            anchors.append(h)
    return anchors, secondary


def _h_passes_anchor(h: search.HitRow, min_cov: float, min_id: float) -> bool:
    return search.hit_passes(h, min_id, min_cov)


def find_pairs(hits: list[search.HitRow], *,
               anchor_min_cov: float, anchor_min_id: float,
               second_min_cov: float, second_min_id: float,
               spacer_min: int, spacer_max: int,
               overlap_merge_bp: int) -> list[RepeatPair]:
    """Group hits by (contig, strand); within each group, find anchor-second pairs."""
    # bucket
    by_group: dict[tuple[str, str], list[search.HitRow]] = {}
    for h in hits:
        if not search.hit_passes(h, second_min_id, second_min_cov):
            continue
        by_group.setdefault((h.sseqid, h.sstrand), []).append(h)

    pairs: list[RepeatPair] = []
    for (contig, strand), group in by_group.items():
        # sort by genomic start on the forward strand
        group_coords = [(search.hit_start0_end0(h), h) for h in group]
        group_coords.sort(key=lambda t: (t[0][0], t[0][1]))

        # merge overlapping hits on the same strand: keep highest bitscore
        merged: list[tuple[tuple[int, int], search.HitRow]] = []
        for (s0, e0), h in group_coords:
            if merged and s0 - merged[-1][0][1] < overlap_merge_bp:
                prev = merged[-1][1]
                if h.bitscore > prev.bitscore:
                    merged[-1] = ((s0, e0), h)
                continue
            merged.append(((s0, e0), h))

        # anchors in this group
        for i, ((s0, e0), hi) in enumerate(merged):
            if not _h_passes_anchor(hi, anchor_min_cov, anchor_min_id):
                continue
            # candidate second repeats: any hit downstream whose spacer gap
            # (distance between hi.end0 and hj.start0) is in [spacer_min, spacer_max]
            for j in range(i + 1, len(merged)):
                (sj0, ej0), hj = merged[j]
                gap = sj0 - e0
                if gap < spacer_min:
                    continue
                if gap > spacer_max:
                    break   # further hits are even farther
                pairs.append(RepeatPair(
                    contig=contig, strand=strand,
                    anchor_qseqid=hi.qseqid,
                    anchor_pident=hi.pident, anchor_cov=search.hit_query_coverage(hi),
                    anchor_start0=s0, anchor_end0=e0,
                    second_qseqid=hj.qseqid,
                    second_pident=hj.pident, second_cov=search.hit_query_coverage(hj),
                    second_start0=sj0, second_end0=ej0,
                    spacer_len=gap, spacer_start0=e0, spacer_end0=sj0,
                    anchor_is_strict=True,
                    second_is_full=_h_passes_anchor(
                        hj, anchor_min_cov, anchor_min_id),
                ))
    return pairs


PAIR_COLS = [
    "contig", "strand",
    "anchor_qseqid", "anchor_pident", "anchor_cov",
    "anchor_start0", "anchor_end0",
    "second_qseqid", "second_pident", "second_cov",
    "second_start0", "second_end0",
    "spacer_len", "spacer_start0", "spacer_end0",
    "second_is_full",
]


def write_pairs_tsv(pairs: list[RepeatPair], out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(PAIR_COLS)
        for p in pairs:
            w.writerow([
                p.contig, p.strand,
                p.anchor_qseqid, round(p.anchor_pident, 2), round(p.anchor_cov, 3),
                p.anchor_start0, p.anchor_end0,
                p.second_qseqid, round(p.second_pident, 2), round(p.second_cov, 3),
                p.second_start0, p.second_end0,
                p.spacer_len, p.spacer_start0, p.spacer_end0,
                "true" if p.second_is_full else "false",
            ])
