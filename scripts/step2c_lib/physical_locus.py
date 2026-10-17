"""Phase K1: merge candidates representing the same physical interval.

Two candidate loci (A, B) collapse to one `physical_locus_id` iff:
  - same genome_id AND same contig_uid
  - their genomic footprints overlap by ≥ 50% of the shorter one, OR
  - start0/end0 differ by ≤ `merge_bp` (default 10) on both ends
Strand is NOT a merger criterion — strand hypotheses are stored underneath.

Important: we do NOT merge truly distinct nearby loci. Candidates with
disjoint footprints or minimal overlap stay separate.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class PhysicalLocus:
    physical_locus_id: str
    genome_id: str
    contig_uid: str
    span_start0: int
    span_end0: int
    strand_hypotheses: set[str] = field(default_factory=set)
    boundary_variants: list[tuple[int, int, str]] = field(default_factory=list)
    original_locus_ids: list[str] = field(default_factory=list)
    attributions: list[dict] = field(default_factory=list)


def _overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def dedup_candidates(candidates: list[dict], merge_bp: int = 10) -> list[PhysicalLocus]:
    """Candidates must have genome_id, contig_uid, candidate_start0,
    candidate_end0, strand, locus_id — all other fields are preserved
    in the attribution list."""
    by_key: dict[tuple[str, str], list[dict]] = {}
    for c in candidates:
        key = (c["genome_id"], c["contig_uid"])
        by_key.setdefault(key, []).append(c)

    phys: list[PhysicalLocus] = []
    for (gid, cuid), lst in by_key.items():
        lst.sort(key=lambda x: (int(x["candidate_start0"]),
                                int(x["candidate_end0"])))
        cur: PhysicalLocus | None = None
        for c in lst:
            s0 = int(c["candidate_start0"])
            e0 = int(c["candidate_end0"])
            if cur is None:
                phys.append(_new_physical(gid, cuid, s0, e0, c, len(phys)))
                cur = phys[-1]
                continue
            # Decide: merge or new?
            ov = _overlap(cur.span_start0, cur.span_end0, s0, e0)
            shorter = min(cur.span_end0 - cur.span_start0, e0 - s0)
            close_both_ends = (abs(cur.span_start0 - s0) <= merge_bp
                                and abs(cur.span_end0 - e0) <= merge_bp)
            if (shorter > 0 and ov / shorter >= 0.50) or close_both_ends:
                cur.span_start0 = min(cur.span_start0, s0)
                cur.span_end0 = max(cur.span_end0, e0)
                cur.strand_hypotheses.add(c["strand"])
                cur.boundary_variants.append((s0, e0, c["strand"]))
                cur.original_locus_ids.append(c["locus_id"])
                cur.attributions.append(c)
            else:
                phys.append(_new_physical(gid, cuid, s0, e0, c, len(phys)))
                cur = phys[-1]
    return phys


def _new_physical(gid, cuid, s0, e0, c, idx) -> PhysicalLocus:
    L = PhysicalLocus(
        physical_locus_id=f"phys_{idx+1:05d}",
        genome_id=gid, contig_uid=cuid,
        span_start0=s0, span_end0=e0,
    )
    L.strand_hypotheses.add(c["strand"])
    L.boundary_variants.append((s0, e0, c["strand"]))
    L.original_locus_ids.append(c["locus_id"])
    L.attributions.append(c)
    return L
