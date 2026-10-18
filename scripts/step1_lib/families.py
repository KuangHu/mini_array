"""NR100 catalog, repeat_sources, and deterministic greedy-centroid families.

Inputs: list of (ParsedArray, qc_tier) from upstream.
Outputs:
  repeats.hq.nr100.fna          — forward-canonical seqs for every HQ repeat_id
  repeat_sources.tsv            — repeat_id → array_id → genome/contig/strand
  repeat_family_members.tsv     — per-NR100-seq cluster assignment
  repeat_families.tsv           — per-family summary

Family rule (step1.md §7):
  d(a,b) = min(global_Levenshtein(a,b), global_Levenshtein(a,RC(b)))
  pass iff (1 - d/max(len)) >= 0.90 AND min_len/max_len >= 0.90
  deterministic greedy by (length desc, canonical seq asc)
  new seq joins the lowest-distance passing centroid; ties by centroid_id
  NO transitive merging
"""
from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

try:
    import edlib
except ImportError:  # pragma: no cover
    edlib = None

from . import parse as parse_mod


def _edit_distance(a: str, b: str) -> int:
    if edlib is not None:
        return edlib.align(a, b, task="distance", mode="NW")["editDistance"]
    # Exact Levenshtein fallback; SequenceMatcher is not an edit distance and
    # can split near-identical repetitive sequences into different families.
    previous = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        current = [i]
        for j, cb in enumerate(b, 1):
            current.append(min(current[-1] + 1, previous[j] + 1,
                               previous[j - 1] + (ca != cb)))
        previous = current
    return previous[-1]


def _distance(a: str, b: str) -> int:
    rc_b = parse_mod.revcomp(b)
    return min(_edit_distance(a, b), _edit_distance(a, rc_b))


def _similarity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    return 1 - _distance(a, b) / max(len(a), len(b))


def _length_ratio(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    return min(len(a), len(b)) / max(len(a), len(b))


def _canonical(seq: str) -> str:
    rc = parse_mod.revcomp(seq)
    return seq if seq <= rc else rc


def _family_id(centroid_canonical: str, algo_version: str = "greedy-centroid-v1") -> str:
    key = f"{algo_version}|{centroid_canonical}".encode()
    return "fam_" + hashlib.sha256(key).hexdigest()[:20]


@dataclass
class RepeatCatalogEntry:
    repeat_id: str
    canonical_seq: str
    forward_consensus: str   # one chosen forward representative
    sources: list[dict] = field(default_factory=list)


def build_nr100(hq_arrays: list) -> dict[str, RepeatCatalogEntry]:
    """Group HQ arrays by canonical repeat_id; preserve all sources."""
    catalog: dict[str, RepeatCatalogEntry] = {}
    for a in hq_arrays:
        rep = a.consensus_forward.upper()
        canon = _canonical(rep)
        rid = parse_mod.repeat_id(rep)
        entry = catalog.get(rid)
        if entry is None:
            entry = RepeatCatalogEntry(
                repeat_id=rid,
                canonical_seq=canon,
                forward_consensus=rep,
            )
            catalog[rid] = entry
        entry.sources.append({
            "repeat_id": rid,
            "array_id": a.array_id,
            "genome_id": a.genome_id,
            "contig_uid": a.contig_uid,
            "source_strand": a.orientation_prediction,
            "consensus_forward": rep,
        })
    return catalog


def cluster_families(catalog: dict[str, RepeatCatalogEntry],
                     min_similarity: float, min_length_ratio: float,
                     ) -> dict[str, dict]:
    """Deterministic greedy-centroid clustering per step1.md §7."""
    entries = list(catalog.values())
    entries.sort(key=lambda e: (-len(e.canonical_seq), e.canonical_seq))

    families: dict[str, dict] = {}
    centroid_order: list[str] = []
    assignments: dict[str, dict] = {}

    for e in entries:
        best: tuple[float, str, int, float, float] | None = None
        for fid in centroid_order:
            cent = families[fid]
            if _length_ratio(cent["canonical_seq"], e.canonical_seq) < min_length_ratio:
                continue
            sim = _similarity(cent["canonical_seq"], e.canonical_seq)
            if sim < min_similarity:
                continue
            dist = _distance(cent["canonical_seq"], e.canonical_seq)
            # Prefer lowest distance, then earliest centroid id for ties.
            key = (dist, fid)
            if best is None or key < (best[2], best[1]):
                best = (sim, fid, dist, sim, _length_ratio(cent["canonical_seq"], e.canonical_seq))

        if best is None:
            fid = _family_id(e.canonical_seq)
            families[fid] = {
                "family_id": fid,
                "representative_repeat_id": e.repeat_id,
                "canonical_seq": e.canonical_seq,
                "forward_consensus": e.forward_consensus,
                "members": [],
            }
            centroid_order.append(fid)
            families[fid]["members"].append(e.repeat_id)
            assignments[e.repeat_id] = {
                "family_id": fid,
                "edit_distance": 0,
                "similarity": 1.0,
                "aligned_orientation": "forward",
            }
        else:
            _, fid, dist, sim, _lr = best
            families[fid]["members"].append(e.repeat_id)
            # Report aligned_orientation as "forward" if direct match wins; else "revcomp"
            direct = _edit_distance(families[fid]["canonical_seq"], e.canonical_seq)
            rc = _edit_distance(families[fid]["canonical_seq"], parse_mod.revcomp(e.canonical_seq))
            aligned = "forward" if direct <= rc else "revcomp"
            assignments[e.repeat_id] = {
                "family_id": fid,
                "edit_distance": dist,
                "similarity": sim,
                "aligned_orientation": aligned,
            }
    return {"families": families, "assignments": assignments,
            "centroid_order": centroid_order}


def write_repeat_fasta(catalog: dict[str, RepeatCatalogEntry], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        for rid in sorted(catalog):
            e = catalog[rid]
            fh.write(f">{rid}\n{e.forward_consensus}\n")


SOURCES_COLS = ["repeat_id", "array_id", "genome_id", "contig_uid",
                "source_strand", "consensus_forward"]
MEMBERS_COLS = ["repeat_id", "family_id", "edit_distance", "similarity",
                "aligned_orientation", "canonical_seq", "forward_consensus"]
FAMILIES_COLS = ["family_id", "representative_repeat_id", "n_members",
                 "representative_canonical_seq"]


def write_tsv(path: Path, cols: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    tmp.rename(path)


def write_all(catalog: dict[str, RepeatCatalogEntry], cluster_result: dict,
              out_dir: Path, fasta_path: Path) -> dict:
    write_repeat_fasta(catalog, fasta_path)
    source_rows: list[dict] = []
    for e in catalog.values():
        source_rows.extend(e.sources)
    write_tsv(out_dir / "repeat_sources.tsv", SOURCES_COLS, source_rows)

    member_rows: list[dict] = []
    for rid, info in cluster_result["assignments"].items():
        e = catalog[rid]
        member_rows.append({
            "repeat_id": rid,
            "family_id": info["family_id"],
            "edit_distance": info["edit_distance"],
            "similarity": round(info["similarity"], 4),
            "aligned_orientation": info["aligned_orientation"],
            "canonical_seq": e.canonical_seq,
            "forward_consensus": e.forward_consensus,
        })
    write_tsv(out_dir / "repeat_family_members.tsv", MEMBERS_COLS, member_rows)

    family_rows: list[dict] = []
    for fid in cluster_result["centroid_order"]:
        f = cluster_result["families"][fid]
        family_rows.append({
            "family_id": fid,
            "representative_repeat_id": f["representative_repeat_id"],
            "n_members": len(f["members"]),
            "representative_canonical_seq": f["canonical_seq"],
        })
    write_tsv(out_dir / "repeat_families.tsv", FAMILIES_COLS, family_rows)

    return {
        "nr100_repeats": len(catalog),
        "families": len(cluster_result["families"]),
        "sources": len(source_rows),
        "fasta_path": str(fasta_path),
    }
