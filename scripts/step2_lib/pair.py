"""Coherent same-query pairs, with either genomic side allowed as anchor.

Coverage and sequence identity are separate attributes. A mini_array pair has
two full-coverage copies; a half_repeat pair has one full anchor and a shortened
copy (not necessarily exactly half-length). Edge pairs are audit records only.
"""
from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from pathlib import Path
from . import search

METHOD_VERSION = "same-query-bidirectional-v2"


def completeness(coverage: float, full_cov: float = 0.90) -> str:
    return "full" if coverage >= full_cov else "truncated"


def sequence_quality(pident: float, full_identity: float = 0.90) -> str:
    return "conserved" if pident / 100.0 >= full_identity else "degenerate"


@dataclass
class RepeatPair:
    contig: str
    strand: str
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
    anchor_is_strict: bool
    second_is_full: bool  # length completeness only; identity is separate
    anchor_qstart: int = 0
    anchor_qend: int = 0
    anchor_qlen: int = 0
    anchor_gapopen: int = 0
    second_qstart: int = 0
    second_qend: int = 0
    second_qlen: int = 0
    second_gapopen: int = 0
    contig_length: int = 0
    min_flank_bp: int = 100
    full_cov_threshold: float = 0.90
    full_identity_threshold: float = 0.90

    @property
    def pair_kind(self) -> str:
        return "mini_array" if self.second_is_full else "half_repeat"

    @property
    def edge_excluded(self) -> bool:
        s = min(self.anchor_start0, self.second_start0)
        e = max(self.anchor_end0, self.second_end0)
        return s < self.min_flank_bp or self.contig_length - e < self.min_flank_bp


def classify_hits(hits, anchor_min_cov, anchor_min_id, second_min_cov, second_min_id):
    return ([h for h in hits if search.hit_passes(h, anchor_min_id, anchor_min_cov)],
            [h for h in hits if search.hit_passes(h, second_min_id, second_min_cov)])


def find_pairs(hits: list[search.HitRow], *, anchor_min_cov: float,
               anchor_min_id: float, second_min_cov: float, second_min_id: float,
               spacer_min: int, spacer_max: int, overlap_merge_bp: int = 0,
               min_flank_bp: int = 100) -> list[RepeatPair]:
    """Keep query attributions until after pairing; emit edge pairs for audit.

    Only identical genomic hit footprints within one query/strand are collapsed.
    Distinct overlapping HSP boundaries are retained so a stronger partial HSP
    cannot hide another usable full anchor. overlap_merge_bp is a deprecated
    compatibility argument: distance-based pre-pair suppression is disabled.
    """
    if not (0 < second_min_cov <= anchor_min_cov <= 1
            and 0 < second_min_id <= anchor_min_id <= 1
            and 0 <= spacer_min <= spacer_max and min_flank_bp >= 0):
        raise ValueError("Invalid pairing thresholds")
    groups = {}
    for h in hits:
        if search.hit_passes(h, second_min_id, second_min_cov):
            s, e = search.hit_start0_end0(h)
            if not (0 <= s < e <= h.slen) or h.sstrand not in ("plus", "minus"):
                raise ValueError(f"Invalid subject coordinates: {h}")
            groups.setdefault((h.sseqid, h.sstrand, h.qseqid), []).append(h)

    def is_anchor(h):
        return search.hit_passes(h, anchor_min_id, anchor_min_cov)

    def rank(h):
        return (is_anchor(h), h.bitscore, search.hit_query_coverage(h), h.pident,
                -h.gapopen, -h.qstart, -h.qend)

    pairs = []
    for (contig, strand, query), group in sorted(groups.items()):
        if len({h.slen for h in group}) != 1:
            raise ValueError(f"Inconsistent contig lengths: {contig}")
        if len({h.qlen for h in group}) != 1:
            raise ValueError(f"Inconsistent query lengths: {query}")
        unique = {}
        for h in group:
            coords = search.hit_start0_end0(h)
            if coords not in unique or rank(h) > rank(unique[coords]):
                unique[coords] = h
        ordered = sorted(unique.items())
        for i, ((left_s, left_e), left) in enumerate(ordered):
            for j in range(i + 1, len(ordered)):
                (right_s, right_e), right = ordered[j]
                gap = right_s - left_e
                if gap > spacer_max:
                    break
                if gap < spacer_min or not (is_anchor(left) or is_anchor(right)):
                    continue
                # Both full: deterministic left anchor; otherwise full hit wins.
                anchor, second = (left, right) if is_anchor(left) else (right, left)
                a0, a1 = search.hit_start0_end0(anchor)
                b0, b1 = search.hit_start0_end0(second)
                pairs.append(RepeatPair(
                    contig, strand, query, anchor.pident, search.hit_query_coverage(anchor),
                    a0, a1, query, second.pident, search.hit_query_coverage(second), b0, b1,
                    gap, left_e, right_s, True,
                    search.hit_query_coverage(second) >= anchor_min_cov,
                    anchor.qstart, anchor.qend, anchor.qlen, anchor.gapopen,
                    second.qstart, second.qend, second.qlen, second.gapopen,
                    anchor.slen, min_flank_bp, anchor_min_cov, anchor_min_id,
                ))
    return pairs


PAIR_COLS = list(RepeatPair.__dataclass_fields__) + [
    "method_version", "pair_kind", "second_completeness", "second_sequence_quality",
    "second_passes_anchor", "edge_excluded",
]


def write_pairs_tsv(pairs: list[RepeatPair], out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=PAIR_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for p in pairs:
            row = asdict(p)
            row.update(method_version=METHOD_VERSION, pair_kind=p.pair_kind,
                       second_completeness=completeness(p.second_cov, p.full_cov_threshold),
                       second_sequence_quality=sequence_quality(p.second_pident, p.full_identity_threshold),
                       second_passes_anchor=(p.second_is_full and p.second_pident / 100 >= p.full_identity_threshold),
                       edge_excluded=p.edge_excluded)
            w.writerow({k: str(v).lower() if isinstance(v, bool) else v for k, v in row.items()})


def normalized_record(row: dict) -> dict:
    """Validate a v2 row and label genomic-left/right repeats consistently.

    Old exports lack query coordinates. Reprocess raw hits rather than relabel
    old results as corrected output.
    """
    if row.get("method_version") != METHOD_VERSION:
        raise ValueError("Legacy pair table: reprocess raw BLAST hits with Step 2 v2")
    if row["anchor_qseqid"] != row["second_qseqid"]:
        raise ValueError("Pair must use the same repeat query")
    if int(row["anchor_qlen"]) != int(row["second_qlen"]):
        raise ValueError("Same query must have the same length")
    if (float(row["anchor_cov"]) < float(row["full_cov_threshold"])
            or float(row["anchor_pident"]) / 100 < float(row["full_identity_threshold"])):
        raise ValueError("Pair has no valid full anchor")
    if row["strand"] not in ("plus", "minus"):
        raise ValueError("Invalid alignment strand")
    order = sorted(("anchor", "second"), key=lambda x: int(row[x + "_start0"]))
    out = {"method_version": METHOD_VERSION}
    for i, role in enumerate(order, 1):
        prefix = f"repeat_{i}"
        qstart, qend, qlen = (int(row[role + "_" + k]) for k in ("qstart", "qend", "qlen"))
        if not 1 <= min(qstart, qend) <= max(qstart, qend) <= qlen:
            raise ValueError("Invalid query coordinates")
        cov = (abs(qend - qstart) + 1) / qlen
        if abs(cov - float(row[role + "_cov"])) > 1e-9:
            raise ValueError("Stored coverage differs from query span")
        pid = float(row[role + "_pident"])
        if not 0 <= int(row[role + "_start0"]) < int(row[role + "_end0"]) <= int(row["contig_length"]):
            raise ValueError("Invalid subject coordinates")
        out.update({prefix + "_role_in_pair": role,
                    prefix + "_query_id": row[role + "_qseqid"],
                    prefix + "_start0": int(row[role + "_start0"]),
                    prefix + "_end0": int(row[role + "_end0"]),
                    prefix + "_pident": pid, prefix + "_cov": cov,
                    prefix + "_qstart": qstart, prefix + "_qend": qend,
                    prefix + "_qlen": qlen, prefix + "_gapopen": int(row[role + "_gapopen"]),
                    prefix + "_completeness": completeness(cov, float(row["full_cov_threshold"])),
                    prefix + "_sequence_quality": sequence_quality(pid, float(row["full_identity_threshold"]))})
    out["spacer_start0"] = out["repeat_1_end0"]
    out["spacer_end0"] = out["repeat_2_start0"]
    out["spacer_len"] = out["spacer_end0"] - out["spacer_start0"]
    if (out["spacer_len"] != int(row["spacer_len"]) or out["spacer_len"] < 0
            or out["spacer_start0"] != int(row["spacer_start0"])
            or out["spacer_end0"] != int(row["spacer_end0"])):
        raise ValueError("Inconsistent spacer coordinates")
    out["pair_kind"] = "half_repeat" if any(out[f"repeat_{i}_completeness"] == "truncated" for i in (1, 2)) else "mini_array"
    out["edge_excluded"] = (out["repeat_1_start0"] < int(row["min_flank_bp"])
                            or int(row["contig_length"]) - out["repeat_2_end0"] < int(row["min_flank_bp"]))
    if (out["pair_kind"] != row["pair_kind"]
            or str(out["edge_excluded"]).lower() != str(row["edge_excluded"]).lower()):
        raise ValueError("Stored classification disagrees with coordinates")
    return out
