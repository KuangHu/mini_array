"""Pull ±10 kb CDS + annotations around candidates (reuses cctyper outputs).

cctyper's `genes.tab` has every prodigal ORF (Contig, Start, End, Strand, Pos).
`hmmer.tab` has only the subset of ORFs that hit a Cas-related HMM. For each
candidate we emit:
  - a flat list of ORFs within ±neighborhood_bp (position-ordered)
  - for each ORF: pos, strand, start0, end0, length, and ALL hmm hits
  - a "signature" of HMM families seen in the window (for module clustering)
  - a flag for whether any Cas operon elsewhere on the genome exists
    (via cas_operons.tab, cas_operons_putative.tab, cas_operons_orphan.tab)
"""
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


def load_cas_operons_for_genome(genome_uid: str, raw_cctyper_root: Path) -> list[dict]:
    out = []
    for name, conf in (("cas_operons.tab", "high"),
                        ("cas_operons_putative.tab", "putative"),
                        ("cas_operons_orphan.tab", "orphan")):
        p = raw_cctyper_root / genome_uid / name
        if not p.is_file():
            continue
        with p.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                out.append({"contig": r.get("Contig"),
                             "operon": r.get("Operon"),
                             "prediction": r.get("Prediction"),
                             "best_type": r.get("Best_type"),
                             "confidence": conf,
                             "start1": int(r.get("Start", 0) or 0),
                             "end1": int(r.get("End", 0) or 0)})
    return out


def neighborhood_for_candidate(contig_uid: str, candidate_start0: int,
                                candidate_end0: int, *,
                                neighborhood_bp: int,
                                orfs: list[dict],
                                hmm_by_pos: dict[int, list[dict]]
                                ) -> dict:
    window_s = candidate_start0 - neighborhood_bp
    window_e = candidate_end0 + neighborhood_bp
    local = []
    hmm_families: list[str] = []
    for o in orfs:
        o_s0 = o["Start"] - 1
        o_e0 = o["End"]
        if o_e0 < window_s or o_s0 > window_e:
            continue
        hmm_hits = hmm_by_pos.get(o["Pos"], [])
        local.append({
            "pos": o["Pos"],
            "strand": o["Strand"],
            "start0": o_s0,
            "end0": o_e0,
            "length": o_e0 - o_s0,
            "distance_bp": max(0,
                               min(abs(o_s0 - candidate_end0),
                                   abs(candidate_start0 - o_e0))),
            "hmm_hits": hmm_hits,
        })
        hmm_families.extend(h["hmm"] for h in hmm_hits)
    local.sort(key=lambda x: x["start0"])
    sig = _signature(hmm_families)
    return {
        "window_start0": window_s,
        "window_end0": window_e,
        "n_orfs_in_window": len(local),
        "n_orfs_with_hmm": sum(1 for L in local if L["hmm_hits"]),
        "hmm_families_in_window": sorted(set(hmm_families)),
        "hmm_signature": sig,
        "orfs": local,
    }


def _signature(hmm_families: list[str]) -> str:
    """Sorted, dedup list of HMM families with multiplicity: 'A×2+B+C'."""
    if not hmm_families:
        return "no_hmm_hits"
    c = Counter(hmm_families)
    parts = []
    for name in sorted(c):
        n = c[name]
        parts.append(name if n == 1 else f"{name}×{n}")
    return "+".join(parts)
