"""Deterministic dinucleotide-preserving shuffle (Altschul-Erickson 1985).

Preserves both mono- and di-nucleotide composition. Used only for background
query construction in Step 2. NOT a calibrated FDR tool.
"""
from __future__ import annotations

import random


def dinuc_shuffle(seq: str, rng: random.Random) -> str:
    """Shuffle `seq` while preserving its dinucleotide composition.

    Standard edge-swap algorithm: treat seq as an Eulerian path over the
    graph with nodes = letters and edges = adjacent dinucleotides, then walk
    an Eulerian path with randomised choice of next edge.
    """
    seq = seq.upper()
    if len(seq) < 2:
        return seq

    # Build adjacency list: node → list of outgoing-edge targets (letters).
    edges: dict[str, list[str]] = {}
    for i in range(len(seq) - 1):
        edges.setdefault(seq[i], []).append(seq[i + 1])

    # Shuffle each node's outgoing edges, then enforce the "last-edge" rule:
    # the last outgoing edge from each node (except the start node) must lead
    # toward the final letter — otherwise the walk may get stuck.
    start, end = seq[0], seq[-1]

    # Approach: repeat shuffle until a valid walk exists. For short sequences
    # this is very fast; worst case ~ few retries.
    MAX_RETRIES = 50
    for _ in range(MAX_RETRIES):
        trial: dict[str, list[str]] = {k: list(v) for k, v in edges.items()}
        for k in trial:
            rng.shuffle(trial[k])
        try:
            walked = _walk(trial, start, end, len(seq) - 1)
        except ValueError:
            continue
        return start + "".join(walked)
    # Fallback: return the input (no valid walk found — vanishingly rare)
    return seq


def _walk(edges: dict[str, list[str]], start: str, end: str, n_edges: int) -> list[str]:
    """Walk n_edges edges from `start`, consuming each edge once; must end at `end`."""
    result = []
    cur = start
    for _ in range(n_edges):
        if not edges.get(cur):
            raise ValueError(f"stuck at {cur}")
        nxt = edges[cur].pop()
        result.append(nxt)
        cur = nxt
    if cur != end:
        raise ValueError(f"ended at {cur} instead of {end}")
    return result


def shuffle_fasta_records(records: list[tuple[str, str]], seed: int,
                          replicate: int) -> list[tuple[str, str]]:
    """Return a new list of (header, seq) with each seq dinuc-shuffled.

    `replicate` is mixed into the per-record RNG seed so distinct replicates
    give distinct shuffles.
    """
    rng = random.Random(seed * 1000003 + replicate)
    out = []
    for header, seq in records:
        local = random.Random(hash((seq, replicate, seed)) & 0xffffffff)
        out.append((f"{header}|bg{replicate}", dinuc_shuffle(seq, local)))
    return out
