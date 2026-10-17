"""Extract candidate region + flanks from staged FASTA.

For each candidate locus we pull:
  - the candidate interval [start0, end0)
  - up to `flank_bp` on each side, clipped to the real contig length (so a
    candidate near a contig edge gets whatever flank is actually available)

Lookups are driven by `contigs.1000.tsv` which maps `contig_uid → genome_id`
+ which staged FASTA file it lives in.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass
class ContigSeq:
    contig_uid: str
    genome_id: str
    seq: str


def load_contig_lookup(contigs_tsv: Path) -> dict[str, dict]:
    """{contig_uid: {genome_id, genome_uid, length}} — no sequence yet."""
    out = {}
    with contigs_tsv.open(newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["contig_uid"]] = {
                "genome_id": r["genome_id"],
                "genome_uid": r["genome_uid"],
                "length": int(r["length"]),
            }
    return out


def _read_staged_fasta(path: Path) -> dict[str, str]:
    seqs = {}
    cur = None
    parts = []
    with path.open() as fh:
        for ln in fh:
            if ln.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(parts)
                cur = ln[1:].split()[0]
                parts = []
            else:
                parts.append(ln.strip())
        if cur is not None:
            seqs[cur] = "".join(parts)
    return seqs


class StagedGenomeCache:
    """Lazy per-genome FASTA cache so we don't reload a staged file more than
    once even when we visit many candidates on the same genome."""

    def __init__(self, staged_dir: Path):
        self.staged_dir = staged_dir
        self._cache: dict[str, dict[str, str]] = {}

    def get_contig_seq(self, genome_uid: str, contig_uid: str) -> str:
        if genome_uid not in self._cache:
            self._cache[genome_uid] = _read_staged_fasta(
                self.staged_dir / f"{genome_uid}.fna"
            )
        return self._cache[genome_uid][contig_uid]


def extract_with_flanks(contig_seq: str, start0: int, end0: int,
                         flank_bp: int) -> dict:
    L = len(contig_seq)
    left0 = max(0, start0 - flank_bp)
    right1 = min(L, end0 + flank_bp)
    return {
        "contig_len": L,
        "window_start0": left0,
        "window_end0": right1,
        "left_flank_bp": start0 - left0,
        "right_flank_bp": right1 - end0,
        "candidate_seq": contig_seq[start0:end0],
        "window_seq": contig_seq[left0:right1],
    }
