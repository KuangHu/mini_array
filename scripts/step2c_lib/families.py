"""Phase K4: annotation-independent protein-family clustering.

Pipeline:
  1. Concatenate all neighborhood proteins into one FASTA.
  2. `makeblastdb -dbtype prot`.
  3. `blastp` all-vs-all with permissive settings.
  4. Build an undirected graph: edge if (ident ≥ 30% AND cov ≥ 70%).
  5. Connected components → families.
  6. Each component gets a representative = longest protein.

Important: proteins without any neighbors (singletons) stay as their own
family — do NOT collapse unknowns into one group.
"""
from __future__ import annotations

import subprocess
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class ProteinFamily:
    family_id: str
    representative: str
    members: list[str] = field(default_factory=list)
    n_members: int = 0


def build_db(faa: Path, db_prefix: Path, env_bin: Path) -> None:
    cmd = [str(env_bin / "makeblastdb"), "-in", str(faa), "-dbtype", "prot",
           "-out", str(db_prefix), "-title", db_prefix.name]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"makeblastdb (prot) failed: {r.stderr[-500:]}")


def blastp(query: Path, db_prefix: Path, out_tsv: Path, env_bin: Path, *,
            threads: int, evalue: float = 1e-5) -> None:
    cmd = [str(env_bin / "blastp"), "-query", str(query), "-db", str(db_prefix),
           "-out", str(out_tsv), "-num_threads", str(threads),
           "-evalue", str(evalue),
           "-outfmt",
           "6 qseqid sseqid pident length mismatch gapopen qstart qend "
           "sstart send evalue bitscore qlen slen"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"blastp failed: {r.stderr[-500:]}")


def cluster_connected_components(blastp_tsv: Path, protein_ids: list[str], *,
                                   min_pident: float, min_cov: float
                                   ) -> list[ProteinFamily]:
    """Union-find clustering of protein_ids using edges from the blastp TSV."""
    parent = {p: p for p in protein_ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    with blastp_tsv.open() as fh:
        for ln in fh:
            if not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 14:
                continue
            qid, sid = p[0], p[1]
            if qid == sid:
                continue
            try:
                pident = float(p[2])
                length = int(p[3])
                qlen = int(p[12])
                slen = int(p[13])
            except ValueError:
                continue
            min_len = min(qlen, slen)
            cov = length / min_len if min_len else 0
            if pident / 100.0 >= min_pident and cov >= min_cov:
                union(qid, sid)

    # Group
    groups: dict[str, list[str]] = {}
    for p in protein_ids:
        r = find(p)
        groups.setdefault(r, []).append(p)

    families = []
    for i, (rep, members) in enumerate(sorted(groups.items())):
        families.append(ProteinFamily(
            family_id=f"pfam_{i+1:05d}",
            representative=rep,
            members=sorted(members),
            n_members=len(members),
        ))
    return families
