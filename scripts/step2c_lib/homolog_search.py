"""Phase K5: find homologs of each protein family across all 1000 proteomes.

Approach:
  - Select each family's representative sequence.
  - BLAST all representatives against a concatenated DB of all 1000
    recovered proteomes.
  - For each HSP: record subject contig+pos, subject genome, identity,
    coverage, bitscore.
"""
from __future__ import annotations

import subprocess
from pathlib import Path


HOMOLOG_COLS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen",
                 "evalue", "bitscore"]


def blastp_homologs(reps_faa: Path, db_prefix: Path, out_tsv: Path,
                     env_bin: Path, *, threads: int,
                     evalue: float = 1e-5) -> None:
    cmd = [str(env_bin / "blastp"), "-query", str(reps_faa),
           "-db", str(db_prefix), "-out", str(out_tsv),
           "-num_threads", str(threads), "-evalue", str(evalue),
           "-max_target_seqs", "5000",
           "-outfmt",
           "6 qseqid sseqid pident length qlen slen evalue bitscore"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"blastp homolog search failed: {r.stderr[-500:]}")


def filter_hits(tsv: Path, min_pident: float, min_cov: float) -> list[dict]:
    out = []
    with tsv.open() as fh:
        for ln in fh:
            if not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            try:
                pident = float(p[2])
                length = int(p[3])
                qlen = int(p[4])
                slen = int(p[5])
                evalue = float(p[6])
                bitscore = float(p[7])
            except ValueError:
                continue
            cov = length / qlen if qlen else 0
            if pident / 100.0 < min_pident or cov < min_cov:
                continue
            out.append({
                "qseqid": p[0], "sseqid": p[1],
                "pident": pident, "length": length,
                "qlen": qlen, "slen": slen,
                "evalue": evalue, "bitscore": bitscore,
            })
    return out
