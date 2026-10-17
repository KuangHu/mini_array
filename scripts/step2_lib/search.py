"""BLAST the Step-1 HQ repeat queries against all 1000 staged contigs.

Expensive steps (building the subject DB, running blastn) are run via the
patched env `envs/identify`. This module just drives subprocess calls and
parses the tabular output.
"""
from __future__ import annotations

import csv
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen", "sstrand",
]


@dataclass
class HitRow:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int   # 1-based
    qend: int     # 1-based inclusive
    sstart: int   # 1-based; cf. sstrand ("plus" → sstart<=send, "minus" → sstart>=send)
    send: int
    evalue: float
    bitscore: float
    qlen: int
    slen: int
    sstrand: str  # "plus" | "minus"


def build_subject_db(staged_dir: Path, out_db: Path, blastdb_env_bin: Path) -> None:
    """Concatenate every staged FASTA into one combined DB file and makeblastdb."""
    tmp = out_db.with_suffix(".combined.fna")
    out_db.parent.mkdir(parents=True, exist_ok=True)
    with tmp.open("w") as outfh:
        for f in sorted(staged_dir.glob("*.fna")):
            with f.open() as inp:
                for line in inp:
                    outfh.write(line)
    cmd = [str(blastdb_env_bin / "makeblastdb"),
           "-in", str(tmp), "-dbtype", "nucl",
           "-out", str(out_db), "-title", out_db.name]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"makeblastdb failed: {r.stderr}")
    # keep combined FASTA for other tools; could remove to save disk
    # tmp.unlink(missing_ok=True)


def run_blastn(query_fa: Path, db_prefix: Path, out_tsv: Path,
                blastdb_env_bin: Path, *, threads: int,
                word_size: int, evalue: float, dust: str,
                outfmt: str) -> None:
    cmd = [
        str(blastdb_env_bin / "blastn"),
        "-query", str(query_fa),
        "-db", str(db_prefix),
        "-out", str(out_tsv),
        "-num_threads", str(threads),
        "-word_size", str(word_size),
        "-evalue", str(evalue),
        "-dust", str(dust),
        "-outfmt", outfmt,
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"blastn failed: {r.stderr}")


def read_hits(tsv: Path) -> list[HitRow]:
    rows = []
    with tsv.open() as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < len(BLAST_COLUMNS):
                continue
            rows.append(HitRow(
                qseqid=parts[0],
                sseqid=parts[1],
                pident=float(parts[2]),
                length=int(parts[3]),
                mismatch=int(parts[4]),
                gapopen=int(parts[5]),
                qstart=int(parts[6]),
                qend=int(parts[7]),
                sstart=int(parts[8]),
                send=int(parts[9]),
                evalue=float(parts[10]),
                bitscore=float(parts[11]),
                qlen=int(parts[12]),
                slen=int(parts[13]),
                sstrand=parts[14],
            ))
    return rows


def hit_start0_end0(h: HitRow) -> tuple[int, int]:
    """Return (0-based, exclusive) interval on the subject's forward strand."""
    if h.sstrand == "plus":
        return h.sstart - 1, h.send
    else:
        # BLAST reports sstart > send for the minus strand; flip for forward coords.
        return h.send - 1, h.sstart


def hit_query_coverage(h: HitRow) -> float:
    return h.length / h.qlen if h.qlen else 0.0


def hit_passes(h: HitRow, min_identity: float, min_cov: float) -> bool:
    return (h.pident / 100.0) >= min_identity and hit_query_coverage(h) >= min_cov
