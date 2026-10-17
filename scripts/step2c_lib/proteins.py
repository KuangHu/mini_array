"""Phase K3: recover ALL CDS by re-running prodigal on staged FASTA.

cctyper's protein FASTA is not persisted (--keep_tmp off). We re-run
prodigal under envs/identify (which already has prodigal 2.6.3). Mode
matches Step 1 (`isolate → single`, `MAG → meta`).

Output per genome: `<out>/proteins/<guid>.faa` + `<out>/proteins/<guid>.gff`.
Deterministic (prodigal has no RNG).
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path


def prodigal_mode(source_type: str) -> str:
    return "single" if source_type.strip().lower() == "isolate" else "meta"


def recover_proteome(staged_fna: Path, out_faa: Path, out_gff: Path,
                      mode: str, prodigal_bin: str) -> None:
    """Run prodigal on `staged_fna`. Returns silently on success."""
    out_faa.parent.mkdir(parents=True, exist_ok=True)
    if out_faa.is_file() and out_faa.stat().st_size > 0:
        return
    cmd = [prodigal_bin, "-i", str(staged_fna),
           "-p", mode, "-a", str(out_faa), "-f", "gff",
           "-o", str(out_gff), "-q"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"prodigal failed on {staged_fna}: {r.stderr[-500:]}")


def read_fasta(path: Path) -> list[tuple[str, str]]:
    rows = []
    cur_h, cur_s = None, []
    with path.open() as fh:
        for ln in fh:
            if ln.startswith(">"):
                if cur_h is not None:
                    rows.append((cur_h, "".join(cur_s)))
                cur_h = ln[1:].rstrip("\n")
                cur_s = []
            else:
                cur_s.append(ln.strip())
        if cur_h is not None:
            rows.append((cur_h, "".join(cur_s)))
    return rows


def parse_prodigal_header(header: str) -> dict:
    """Prodigal faa headers are like
      '<contig_id>_<N> # start # end # strand # ID=...'
    Return {contig, pos, start, end, strand, length}."""
    parts = header.split(" # ")
    first = parts[0]
    # first = contig_<X>_N
    if "_" in first:
        # Find trailing _N
        idx = first.rfind("_")
        contig = first[:idx]
        pos = int(first[idx + 1:])
    else:
        contig = first
        pos = 0
    out = {"contig": contig, "pos": pos}
    if len(parts) >= 4:
        try:
            out["start"] = int(parts[1])
            out["end"] = int(parts[2])
            out["strand"] = int(parts[3])
        except ValueError:
            pass
    return out
