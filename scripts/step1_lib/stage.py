"""Stage FASTAs from the manifest into input/staged/ with safe headers.

For each manifest row:
  1. Open the FASTA (plain or gzip).
  2. Uppercase sequences, validate characters against the IUPAC nucleotide set.
  3. Assign a filesystem-safe `genome_uid` (sanitized from `genome_id`), and a
     globally unique `contig_uid` derived from genome_id + record index +
     sequence SHA256 (step1.md §1).
  4. Rewrite each contig under the header `>contig_uid` into
     `input/staged/<genome_uid>.fna`, 80 cols/line.
  5. Record both source and staged file SHA256 so repeated staging is a no-op.

Returns per-contig rows for aggregation into `contigs.tsv`.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import re
from collections.abc import Iterator
from dataclasses import asdict, dataclass
from pathlib import Path

SAFE_RE = re.compile(r"[^A-Za-z0-9_-]")

# IUPAC nucleotides + gap. U is kept (RNA input would be rare but shouldn't be
# silently rewritten per step1.md §1).
ALLOWED_BASES = set("ACGTRYSWKMBDHVNU-")

CONTIG_UID_PREFIX = "c_"
CONTIG_UID_HEX_LEN = 20


@dataclass
class ContigRow:
    genome_id: str
    genome_uid: str
    contig_uid: str
    original_id: str
    original_header: str
    record_index: int
    length: int
    sequence_sha256: str
    topology: str


def sanitize(name: str) -> str:
    """Replace any char outside [A-Za-z0-9_-] with underscore."""
    return SAFE_RE.sub("_", name)


def _open_fasta(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def _iter_records(path: Path) -> Iterator[tuple[int, str, str]]:
    """Yield (record_index, header_without_gt, sequence_concatenated_upper)."""
    header: str | None = None
    parts: list[str] = []
    idx = -1
    with _open_fasta(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if header is not None:
                    yield idx, header, "".join(parts).upper()
                idx += 1
                header = line[1:]
                parts = []
            else:
                parts.append(line.strip())
        if header is not None:
            yield idx, header, "".join(parts).upper()


def _validate_sequence(seq: str, genome_id: str, record_index: int) -> None:
    bad = set(seq) - ALLOWED_BASES
    if bad:
        sample = sorted(bad)[:5]
        raise ValueError(
            f"{genome_id} contig {record_index}: illegal characters {sample}"
        )


def _contig_uid(genome_id: str, record_index: int, seq_sha256: str) -> str:
    key = f"{genome_id}|{record_index}|{seq_sha256}".encode()
    return CONTIG_UID_PREFIX + hashlib.sha256(key).hexdigest()[:CONTIG_UID_HEX_LEN]


def stage_one(genome_id: str, fasta_path: Path, staged_dir: Path,
              line_width: int = 80) -> tuple[list[ContigRow], Path, str, str]:
    """Stage one manifest row into `<staged_dir>/<genome_uid>.fna`.

    Returns `(contig_rows, staged_path, input_sha256, staged_sha256)`.
    """
    staged_dir.mkdir(parents=True, exist_ok=True)
    input_sha = sha256_file(fasta_path)
    genome_uid = sanitize(genome_id)
    out_path = staged_dir / f"{genome_uid}.fna"
    tmp = out_path.with_suffix(".fna.tmp")

    rows: list[ContigRow] = []
    with tmp.open("w") as out:
        for idx, header, seq in _iter_records(fasta_path):
            _validate_sequence(seq, genome_id, idx)
            seq_sha = hashlib.sha256(seq.encode()).hexdigest()
            cid = _contig_uid(genome_id, idx, seq_sha)
            rows.append(ContigRow(
                genome_id=genome_id,
                genome_uid=genome_uid,
                contig_uid=cid,
                original_id=header.split()[0] if header else "",
                original_header=header,
                record_index=idx,
                length=len(seq),
                sequence_sha256=seq_sha,
                topology="unknown",
            ))
            # Append a trailing " contig" so CRISPRidentify's
            # multiline_fasta_handle_python can safely do header.split(" ")[0]
            # without carrying the terminating newline into the filename.
            # (components_helpers.py:75 bug — a bare ">cid\n" header ends up
            # producing "cid\n.fa" filenames, which breaks the subsequent
            # os.system("mv ...") command.)
            out.write(f">{cid} contig\n")
            for i in range(0, len(seq), line_width):
                out.write(seq[i:i + line_width] + "\n")

    if not rows:
        tmp.unlink(missing_ok=True)
        raise ValueError(f"{genome_id}: no FASTA records in {fasta_path}")

    tmp.rename(out_path)
    staged_sha = sha256_file(out_path)
    return rows, out_path, input_sha, staged_sha


CONTIGS_COLS = [
    "genome_id", "genome_uid", "contig_uid", "original_id", "original_header",
    "record_index", "length", "sequence_sha256", "topology",
    "input_fasta_path", "input_sha256", "staged_fasta_path", "staged_sha256",
]


def write_contigs_tsv(all_rows: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(".tsv.tmp")
    with tmp.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CONTIGS_COLS, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in all_rows:
            w.writerow(r)
    tmp.rename(out_path)


def stage_manifest(manifest_rows: list[dict], staged_dir: Path,
                   contigs_tsv: Path, limit: int | None = None) -> dict:
    """Stage up to `limit` manifest rows; write contigs.tsv; return summary."""
    all_rows: list[dict] = []
    staged: list[str] = []
    failed: list[tuple[str, str]] = []
    seen_uids: dict[str, str] = {}

    for row in manifest_rows[:limit] if limit else manifest_rows:
        gid = row["genome_id"]
        src = Path(row["fasta_path"])
        try:
            contigs, staged_path, input_sha, staged_sha = stage_one(
                gid, src, staged_dir
            )
        except (OSError, ValueError) as exc:
            failed.append((gid, str(exc)))
            continue

        guid = contigs[0].genome_uid
        prev = seen_uids.get(guid)
        if prev is not None and prev != gid:
            failed.append((
                gid,
                f"genome_uid collision with {prev}: both sanitize to {guid}"
            ))
            staged_path.unlink(missing_ok=True)
            continue
        seen_uids[guid] = gid

        for c in contigs:
            d = asdict(c)
            d["input_fasta_path"] = str(src)
            d["input_sha256"] = input_sha
            d["staged_fasta_path"] = str(staged_path)
            d["staged_sha256"] = staged_sha
            all_rows.append(d)
        staged.append(gid)

    write_contigs_tsv(all_rows, contigs_tsv)
    return {
        "staged_count": len(staged),
        "failed_count": len(failed),
        "failed": failed,
        "contig_rows": len(all_rows),
        "contigs_tsv": str(contigs_tsv),
    }


def read_manifest(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))
