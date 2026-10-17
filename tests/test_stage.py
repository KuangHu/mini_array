"""Staging unit tests covering step1.md §12 cases that are local to prepare_inputs."""
from __future__ import annotations

import csv
import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from step1_lib import stage as stage_mod  # noqa: E402


def _write_fasta(path: Path, records: list[tuple[str, str]], gzipped=False) -> Path:
    body = "".join(f">{h}\n{s}\n" for h, s in records)
    if gzipped:
        with gzip.open(path, "wt") as fh:
            fh.write(body)
    else:
        path.write_text(body)
    return path


def test_sanitize_removes_unsafe_chars():
    assert stage_mod.sanitize("GCA_001.1") == "GCA_001_1"
    assert stage_mod.sanitize("foo bar|baz") == "foo_bar_baz"
    assert stage_mod.sanitize("ok-stuff_123") == "ok-stuff_123"


def test_stage_one_rewrites_header_and_hashes(tmp_path: Path):
    fa = _write_fasta(tmp_path / "g.fna", [
        ("contig_1 some descr", "ACGTacgtNnRy"),
        ("contig_2|other", "AAAAACCCCCGGGGGTTTTT"),
    ])
    staged = tmp_path / "staged"
    rows, out, in_sha, out_sha = stage_mod.stage_one("GCA_1.1", fa, staged)

    assert out.is_file()
    text = out.read_text().splitlines()
    # Header is `>c_<20 hex> contig` — the trailing " contig" is required so
    # upstream CRISPRidentify's header.split(" ")[0] drops the newline.
    assert text[0].startswith(">c_"), text[0]
    first_token = text[0].split(" ")[0]
    assert len(first_token) == 1 + len("c_") + 20
    assert text[0].endswith(" contig")
    # upper-cased, IUPAC preserved
    joined_seq = text[1]
    assert joined_seq == "ACGTACGTNNRY"

    assert [r.record_index for r in rows] == [0, 1]
    assert rows[0].original_id == "contig_1"
    assert rows[0].original_header == "contig_1 some descr"
    assert rows[0].length == 12
    assert rows[0].topology == "unknown"
    assert in_sha != out_sha  # staging rewrites, so hashes differ by design


def test_stage_one_rejects_illegal_chars(tmp_path: Path):
    fa = _write_fasta(tmp_path / "bad.fna", [("c1", "ACGT1234")])
    with pytest.raises(ValueError, match="illegal"):
        stage_mod.stage_one("GCA_bad.1", fa, tmp_path / "staged")


def test_contig_uid_is_deterministic(tmp_path: Path):
    fa = _write_fasta(tmp_path / "g.fna", [("c", "ACGTACGTACGT")])
    r1, *_ = stage_mod.stage_one("GCA_x.1", fa, tmp_path / "s1")
    r2, *_ = stage_mod.stage_one("GCA_x.1", fa, tmp_path / "s2")
    assert r1[0].contig_uid == r2[0].contig_uid
    assert r1[0].sequence_sha256 == r2[0].sequence_sha256


def test_header_collision_across_genomes_is_safe(tmp_path: Path):
    # Two different genomes both call their first contig "contig_1".
    g1 = _write_fasta(tmp_path / "g1.fna", [("contig_1", "ACGTACGTAC")])
    g2 = _write_fasta(tmp_path / "g2.fna", [("contig_1", "GGCCGGCCGG")])
    staged = tmp_path / "staged"
    r1, out1, *_ = stage_mod.stage_one("GCA_a.1", g1, staged)
    r2, out2, *_ = stage_mod.stage_one("GCA_b.1", g2, staged)
    assert out1 != out2  # different genome_uid → different staged files
    assert r1[0].contig_uid != r2[0].contig_uid
    # Both record the same original header, which is OK — contig_uid disambiguates.
    assert r1[0].original_id == r2[0].original_id == "contig_1"


def test_stage_one_handles_gzip_input(tmp_path: Path):
    fa = tmp_path / "g.fna.gz"
    _write_fasta(fa, [("c1", "ACGT")], gzipped=True)
    rows, *_ = stage_mod.stage_one("GCA_g.1", fa, tmp_path / "staged")
    assert rows[0].length == 4


def test_stage_manifest_writes_contigs_tsv(tmp_path: Path):
    f1 = _write_fasta(tmp_path / "g1.fna", [("c1", "ACGTACGT")])
    f2 = _write_fasta(tmp_path / "g2.fna", [("c1", "AAACCCGGGTTT"), ("c2", "ACGT")])
    manifest = [
        {"genome_id": "GCA_a.1", "fasta_path": str(f1)},
        {"genome_id": "GCA_b.1", "fasta_path": str(f2)},
    ]
    staged = tmp_path / "staged"
    contigs = tmp_path / "contigs.tsv"
    summary = stage_mod.stage_manifest(manifest, staged, contigs)
    assert summary["staged_count"] == 2
    assert summary["failed_count"] == 0
    assert summary["contig_rows"] == 3

    with contigs.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 3
    assert {r["genome_id"] for r in rows} == {"GCA_a.1", "GCA_b.1"}
    for r in rows:
        for col in stage_mod.CONTIGS_COLS:
            assert col in r


def test_stage_manifest_captures_failure_without_aborting(tmp_path: Path):
    good = _write_fasta(tmp_path / "g.fna", [("c1", "ACGT")])
    bad = _write_fasta(tmp_path / "b.fna", [("c1", "ACGT123")])
    manifest = [
        {"genome_id": "GCA_a.1", "fasta_path": str(good)},
        {"genome_id": "GCA_b.1", "fasta_path": str(bad)},
    ]
    summary = stage_mod.stage_manifest(manifest, tmp_path / "s", tmp_path / "c.tsv")
    assert summary["staged_count"] == 1
    assert summary["failed_count"] == 1
    assert "GCA_b.1" == summary["failed"][0][0]
