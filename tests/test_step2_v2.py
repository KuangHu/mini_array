"""Regression and synthetic end-to-end tests for v2 scientific semantics."""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from dataclasses import asdict, replace
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import step2
import step2b
import step2c
import build_candidate_catalog as catalog
from step2_lib import pair, search, context
from step2b_lib import classify as review_b
from step2c_lib import structural_review as review_c


def hit(start, *, query="q", span=30, qstart=1, pid=100, strand="plus", slen=10000, **kwargs):
    sstart, send = start + 1, start + span
    if strand == "minus":
        sstart, send = send, sstart
    return search.HitRow(query, "ctg", pid, kwargs.get("length", span), 0,
                         kwargs.get("gapopen", 0), qstart, qstart + span - 1,
                         sstart, send, 1e-5, kwargs.get("bitscore", 100), 30, slen, strand)


def find(hits):
    return pair.find_pairs(hits, anchor_min_cov=.9, anchor_min_id=.9,
                           second_min_cov=.6, second_min_id=.7, spacer_min=15, spacer_max=60)


def table(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def serialized(p, tmp_path):
    dest = tmp_path / "pairs.tsv"
    pair.write_pairs_tsv([p], dest)
    return table(dest)[0]


@pytest.mark.parametrize("strand", ["plus", "minus"])
@pytest.mark.parametrize("short_side", ["left", "right"])
def test_short_copy_on_either_side_and_strand(strand, short_side, tmp_path):
    left = hit(3000, span=21 if short_side == "left" else 30, strand=strand)
    right = hit(3060, span=21 if short_side == "right" else 30, strand=strand)
    result = find([right, left])
    assert len(result) == 1
    p = result[0]
    assert p.pair_kind == "half_repeat"
    assert p.anchor_start0 == (3060 if short_side == "left" else 3000)
    nr = pair.normalized_record(serialized(p, tmp_path))
    assert nr["repeat_1_completeness"] == ("truncated" if short_side == "left" else "full")
    assert nr["repeat_2_completeness"] == ("truncated" if short_side == "right" else "full")
    assert nr["spacer_start0"] == 3000 + left.qend
    assert nr["spacer_end0"] == 3060


def test_reverse_complement_invariant():
    hits = [hit(3000, span=21), hit(3060)]
    reverse = [replace(h, sstart=h.slen-h.sstart+1, send=h.slen-h.send+1, sstrand="minus") for h in hits]
    forward, backward = find(hits)[0], find(reverse)[0]
    assert forward.pair_kind == backward.pair_kind == "half_repeat"
    assert forward.spacer_len == backward.spacer_len
    assert forward.anchor_cov == backward.anchor_cov


def test_cross_query_not_a_pair_and_same_query_alternatives_survive():
    assert find([hit(3000, query="A"), hit(3060, query="B")]) == []
    # A high-scoring hit from another query must not suppress either true pair.
    results = find([hit(3000, query="A", bitscore=80), hit(3060, query="A"),
                    hit(3000, query="B", bitscore=200), hit(3060, query="B")])
    assert len(results) == 2
    assert {p.anchor_qseqid for p in results} == {"A", "B"}
    assert all(p.anchor_qseqid == p.second_qseqid for p in results)


def test_gap_columns_cannot_create_an_anchor():
    gapped = hit(3000, span=24, length=30, gapopen=1)
    assert search.hit_query_coverage(gapped) == .8
    assert not search.hit_passes(gapped, .9, .9)
    assert find([gapped, hit(3060, span=24)]) == []
    assert find([gapped, hit(3060)])[0].pair_kind == "half_repeat"


def test_full_coverage_degenerate_copy_is_mini_array(tmp_path):
    p = find([hit(3000, pid=75), hit(3060)])[0]
    assert p.pair_kind == "mini_array" and p.second_is_full
    row = serialized(p, tmp_path)
    assert row["second_passes_anchor"] == "false"
    nr = pair.normalized_record(row)
    assert nr["repeat_1_completeness"] == "full"
    assert nr["repeat_1_sequence_quality"] == "degenerate"
    assert nr["repeat_2_sequence_quality"] == "conserved"


def test_overlapping_hsp_does_not_hide_full_anchor():
    # A higher-scoring partial overlapping the anchor must not suppress it.
    results = find([hit(3000, span=24, bitscore=200), hit(3002), hit(3062, span=21)])
    assert len(results) == 1
    assert results[0].anchor_start0 == 3002


@pytest.mark.parametrize("start,excluded", [(99, True), (100, False), (9810, False), (9811, True)])
def test_edge_boundaries(start, excluded):
    p = find([hit(start), hit(start+60)])[0]
    assert p.edge_excluded is excluded


@pytest.mark.parametrize("review", [review_b.classify, review_c.classify])
def test_edge_always_wins_even_known_array_context(review):
    ov = {"nearby_identify_arrays": [{"array_id": "arr", "category_raw": "Bona-fide", "overlap_bp": 90,
                                      "start0": 0, "end0": 1000}],
          "nearby_identify_units": [{"inside_candidate": True, "array_id": "arr", "unit_kind": "repeat"}],
          "nearby_cctyper_crisprs": []}
    assert review(1500, 1590, ov, 1500, 2000, 2000)["structural_class"] == "edge_excluded"
    assert context.classify_locus("ctg", 10, 90, {"ctg": [(0, 1000)]}, {}, {"ctg": 10000}) == "contig_edge"


def test_spacer_units_are_not_repeat_evidence():
    ov = {"nearby_identify_arrays": [], "nearby_cctyper_crisprs": [],
          "nearby_identify_units": [{"inside_candidate": True, "array_id": "a", "unit_kind": "spacer"}] * 2}
    assert review_c.classify(3000, 3090, ov, 2000, 2000, 2000)["structural_class"] == "candidate_mini_array"


def test_legacy_or_corrupted_pair_records_rejected(tmp_path):
    with pytest.raises(ValueError, match="Legacy"):
        pair.normalized_record({})
    row = serialized(find([hit(3000), hit(3060)])[0], tmp_path)
    with pytest.raises(ValueError, match="same repeat query"):
        pair.normalized_record(dict(row, second_qseqid="other"))
    with pytest.raises(ValueError, match="coverage"):
        pair.normalized_record(dict(row, second_cov="1.1"))
    with pytest.raises(ValueError, match="classification"):
        pair.normalized_record(dict(row, pair_kind="half_repeat"))


def fixture_project(tmp_path, hits, background=True):
    cfg = yaml.safe_load((Path(__file__).resolve().parents[1] / "config/step2.yaml").read_text())
    cfg["project_root"], cfg["code_root"] = tmp_path, Path(__file__).resolve().parents[1]
    cfg["background"]["n_replicates"] = 1
    cfg["queries_fna"], cfg["strict_queries_fna"] = "queries.fna", "strict.fna"
    cfg["known_array_bed"], cfg["arrays_all_tsv"] = "known.bed", "arrays.tsv"
    for path in ("input/staged", "results/staging", "raw/cctyper", "raw_hits", "results/step1_patched_min_score_0.75"):
        (tmp_path / path).mkdir(parents=True)
    seq = ("ACGTTGCA" * 1250)
    (tmp_path / "input/staged/G.fna").write_text(f">ctg\n{seq}\n")
    (tmp_path / "queries.fna").write_text(">q\n" + "ACG" * 10 + "\n")
    (tmp_path / "strict.fna").write_text(">q\n" + "ACG" * 10 + "\n")
    (tmp_path / "known.bed").write_text("")
    (tmp_path / "arrays.tsv").write_text("coordinate_check_pass\tcategory\tcontig_uid\tstart0\tend0\n")
    (tmp_path / "results/step1_patched_min_score_0.75/array_units.tsv").write_text("contig_uid\tstart0\tend0\n")
    (tmp_path / "results/staging/contigs.1000.tsv").write_text("genome_id\tgenome_uid\tcontig_uid\tlength\nG\tG\tctg\t10000\n")
    (tmp_path / "input/genomes.tsv").write_text("genome_id\tindependent_group_id\tsource_type\nG\tTest\tisolate\n")
    for name in (["real", "bg_rep0"] if background else ["real"]):
        (tmp_path / f"raw_hits/hits_{name}.tsv").write_text(
            "".join("\t".join(str(asdict(h)[c]) for c in search.BLAST_COLUMNS) + "\n" for h in hits))
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump({k: str(v) if isinstance(v, Path) else v for k, v in cfg.items()}))
    return cfg, config_path


def run_review(cfg, config_path, monkeypatch):
    import _slurm_gate
    # Tests run tiny fixtures, not production workloads or fake SLURM jobs.
    monkeypatch.setattr(_slurm_gate, "require_slurm", lambda: "test-fixture")
    monkeypatch.setattr(step2c, "require_slurm", lambda: "test-fixture")
    monkeypatch.delenv("PROJECT_ROOT", raising=False)
    assert step2b.main(["--config", str(config_path)]) == 0
    assert step2c.main(["--config", str(config_path), "--phases", "dedup", "review"]) == 0


def test_end_to_end_replay_review_and_split_exports(tmp_path, monkeypatch):
    hits = [hit(3000), hit(3060), hit(3005, query="qhalf", span=21, qstart=10), hit(3060, query="qhalf"),
            hit(5000), hit(5060, span=21), hit(1500), hit(1560), hit(50), hit(110),
            hit(7000, pid=75, strand="minus"), hit(7060, strand="minus")]
    cfg, path = fixture_project(tmp_path, hits)
    step2.do_run(cfg, hits_dir=tmp_path / "raw_hits")
    stats = json.loads((tmp_path / cfg["output_dir"] / "summary.json").read_text())["pair_counts"]
    assert stats["real"] == stats["bg_rep0"]
    assert stats["real"]["edge_excluded"] == 1
    assert stats["real"]["mini_array"] == 3
    assert stats["real"]["half_repeat"] == 2
    run_review(cfg, path, monkeypatch)
    summary = catalog.build(cfg)
    assert summary["candidate_physical_loci_union"] == 3
    assert summary["mini_array_loci"] == summary["half_repeat_loci"] == 2
    assert summary["loci_with_both_kinds"] == 1
    assert summary["candidate_pair_counts"] == {"mini_array": 2, "half_repeat": 2}
    assert summary["review_edge_excluded"] == summary["step2_edge_excluded"] == 1
    out = tmp_path / cfg["catalog_dir"]
    halves = table(out / "half_repeats.pairs.tsv")
    assert any(r["repeat_1_completeness"] == "truncated" for r in halves)
    assert any(r["repeat_2_completeness"] == "truncated" for r in halves)
    for stem in ("all_pairs", "mini_arrays", "half_repeats"):
        fasta = dict(step2._read_fasta(out / f"{stem}.fna"))
        index = table(out / f"{stem}.index.tsv")
        assert len(index) == len(fasta)
        seq = "ACGTTGCA" * 1250
        for r in index:
            dna = fasta[r["record_id"]]
            assert dna == seq[int(r["start0"]):int(r["end0"])]
            assert len(dna) == int(r["length"])
            assert hashlib.sha256(dna.encode()).hexdigest() == r["sha256"]
    with pytest.raises(ValueError, match="empty"):
        step2.do_run(cfg, hits_dir=tmp_path / "raw_hits")
    with pytest.raises(ValueError, match="Config differs"):
        step2.require_v2_run(tmp_path / cfg["output_dir"], dict(cfg, spacer={"min_len": 1, "max_len": 60}))
    pair_file = tmp_path / cfg["output_dir"] / "pairs_real.tsv"
    pair_file.write_text(pair_file.read_text() + "\n")
    with pytest.raises(ValueError, match="output changed"):
        step2.require_v2_run(tmp_path / cfg["output_dir"])


def test_empty_run_still_exports_headers(tmp_path, monkeypatch):
    cfg, path = fixture_project(tmp_path, [], background=False)
    step2.do_run(cfg, hits_dir=tmp_path / "raw_hits", run_background=False)
    run_review(cfg, path, monkeypatch)
    assert catalog.build(cfg)["candidate_physical_loci_union"] == 0
    assert (tmp_path / cfg["catalog_dir"] / "half_repeats.pairs.tsv").read_text().startswith("pair_id\t")


def test_missing_requested_background_fails_closed(tmp_path):
    cfg, _ = fixture_project(tmp_path, [], background=False)
    with pytest.raises(FileNotFoundError):
        step2.do_run(cfg, hits_dir=tmp_path / "raw_hits")
    with pytest.raises(ValueError, match="completed"):
        step2.require_v2_run(tmp_path / cfg["output_dir"])


def test_exact_distance_fallback(monkeypatch):
    from step1_lib import families, qc
    monkeypatch.setattr(families, "edlib", None)
    monkeypatch.setattr(qc, "edlib", None)
    assert families._edit_distance("ACGTACGT", "ACGTTCGT") == 1
    assert families._edit_distance("", "ACG") == 3
    assert families._edit_distance("ACGT", "ACT") == 1
    assert qc._normalized_edit_similarity("ACGTACGT", "ACGTTCGT") == .875
