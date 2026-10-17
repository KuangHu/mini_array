"""Unit tests for Step 2 modules."""
from __future__ import annotations

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from step2_lib import context as context_mod  # noqa: E402
from step2_lib import dedup as dedup_mod  # noqa: E402
from step2_lib import dinuc_shuffle  # noqa: E402
from step2_lib import pair as pair_mod  # noqa: E402
from step2_lib import search  # noqa: E402


def test_dinuc_shuffle_preserves_mono_composition():
    rng = random.Random(0)
    s = "ACGTACGTACGTACGTACGTACGTACGTACGT"
    shuffled = dinuc_shuffle.dinuc_shuffle(s, rng)
    assert len(shuffled) == len(s)
    for base in "ACGT":
        assert shuffled.count(base) == s.count(base)


def test_dinuc_shuffle_preserves_dinuc_composition():
    rng = random.Random(1)
    s = "ACGTACGTACAGCTGACCGTGATAATCGGCAGA"
    shuffled = dinuc_shuffle.dinuc_shuffle(s, rng)
    def _dinuc(seq):
        return sorted(seq[i:i+2] for i in range(len(seq) - 1))
    assert _dinuc(shuffled) == _dinuc(s)


def _make_hit(q, s, pident, length, qlen, slen, sstart, send, strand,
              bitscore=100.0, evalue=1e-20):
    return search.HitRow(
        qseqid=q, sseqid=s, pident=pident, length=length, mismatch=0, gapopen=0,
        qstart=1, qend=length, sstart=sstart, send=send, evalue=evalue,
        bitscore=bitscore, qlen=qlen, slen=slen, sstrand=strand,
    )


def test_hit_start0_end0_both_strands():
    h_plus = _make_hit("q1", "ctg", 95.0, 30, 30, 1000, 101, 130, "plus")
    assert search.hit_start0_end0(h_plus) == (100, 130)
    h_minus = _make_hit("q1", "ctg", 95.0, 30, 30, 1000, 130, 101, "minus")
    assert search.hit_start0_end0(h_minus) == (100, 130)


def test_hit_passes_coverage_identity():
    h = _make_hit("q", "s", 92.0, 29, 30, 1000, 1, 29, "plus")
    assert search.hit_passes(h, 0.90, 0.90)   # ident 92, cov 29/30
    assert not search.hit_passes(h, 0.95, 0.90)


def test_pair_finder_rejects_overlap_and_gap_out_of_range():
    repeat_len = 30
    # anchor at 100-130 plus strand; second at 200-230 (gap 70 — too long)
    h1 = _make_hit("q1", "ctg", 100.0, 30, 30, 1000, 101, 130, "plus")
    h2 = _make_hit("q1", "ctg", 100.0, 30, 30, 1000, 201, 230, "plus")
    pairs = pair_mod.find_pairs([h1, h2],
                                 anchor_min_cov=0.9, anchor_min_id=0.9,
                                 second_min_cov=0.6, second_min_id=0.7,
                                 spacer_min=15, spacer_max=60,
                                 overlap_merge_bp=10)
    assert pairs == []
    # second at 170-200: gap = 170 - 130 = 40 bp → valid
    h3 = _make_hit("q1", "ctg", 100.0, 30, 30, 1000, 171, 200, "plus")
    pairs = pair_mod.find_pairs([h1, h3],
                                 anchor_min_cov=0.9, anchor_min_id=0.9,
                                 second_min_cov=0.6, second_min_id=0.7,
                                 spacer_min=15, spacer_max=60,
                                 overlap_merge_bp=10)
    assert len(pairs) == 1
    assert pairs[0].spacer_len == 40
    assert pairs[0].anchor_start0 == 100 and pairs[0].second_end0 == 200


def test_pair_finder_strand_specific():
    """A plus hit and a minus hit should NOT pair even on same contig."""
    h1 = _make_hit("q1", "ctg", 100.0, 30, 30, 1000, 101, 130, "plus")
    h2 = _make_hit("q1", "ctg", 100.0, 30, 30, 1000, 200, 171, "minus")
    pairs = pair_mod.find_pairs([h1, h2],
                                 anchor_min_cov=0.9, anchor_min_id=0.9,
                                 second_min_cov=0.6, second_min_id=0.7,
                                 spacer_min=15, spacer_max=60,
                                 overlap_merge_bp=10)
    assert pairs == []


def test_context_bucket_known_interior_vs_terminal(tmp_path: Path):
    bed = tmp_path / "known.bed"
    bed.write_text("ctg\t100\t200\tA\t500\t+\n")
    arrays_all = tmp_path / "arrays.tsv"
    arrays_all.write_text("coordinate_check_pass\tcategory\tcontig_uid\tstart0\tend0\n")
    hard = context_mod.load_known_array_bed(bed)
    soft = {}
    lengths = {"ctg": 1000}
    # pair interior to known array
    assert context_mod.classify_locus("ctg", 110, 190, hard, soft, lengths) == "known_array_interior"
    # pair abutting right side within 50 bp
    assert context_mod.classify_locus("ctg", 220, 250, hard, soft, lengths) == "suspected_terminal_extension"
    # pair far away
    assert context_mod.classify_locus("ctg", 500, 580, hard, soft, lengths) == "standalone"
    # pair at contig edge
    assert context_mod.classify_locus("ctg", 5, 95, hard, soft, lengths) == "suspected_terminal_extension"
    assert context_mod.classify_locus("ctg", 920, 999, hard, soft, lengths) == "contig_edge"


def _dummy_pair(contig, strand, a0, a1, s0, s1, q):
    return pair_mod.RepeatPair(
        contig=contig, strand=strand,
        anchor_qseqid=q, anchor_pident=100.0, anchor_cov=1.0,
        anchor_start0=a0, anchor_end0=a1,
        second_qseqid=q, second_pident=100.0, second_cov=1.0,
        second_start0=s0, second_end0=s1,
        spacer_len=s0 - a1, spacer_start0=a1, spacer_end0=s0,
        anchor_is_strict=True, second_is_full=True,
    )


def test_dedup_collapses_overlapping_pairs_preserving_queries():
    # Two queries hitting the same locus should collapse to one, keep both query ids.
    pairs = [
        _dummy_pair("ctg", "plus", 100, 130, 170, 200, "qA"),
        _dummy_pair("ctg", "plus", 102, 132, 172, 202, "qB"),
    ]
    loci = dedup_mod.dedup(pairs, merge_bp=60)
    assert len(loci) == 1
    L = loci[0]
    assert {"qA", "qB"} == L.queries
    assert L.n_pairs == 2
