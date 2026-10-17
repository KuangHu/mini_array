"""Unit tests for parse / qc / families / cas modules (step1.md §12 cases)."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from step1_lib import cas as cas_mod  # noqa: E402
from step1_lib import families as fam_mod  # noqa: E402
from step1_lib import parse as parse_mod  # noqa: E402
from step1_lib import qc as qc_mod  # noqa: E402


def test_revcomp_and_canonical():
    assert parse_mod.revcomp("ACGT") == "ACGT"
    assert parse_mod.revcomp("ACGTN") == "NACGT"
    # canonical should pick the lexicographically smaller of (seq, RC(seq))
    assert parse_mod.canonical("CCCCGGGG") == "CCCCGGGG"
    assert parse_mod.canonical("GGGGCCCC") == "GGGGCCCC"  # self-RC
    assert parse_mod.canonical("GTTTAAAC") == "GTTTAAAC"  # self-RC


def test_repeat_id_rc_equivalent_collapses():
    r1 = "ACGTACGTACGT"
    r2 = parse_mod.revcomp(r1)
    assert parse_mod.repeat_id(r1) == parse_mod.repeat_id(r2)


def _make_cand(starts, repeats, spacers, consensus, strand, score=0.95,
               array_index=1):
    return {
        "list_repeats": repeats,
        "list_repeat_starts": starts,
        "list_spacers": spacers,
        "consensus_repeat": consensus,
        "strand": strand,
        "start": starts[0] + 1,
        "end": starts[-1] + len(repeats[-1]),
        "length": starts[-1] + len(repeats[-1]) - starts[0],
        "array_index": array_index,
        "certainty_score": score,
    }


def test_parse_forward_array_matches_reference():
    repeat = "GTTTCAATCCCCTCCTGCGGGGAACGGCGAAC"
    sp = ["AAAAAAAA", "GGGGGGGG"]
    seq = "N" * 50 + repeat + sp[0] + repeat + sp[1] + repeat + "N" * 50
    cand = _make_cand(
        starts=[50, 50 + len(repeat) + len(sp[0]),
                50 + 2 * len(repeat) + len(sp[0]) + len(sp[1])],
        repeats=[repeat, repeat, repeat],
        spacers=sp,
        consensus=repeat,
        strand="Forward (Orientation was not computed)",
    )
    arr = parse_mod.parse_candidate(
        "Bona-fide", "1", cand, "ctg1", seq, "G1",
    )
    assert arr.coordinate_check_pass
    assert arr.qc_reasons == []
    assert arr.n_repeats == 3
    assert arr.n_spacers == 2
    assert arr.consensus_forward == repeat
    assert arr.start0 == 50
    assert arr.end0 == 50 + 3 * len(repeat) + len(sp[0]) + len(sp[1])


def test_parse_reversed_array_matches_reference():
    forward_repeat = "GTTTCAATCCACGCCCCTGCGCGAGGGGCGAC"
    sp_fwd = ["TTACAAAAAAAATTAGCAGAATCAATGATAAATTT",
              "CAGGAGATGAGGAGCGCAAAGGAACCGTTACTCTT",
              "GAGGAAGGACTAATAATTGAATACAAGAGATTTA"]
    seq = "N" * 59 + forward_repeat
    for s in sp_fwd:
        seq += s + forward_repeat
    seq += "N" * 100

    # In reversed display, repeats are RC(forward) and starts descend in list order.
    display_repeat = parse_mod.revcomp(forward_repeat)
    display_spacers = [parse_mod.revcomp(s) for s in reversed(sp_fwd)]
    offsets = []
    cursor = 59
    for i in range(4):
        offsets.append(cursor)
        cursor += len(forward_repeat)
        if i < 3:
            cursor += len(sp_fwd[i])
    display_starts = list(reversed(offsets))

    cand = _make_cand(
        starts=display_starts,
        repeats=[display_repeat] * 4,
        spacers=display_spacers,
        consensus=display_repeat,
        strand="Reversed",
    )
    arr = parse_mod.parse_candidate(
        "Bona-fide", "1", cand, "ctg_rev", seq, "G1",
    )
    assert arr.coordinate_check_pass, arr.qc_reasons
    assert arr.consensus_forward == forward_repeat
    assert arr.n_repeats == 4
    assert arr.n_spacers == 3
    assert arr.start0 == 59
    assert arr.end0 == cursor


def test_parse_rejects_off_by_one_without_fixing():
    repeat = "GTTTCAATCCCCTCCTGCGGGGAACGGCGAAC"
    sp = "AAAAAAAA"
    seq = "N" * 50 + repeat + sp + repeat + "N" * 50
    # Inject a 1-bp error on the second start.
    cand = _make_cand(
        starts=[50, 50 + len(repeat) + len(sp) + 1],
        repeats=[repeat, repeat],
        spacers=[sp],
        consensus=repeat,
        strand="Forward (Orientation was not computed)",
    )
    arr = parse_mod.parse_candidate(
        "Bona-fide", "1", cand, "ctg_bad", seq, "G1",
    )
    assert not arr.coordinate_check_pass
    assert any("does not match reference" in r for r in arr.qc_reasons)


# --- QC tests ---

class _StubRepeat:
    def __init__(self, s0, e0, seq):
        self.start0 = s0
        self.end0 = e0
        self.sequence_forward = seq
        self.unit_index = 0


class _StubArray:
    def __init__(self, *, category, score, n_rep, n_unique_spacer, consensus,
                 start0, end0, repeats, coord_ok=True, qc_reasons=None):
        self.category_raw = category
        self.score = score
        self.n_repeats = n_rep
        self.n_spacers = max(0, n_rep - 1)
        self.n_unique_spacers = n_unique_spacer
        self.consensus_forward = consensus
        self.start0 = start0
        self.end0 = end0
        self.repeats = repeats
        self.coordinate_check_pass = coord_ok
        self.qc_reasons = qc_reasons or []


def _hq_passing_array():
    cons = "ACGTACGTACGTACGTACGTACGTACGTACGT"
    repeats = [_StubRepeat(100 + i * 60, 100 + i * 60 + 32, cons) for i in range(3)]
    return _StubArray(category="Bona-fide", score=0.95, n_rep=3,
                      n_unique_spacer=2, consensus=cons,
                      start0=100, end0=100 + 3 * 32 + 2 * 28,
                      repeats=repeats)


def test_qc_classify_hq_pass():
    a = _hq_passing_array()
    contig = "A" * 500  # clean flanks
    r = qc_mod.classify(a, contig, qc_mod.default_config())
    assert r["qc_tier"] == "hq", r
    assert r["n_full_repeats"] == 3
    assert r["ambiguous_bases_in_context"] == 0


def test_qc_classify_demotes_low_score():
    a = _hq_passing_array()
    a.score = 0.80
    r = qc_mod.classify(a, "A" * 1000, qc_mod.default_config())
    assert r["qc_tier"] == "reserve"
    assert any("score<" in reason for reason in r["qc_reasons"])


def test_qc_classify_context_only_for_low_score():
    a = _hq_passing_array()
    a.category_raw = "Low score"
    a.score = 0.50
    r = qc_mod.classify(a, "A" * 1000, qc_mod.default_config())
    assert r["qc_tier"] == "context_only"


def test_qc_possible_coord_valid_stays_in_reserve():
    """Phase E: a coord-valid primary Possible with sub-HQ score must stay in
    reserve, not drop to context_only (score alone is insufficient to downgrade
    a Possible — only Alternative/Possible Discarded/Low score go to
    context_only)."""
    a = _hq_passing_array()
    a.category_raw = "Possible"
    a.score = 0.60
    r = qc_mod.classify(a, "A" * 1000, qc_mod.default_config())
    assert r["qc_tier"] == "reserve"
    # structural fields are still populated for audit
    assert r["n_full_repeats"] == 3


def test_qc_alternative_and_possible_discarded_go_context_only():
    for cat in ("Alternative", "Possible Discarded"):
        a = _hq_passing_array()
        a.category_raw = cat
        a.score = 0.80
        r = qc_mod.classify(a, "A" * 1000, qc_mod.default_config())
        assert r["qc_tier"] == "context_only", (cat, r)


def test_qc_classify_invalid_when_coord_check_fails():
    a = _hq_passing_array()
    a.coordinate_check_pass = False
    r = qc_mod.classify(a, "A" * 1000, qc_mod.default_config())
    assert r["qc_tier"] == "invalid"


def test_qc_classify_terminal_degenerate_allowed_but_flagged():
    cons = "ACGTACGTACGTACGTACGTACGTACGTACGT"
    bad_repeat = "N" * 32  # too divergent, will fail similarity threshold
    repeats = [
        _StubRepeat(100, 132, cons),
        _StubRepeat(160, 192, cons),
        _StubRepeat(220, 252, cons),
        _StubRepeat(280, 312, bad_repeat),
    ]
    a = _StubArray(category="Bona-fide", score=0.95, n_rep=4,
                   n_unique_spacer=3, consensus=cons,
                   start0=100, end0=312, repeats=repeats)
    r = qc_mod.classify(a, "A" * 800, qc_mod.default_config())
    assert r["n_full_repeats"] == 3
    assert r["terminal_repeat_degenerate"]
    assert r["qc_tier"] == "hq"  # 3 full support repeats satisfy the threshold


# --- Families tests ---

class _CatArray:
    """Minimal object compatible with families.build_nr100."""
    def __init__(self, consensus_forward, array_id="arr", genome_id="G",
                 contig_uid="C", orientation="forward"):
        self.consensus_forward = consensus_forward
        self.array_id = array_id
        self.genome_id = genome_id
        self.contig_uid = contig_uid
        self.orientation_prediction = orientation


def test_nr100_collapses_rc_pairs():
    seq = "ACGTACGTACGTACGTACGTACGTAC"
    rc = parse_mod.revcomp(seq)
    arrays = [_CatArray(seq, "a1"), _CatArray(rc, "a2")]
    catalog = fam_mod.build_nr100(arrays)
    assert len(catalog) == 1
    rid = next(iter(catalog))
    assert {s["array_id"] for s in catalog[rid].sources} == {"a1", "a2"}


def test_family_cluster_no_transitive_merge():
    # A-B pass, B-C pass, A-C fail → C must get its own family.
    base = "ACGTACGTACGTACGTACGTACGTACGT"      # 28 bp
    a = base
    b = base[:14] + "T" + base[15:]            # 1 mismatch vs A  (sim ~0.964)
    c = base[:14] + "T" + "GCGCGCGCGCGCGC"     # ~14 mismatches vs A, ~13 vs B
    catalog = fam_mod.build_nr100([
        _CatArray(a, "a1"), _CatArray(b, "a2"), _CatArray(c, "a3"),
    ])
    result = fam_mod.cluster_families(catalog, 0.90, 0.90)
    fam_of = {rid: info["family_id"]
              for rid, info in result["assignments"].items()}
    rid_a = parse_mod.repeat_id(a)
    rid_b = parse_mod.repeat_id(b)
    rid_c = parse_mod.repeat_id(c)
    assert fam_of[rid_a] == fam_of[rid_b]
    assert fam_of[rid_c] != fam_of[rid_a]


# --- Cas tests ---

def test_cas_gap_bp_within_window():
    # Array [100, 200), operon 1-based [300, 500] → 0-based [299, 500). Gap
    # is the number of positions BETWEEN the two intervals = 299 - 200 = 99.
    assert cas_mod._gap_bp(100, 200, 300, 500) == 99
    assert cas_mod._gap_bp(100, 200, 100, 150) == 0   # overlap
    # Array [300, 400), operon 1-based [100, 200] → 0-based [99, 200). Gap
    # = 300 - 200 = 100.
    assert cas_mod._gap_bp(300, 400, 100, 200) == 100


def test_cas_link_no_cross_contig():
    op = cas_mod.Operon(contig_uid="ctg_other", operon_id="op1",
                        start1=1, end1=10000, prediction="I-C",
                        best_type="I-C", best_score=10.0, confidence="high",
                        genes="", source_table="cas_operons.tab")

    class A:
        array_id = "arr1"
        genome_id = "G1"
        contig_uid = "ctg_this"
        start0 = 100
        end0 = 200
    links = cas_mod.link_arrays_to_operons([A()], {"G1": [op]}, 10000)
    assert len(links) == 1
    assert links[0]["cas_link_status"] == "no_proximal_operon"
    assert links[0]["operon_id"] == ""
