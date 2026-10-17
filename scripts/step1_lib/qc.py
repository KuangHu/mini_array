"""Apply seed_qc tiers (step1.md §5).

Takes ParsedArray objects + loaded staged FASTA sequences and classifies each
array into one of:

  hq               — primary seed library candidate (meets all defaults)
  reserve          — bona-fide 0.75-0.90, or structurally weak but validated
  context_only     — Alternative / Possible / Possible Discarded / Low score
                     calls kept only for context (soft mask)
  invalid          — coordinate verification failed; stays in arrays.all.tsv
                     with coordinate_check_pass=false, excluded from HQ/reserve

Flank / ambiguous-base checks need the staged FASTA, which we load via a
callable so the caller controls IO.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

try:
    import edlib
except ImportError:  # pragma: no cover
    edlib = None


ALLOWED_BASES = set("ACGT")  # for flank-cleanliness check


@dataclass
class QcConfig:
    min_score: float
    min_full_repeats: int
    min_unique_spacers: int
    context_bp: int
    max_ambiguous_bases_in_context: int
    min_copy_normalized_edit_similarity: float
    min_copy_length_ratio: float


def default_config() -> QcConfig:
    return QcConfig(
        min_score=0.90,
        min_full_repeats=3,
        min_unique_spacers=2,
        context_bp=100,
        max_ambiguous_bases_in_context=0,
        min_copy_normalized_edit_similarity=0.90,
        min_copy_length_ratio=0.90,
    )


def _normalized_edit_similarity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    if edlib is not None:
        d = edlib.align(a, b, task="distance", mode="NW")["editDistance"]
    else:  # pragma: no cover — fallback
        import difflib
        d = int(max(len(a), len(b)) * (1 - difflib.SequenceMatcher(None, a, b).ratio()))
    return 1 - d / max(len(a), len(b))


def _length_ratio(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    return min(len(a), len(b)) / max(len(a), len(b))


def count_full_support_repeats(array, cfg: QcConfig) -> int:
    """Count repeats whose normalized edit similarity and length ratio vs the
    forward consensus meet both thresholds (step1.md §5 definition)."""
    cons = array.consensus_forward
    if not cons:
        return 0
    ok = 0
    for r in array.repeats:
        if (_normalized_edit_similarity(r.sequence_forward, cons)
                >= cfg.min_copy_normalized_edit_similarity
                and _length_ratio(r.sequence_forward, cons)
                >= cfg.min_copy_length_ratio):
            ok += 1
    return ok


def _flank_stats(contig_seq: str, start0: int, end0: int, context_bp: int
                ) -> tuple[int, int, int, bool]:
    """Return (left_bp, right_bp, ambiguous_in_context, bounded_both_sides)."""
    left_bp = max(0, min(start0, context_bp))
    right_bp = max(0, min(len(contig_seq) - end0, context_bp))
    left_slice = contig_seq[max(0, start0 - context_bp):start0]
    right_slice = contig_seq[end0:end0 + context_bp]
    array_slice = contig_seq[start0:end0]
    full_context = left_slice + array_slice + right_slice
    ambiguous = sum(1 for ch in full_context if ch not in ALLOWED_BASES)
    bounded = left_bp >= context_bp and right_bp >= context_bp
    return left_bp, right_bp, ambiguous, bounded


def classify(array, contig_seq: str | None, cfg: QcConfig) -> dict:
    """Return QC info for one ParsedArray; attach reasons list for audit.

    Tier assignment per step1.md §5 with the Phase E correction:
      - coord-check fail → `invalid`
      - category == Bona-fide, all checks pass → `hq`
      - category == Bona-fide, one or more checks fail → `reserve`
      - category == Possible, coord-valid → `reserve` (kept as independent
        soft-support seeds; a sub-HQ score alone is NOT sufficient reason
        to drop them to context_only)
      - Alternative / Possible Discarded / Low score → `context_only`
    """
    reasons: list[str] = []
    if not array.coordinate_check_pass:
        return {
            "qc_tier": "invalid",
            "qc_reasons": ["coordinate_check_failed"] + array.qc_reasons,
            "n_full_repeats": 0,
            "left_flank_bp": 0,
            "right_flank_bp": 0,
            "ambiguous_bases_in_context": -1,
            "terminal_repeat_degenerate": False,
        }

    cat = array.category_raw
    if cat == "Bona-fide":
        tier = "hq"
    elif cat == "Possible":
        # Phase E: coord-valid Possibles stay in reserve regardless of score;
        # subsequent structural/flank/ambiguous-base failures stay in reserve
        # too (they only disqualify Bona-fide candidates from HQ, never
        # downgrade Possibles further).
        tier = "reserve"
    else:  # Alternative, Possible Discarded, Low score
        tier = "context_only"
        reasons.append(f"category={cat}")

    # Structural / contextual checks — always computed for the delivery
    # columns, even for context_only, so the audit tables stay complete.
    n_full = count_full_support_repeats(array, cfg)
    terminal_degenerate = (array.n_repeats >= cfg.min_full_repeats and
                           n_full < array.n_repeats)

    if cat == "Bona-fide":
        if array.score < cfg.min_score:
            reasons.append(f"score<{cfg.min_score:.2f}")
            tier = "reserve"
        if n_full < cfg.min_full_repeats:
            reasons.append(f"n_full_repeats<{cfg.min_full_repeats} (got {n_full})")
            tier = "reserve"
        if array.n_unique_spacers < cfg.min_unique_spacers:
            reasons.append(f"n_unique_spacers<{cfg.min_unique_spacers}")
            tier = "reserve"

    # Flank measurement uses the real contig length. The _flank_stats helper
    # already pulls len(contig_seq) directly; we just don't cap the reported
    # values so the delivery columns reflect what's actually available.
    if contig_seq is not None:
        left_bp = max(0, array.start0)
        right_bp = max(0, len(contig_seq) - array.end0)
        full_context = contig_seq[max(0, array.start0 - cfg.context_bp):array.end0 + cfg.context_bp]
        amb = sum(1 for ch in full_context if ch not in ALLOWED_BASES)
        bounded = (left_bp >= cfg.context_bp and right_bp >= cfg.context_bp)
        if cat == "Bona-fide":
            if not bounded:
                reasons.append(
                    f"flank<{cfg.context_bp}bp (left={left_bp}, right={right_bp})"
                )
                tier = "reserve"
            if amb > cfg.max_ambiguous_bases_in_context:
                reasons.append(f"ambiguous_bases_in_context={amb}")
                tier = "reserve"
    else:
        left_bp = right_bp = 0
        amb = -1
        if cat == "Bona-fide":
            reasons.append("contig_sequence_unavailable")
            tier = "reserve"

    if cat == "Bona-fide" and (not array.consensus_forward
                               or set(array.consensus_forward.upper()) - ALLOWED_BASES):
        reasons.append("consensus_contains_ambiguous_or_empty")
        tier = "reserve"

    return {
        "qc_tier": tier,
        "qc_reasons": reasons,
        "n_full_repeats": n_full,
        "left_flank_bp": left_bp,
        "right_flank_bp": right_bp,
        "ambiguous_bases_in_context": amb,
        "terminal_repeat_degenerate": terminal_degenerate,
    }
