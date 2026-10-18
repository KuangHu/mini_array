"""Phase K2: structural_class based on actual repeat unit evidence,
not just array-level overlap.

Rules (per the user instruction):
  - ordinary_array_fragment: candidate overlaps ≥2 identify repeat UNITS
    within its footprint, OR overlaps a Bona-fide array that has units
    reaching into the footprint.
  - terminal_extension: abuts a known array (±100 bp) and has unit evidence
    just outside the candidate.
  - candidate_mini_array: has repeat pair evidence (from Step 2) and either
    no overlap with units/arrays, or only overlaps with low-quality calls
    that lack units in the footprint.
  - repeat_like_ambiguous: overlaps a Possible/Low score band but the
    unit arrangement is partial (<2 units) or ambiguous.
  - edge_excluded: either flank below `flank_bp`, regardless of array context.
"""
from __future__ import annotations


def classify(candidate_start0: int, candidate_end0: int,
              overlay: dict,
              left_flank_bp: int, right_flank_bp: int,
              flank_bp: int, terminal_bp: int = 100) -> dict:
    reasons = []
    if left_flank_bp < flank_bp or right_flank_bp < flank_bp:
        return {"structural_class": "edge_excluded",
                "reasons": [f"edge_short_flank left={left_flank_bp} right={right_flank_bp} required={flank_bp}"],
                "n_units_inside": 0, "n_units_in_flank": 0,
                "evidence_arrays_with_units": []}
    units_inside = [u for u in overlay["nearby_identify_units"]
                    if u["inside_candidate"] and u["unit_kind"] == "repeat"]
    units_in_flank = [u for u in overlay["nearby_identify_units"]
                      if not u["inside_candidate"] and u["unit_kind"] == "repeat"]
    arrays = overlay["nearby_identify_arrays"]
    cctyper = overlay["nearby_cctyper_crisprs"]

    # Which arrays actually have units reaching into the candidate footprint?
    arrays_with_units_in_candidate: list[str] = []
    for a in arrays:
        aid = a["array_id"]
        if a["category_raw"] == "Bona-fide" and any(u["array_id"] == aid for u in units_inside):
            arrays_with_units_in_candidate.append(aid)

    # cctyper carries its own repeat annotation; treat any overlapping cctyper
    # CRISPR with ≥2 repeats as unit evidence equivalent.
    cctyper_with_many_repeats = [
        c for c in cctyper
        if c["overlap_bp"] > 0 and (c.get("n_repeats") or "0").isdigit()
        and int(c.get("n_repeats") or "0") >= 2
    ]

    if len(units_inside) >= 2 or arrays_with_units_in_candidate or cctyper_with_many_repeats:
        reasons.append(
            f"units_inside={len(units_inside)} "
            f"identify_arrays_with_units_inside={arrays_with_units_in_candidate} "
            f"cctyper_units_inside={len(cctyper_with_many_repeats)}"
        )
        return {
            "structural_class": "ordinary_array_fragment",
            "reasons": reasons,
            "n_units_inside": len(units_inside),
            "n_units_in_flank": len(units_in_flank),
            "evidence_arrays_with_units": arrays_with_units_in_candidate,
        }

    # Terminal extension: has units in flank + abuts a known array end.
    def _abuts(a):
        return (abs(a["start0"] - candidate_end0) <= terminal_bp
                or abs(a["end0"] - candidate_start0) <= terminal_bp)
    if (units_in_flank and (any(_abuts(a) for a in arrays)
                             or any(_abuts(c) for c in cctyper))):
        reasons.append(f"abuts_known_array units_in_flank={len(units_in_flank)}")
        return {
            "structural_class": "terminal_extension",
            "reasons": reasons,
            "n_units_inside": len(units_inside),
            "n_units_in_flank": len(units_in_flank),
            "evidence_arrays_with_units": [],
        }

    # Repeat-like ambiguous: overlaps a non-bona-fide call AND has 1 unit
    # (either inside or in flank) — partial evidence but not an array.
    non_bf_arrays = [a for a in arrays
                     if a["category_raw"] in ("Possible", "Alternative",
                                               "Possible Discarded", "Low score")
                     and a["overlap_bp"] > 0]
    if non_bf_arrays and (len(units_inside) + len(units_in_flank)) >= 1:
        reasons.append(
            f"non_bona_fide_overlap={[a['array_id'] for a in non_bf_arrays]} "
            f"units_in_window={len(units_inside)+len(units_in_flank)}"
        )
        return {
            "structural_class": "repeat_like_ambiguous",
            "reasons": reasons,
            "n_units_inside": len(units_inside),
            "n_units_in_flank": len(units_in_flank),
            "evidence_arrays_with_units": [],
        }

    reasons.append("no_array_or_unit_overlay; has_pair_evidence_from_step2")
    return {
        "structural_class": "candidate_mini_array",
        "reasons": reasons,
        "n_units_inside": len(units_inside),
        "n_units_in_flank": len(units_in_flank),
        "evidence_arrays_with_units": [],
    }
