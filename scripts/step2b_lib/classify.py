"""Assign a structural_class to each candidate based on overlay evidence.

Classes (user requested):
  - ordinary_array_fragment      — the candidate is inside, or significantly
                                     overlaps, a known identify or cctyper
                                     array (any category, any score)
  - terminal_extension           — the candidate sits within flank_bp of a
                                     known array end (<50% overlap)
  - candidate_mini_array         — ≥2 repeat units detectable in the window
                                     via identify_units AND outside the
                                     identify/cctyper array footprint
  - repeat_like_ambiguous        — detectable repeat similarity but no second
                                     unit visible in flank; or inside a soft
                                     Possible call
  - edge_unresolved              — not enough flank to tell (either side
                                     short of flank_bp) and no identify
                                     context
"""
from __future__ import annotations

from typing import Any


def classify(candidate_start0: int, candidate_end0: int,
             overlay: dict[str, Any],
             left_flank_bp: int, right_flank_bp: int,
             flank_bp: int, terminal_bp: int = 100) -> dict:
    reasons: list[str] = []
    arrays = overlay["nearby_identify_arrays"]
    cctyper = overlay["nearby_cctyper_crisprs"]
    units_inside = [u for u in overlay["nearby_identify_units"]
                    if u["inside_candidate"]]
    units_in_flank = [u for u in overlay["nearby_identify_units"]
                      if not u["inside_candidate"]]

    cand_len = max(1, candidate_end0 - candidate_start0)

    best_overlap_any = max(
        ([a["overlap_bp"] for a in arrays if a["overlap_bp"] > 0] or [0])
    )
    best_overlap_cct = max(
        ([c["overlap_bp"] for c in cctyper if c["overlap_bp"] > 0] or [0])
    )

    # Ordinary array fragment: ≥50% overlap with any identify / cctyper array.
    if best_overlap_any / cand_len >= 0.50 or best_overlap_cct / cand_len >= 0.50:
        reasons.append(f"overlap_identify={best_overlap_any}bp / "
                        f"cctyper={best_overlap_cct}bp "
                        f"(candidate_len={cand_len})")
        return {"structural_class": "ordinary_array_fragment",
                "reasons": reasons,
                "n_units_inside": len(units_inside),
                "n_units_in_flank": len(units_in_flank)}

    # Terminal extension: abuts a known array within terminal_bp but <50% ovl
    def _abuts(a):
        return (candidate_end0 <= a["end0"] + terminal_bp
                and candidate_start0 >= a["start0"] - terminal_bp
                and (abs(candidate_end0 - a["start0"]) <= terminal_bp
                     or abs(a["end0"] - candidate_start0) <= terminal_bp))
    if any(_abuts(a) for a in arrays) or any(_abuts(a) for a in cctyper):
        reasons.append("abuts_known_array_within_terminal_bp")
        return {"structural_class": "terminal_extension",
                "reasons": reasons,
                "n_units_inside": len(units_inside),
                "n_units_in_flank": len(units_in_flank)}

    # Candidate mini-array: at least 2 units detectable in window and we are
    # outside any known array footprint.
    if len(units_inside) + len(units_in_flank) >= 2 and best_overlap_any == 0:
        reasons.append(f"n_units_in_window="
                        f"{len(units_inside)+len(units_in_flank)}")
        return {"structural_class": "candidate_mini_array",
                "reasons": reasons,
                "n_units_inside": len(units_inside),
                "n_units_in_flank": len(units_in_flank)}

    # Repeat-like ambiguous: inside/overlapping a soft (Possible / Alternative
    # / Possible Discarded / Low score) call, or at least one unit visible
    # somewhere nearby but no complete pair.
    soft_overlap_bp = sum(
        a["overlap_bp"] for a in arrays
        if a["category_raw"] in ("Possible", "Alternative",
                                  "Possible Discarded", "Low score")
    )
    if soft_overlap_bp > 0 or len(units_inside) + len(units_in_flank) >= 1:
        reasons.append(f"soft_overlap_bp={soft_overlap_bp}, "
                        f"units_in_window="
                        f"{len(units_inside)+len(units_in_flank)}")
        return {"structural_class": "repeat_like_ambiguous",
                "reasons": reasons,
                "n_units_inside": len(units_inside),
                "n_units_in_flank": len(units_in_flank)}

    # Edge-unresolved: short flank on either side AND no array context.
    if left_flank_bp < flank_bp or right_flank_bp < flank_bp:
        reasons.append(f"flank_left={left_flank_bp} right={right_flank_bp} "
                        f"required={flank_bp}")
        return {"structural_class": "edge_unresolved",
                "reasons": reasons,
                "n_units_inside": len(units_inside),
                "n_units_in_flank": len(units_in_flank)}

    # Default: standalone with no interpretable context.
    reasons.append("no_overlay_context")
    return {"structural_class": "candidate_mini_array",
            "reasons": reasons,
            "n_units_inside": len(units_inside),
            "n_units_in_flank": len(units_in_flank)}
