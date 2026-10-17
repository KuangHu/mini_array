"""Generate report.md answering step1.md §13 questions."""
from __future__ import annotations

import csv
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def _read_tsv(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _fmt_counter(c: Counter, limit: int = 10) -> str:
    items = sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))[:limit]
    return "\n".join(f"- `{k}` — {v}" for k, v in items) or "- (empty)"


def build_report(project_root: Path, summary: dict, out_path: Path) -> None:
    out_dir = project_root / "results" / "step1"
    arrays = _read_tsv(out_dir / "arrays.all.tsv")
    status = _read_tsv(out_dir / "run_status.tsv")
    links = _read_tsv(out_dir / "array_cas_links.tsv")
    manifest = _read_tsv(project_root / "input" / "genomes.tsv")

    # Pilot provenance
    build_json = project_root / "provenance" / "pilot_manifest.build.json"
    build = json.loads(build_json.read_text()) if build_json.is_file() else {}

    identify_status = Counter(r["identify_status"] for r in status)
    cctyper_status = Counter(r["cctyper_status"] for r in status)
    tier_counts = Counter(a.get("qc_tier", "") for a in arrays)
    cat_counts = Counter(a.get("category", "") for a in arrays)
    cas_counts = Counter(l.get("cas_link_status", "") for l in links)

    # Phylum coverage (from manifest)
    phyla = Counter(m.get("independent_group_id", "NA") for m in manifest)

    # Score threshold comparison per step1.md §13 Q3
    hi_bona = sum(1 for a in arrays
                  if a.get("category") == "bona_fide"
                  and a.get("coordinate_check_pass") == "true"
                  and float(a.get("score") or 0) >= 0.90)
    mid_bona = sum(1 for a in arrays
                   if a.get("category") == "bona_fide"
                   and a.get("coordinate_check_pass") == "true"
                   and float(a.get("score") or 0) >= 0.75)

    md: list[str] = []
    md.append("# Step 1 pilot report")
    md.append("")
    md.append(f"Generated: {datetime.now(timezone.utc).isoformat()}")
    md.append("")

    md.append("## 1. Inputs and run status")
    md.append("")
    md.append(f"- Manifest rows: **{len(manifest)}**  "
              f"(seed {build.get('seed')}, GTDB release {build.get('gtdb_release')})")
    md.append(f"- Phyla sampled: {len(phyla)}")
    md.append("")
    md.append("### Identify")
    md.append(_fmt_counter(identify_status))
    md.append("")
    md.append("### CRISPRCasTyper")
    md.append(_fmt_counter(cctyper_status))
    md.append("")

    md.append("## 2. Candidate category breakdown (post-coordinate-check)")
    md.append("")
    md.append(_fmt_counter(cat_counts))
    md.append("")
    md.append("## 3. QC tier counts")
    md.append("")
    md.append(_fmt_counter(tier_counts))
    md.append(f"- Coordinate check failures (upstream-bug isolate): "
              f"**{summary.get('coordinate_check_fail', '?')}**")
    md.append("")
    md.append("## 4. Score threshold diversity (bona-fide, coord-check pass)")
    md.append("")
    md.append(f"- score ≥ 0.75: **{mid_bona}**")
    md.append(f"- score ≥ 0.90: **{hi_bona}**")
    md.append("")
    md.append("## 5. NR100 + families")
    md.append("")
    md.append(f"- HQ arrays: **{summary.get('hq_arrays')}**")
    md.append(f"- Reserve arrays: **{summary.get('reserve_arrays')}**")
    md.append(f"- NR100 repeat sequences: **{summary.get('nr100_repeats')}**")
    md.append(f"- Operational families: **{summary.get('families')}**")
    md.append("")
    md.append("## 6. Cas link distribution")
    md.append("")
    md.append(_fmt_counter(cas_counts))
    md.append("")
    md.append("## 7. Outputs")
    md.append("")
    for name, p in sorted(summary.get("outputs", {}).items()):
        md.append(f"- `{name}` → `{p}`")
    md.append("")
    md.append("## 8. Caveats")
    md.append("")
    md.append("- No raw reads available; `assembly_validation = not_assessed`.")
    md.append("- Score thresholds (0.90 HQ / 0.75 reserve) are project defaults, "
              "not calibrated to a known-truth benchmark.")
    md.append("- Families are an *operational* grouping (90% norm. edit sim, 90% "
              "length ratio), not evolutionary groups.")
    md.append("- GTDB r220 reps are mostly MAGs (91/100 in this pilot). "
              "cctyper runs with `--prodigal meta` for MAGs/phage/plasmid/unknown.")
    md.append("")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(md))
