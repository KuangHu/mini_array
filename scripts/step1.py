#!/usr/bin/env python3
"""Step 1 driver: CRISPRidentify array recovery + repeat seed library.

CLI contract (per step1.md section 11):

  python scripts/step1.py doctor   --config config/step1.yaml
  python scripts/step1.py run      --config config/step1.yaml --limit 100 --resume
  python scripts/step1.py validate --config config/step1.yaml

This is a stub skeleton. Component functions (prepare_inputs, run_identify,
parse_identify, apply_seed_qc, build_repeat_catalog, cluster_repeat_families,
run_cctyper, link_array_cas, export_and_report) will be filled in incrementally
once the identify env is built and a small test fixture is in place.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step1_lib import cas as cas_mod  # noqa: E402
from step1_lib import families as families_mod  # noqa: E402
from step1_lib import parse as parse_mod  # noqa: E402
from step1_lib import qc as qc_mod  # noqa: E402
from step1_lib import report as report_mod  # noqa: E402
from step1_lib import stage as stage_mod  # noqa: E402
from step1_lib import status as status_mod  # noqa: E402
from step1_lib import submit as submit_mod  # noqa: E402

try:
    import yaml
except ImportError:
    sys.stderr.write(
        "PyYAML not installed in the current Python. Install into the driver env "
        "or run inside envs/driver.\n"
    )
    raise


def load_config(path: Path) -> dict:
    with path.open() as fh:
        cfg = yaml.safe_load(fh)
    for key in ("project_root", "code_root"):
        if key not in cfg:
            raise SystemExit(f"config missing required key: {key}")
    cfg["project_root"] = Path(cfg["project_root"]).resolve()
    cfg["code_root"] = Path(cfg["code_root"]).resolve()
    return cfg


def cmd_doctor(cfg: dict) -> int:
    """Report environment, software commits, and whether required paths exist."""
    p = cfg["project_root"]
    checks = []

    def check(label: str, ok: bool, detail: str = "") -> None:
        checks.append((label, ok, detail))

    check("project_root exists", p.is_dir(), str(p))
    check("code_root exists", cfg["code_root"].is_dir(), str(cfg["code_root"]))
    for sub in ("software", "envs", "input/staged", "work", "raw/identify",
                "raw/cctyper", "results/step1", "provenance", "logs"):
        check(f"{sub}/", (p / sub).is_dir(), str(p / sub))

    ident_py = p / "software/CRISPRidentify/CRISPRidentify.py"
    check("CRISPRidentify.py", ident_py.is_file(), str(ident_py))
    cct_dir = p / "software/CRISPRCasTyper"
    check("CRISPRCasTyper repo", cct_dir.is_dir(), str(cct_dir))

    for env in ("identify", "driver"):
        env_prefix = p / "envs" / env
        check(f"envs/{env}", (env_prefix / "conda-meta").is_dir(), str(env_prefix))

    for name in ("CRISPRidentify.commit", "CRISPRCasTyper.commit"):
        f = p / "provenance" / name
        detail = f.read_text().strip() if f.is_file() else str(f)
        check(f"provenance/{name}", f.is_file(), detail)

    manifest = p / cfg["manifest"]
    check("manifest", manifest.is_file(), str(manifest))

    width = max(len(label) for label, _, _ in checks)
    fail = 0
    for label, ok, detail in checks:
        mark = "ok " if ok else "MISS"
        print(f"[{mark}] {label.ljust(width)}  {detail}")
        if not ok:
            fail += 1
    print(f"\n{len(checks) - fail}/{len(checks)} checks passed.")
    return 0 if fail == 0 else 1


def cmd_prepare(cfg: dict, limit: int | None) -> int:
    p = cfg["project_root"]
    manifest_path = p / cfg["manifest"]
    rows = stage_mod.read_manifest(manifest_path)
    summary = stage_mod.stage_manifest(
        manifest_rows=rows,
        staged_dir=p / "input" / "staged",
        contigs_tsv=p / "results" / "step1" / "contigs.tsv",
        limit=limit,
    )
    print(json.dumps(summary, indent=2))
    return 0 if summary["failed_count"] == 0 else 1


def cmd_submit(cfg: dict, skip_existing: bool) -> int:
    p = cfg["project_root"]
    rows = submit_mod.read_manifest(p / cfg["manifest"])
    out_path = cfg["code_root"] / "scripts/sbatch/submit_pilot.sh"
    summary = submit_mod.build_submit_script(
        manifest_rows=rows,
        staged_dir=p / "input" / "staged",
        raw_identify_dir=p / "raw" / "identify",
        raw_cctyper_dir=p / "raw" / "cctyper",
        code_root=cfg["code_root"],
        out_path=out_path,
        skip_existing=skip_existing,
    )
    print(json.dumps(summary, indent=2))
    return 0


def cmd_parse(cfg: dict) -> int:
    p = cfg["project_root"]
    rows = submit_mod.read_manifest(p / cfg["manifest"])
    summary = parse_mod.parse_all(
        project_root=p,
        manifest_rows=rows,
        arrays_tsv=p / "results" / "step1" / "arrays.all.tsv",
        units_tsv=p / "results" / "step1" / "array_units.tsv",
    )
    print(json.dumps(summary, indent=2))
    return 0 if summary["coordinate_check_fail"] == 0 else 1


def cmd_export(cfg: dict, min_score_override: float | None = None,
               out_subdir: str = "step1",
               identify_source: str = "raw/identify") -> int:
    """Full export: parse → QC → Cas links → NR100 + families → tables.

    `identify_source` picks which identify output tree to parse — the
    original "raw/identify" or the patched "raw_patched/identify". The cctyper
    outputs under raw/cctyper are shared (inputs + cctyper version unchanged).
    """
    from collections import Counter

    p = cfg["project_root"]
    out_dir = p / "results" / out_subdir
    manifest_rows = submit_mod.read_manifest(p / cfg["manifest"])

    qc_cfg_dict = cfg.get("seed_qc", {}) or {}
    fam_cfg_dict = cfg.get("families", {}) or {}
    cas_cfg_dict = cfg.get("cas_annotation", {}) or {}
    effective_min_score = (min_score_override
                           if min_score_override is not None
                           else qc_cfg_dict.get("min_score", 0.9))
    qc_cfg = qc_mod.QcConfig(
        min_score=effective_min_score,
        min_full_repeats=qc_cfg_dict.get("min_full_repeats", 3),
        min_unique_spacers=qc_cfg_dict.get("min_unique_spacers", 2),
        context_bp=qc_cfg_dict.get("context_bp", 100),
        max_ambiguous_bases_in_context=qc_cfg_dict.get(
            "max_ambiguous_bases_in_context", 0),
        min_copy_normalized_edit_similarity=qc_cfg_dict.get(
            "min_copy_normalized_edit_similarity", 0.9),
        min_copy_length_ratio=qc_cfg_dict.get("min_copy_length_ratio", 0.9),
    )

    all_arrays: list = []
    arrays_rows: list[dict] = []
    missing_identify: list[str] = []

    identify_root = p / identify_source
    for row in manifest_rows:
        gid = row["genome_id"]
        guid = stage_mod.sanitize(gid)
        raw_identify = identify_root / guid
        staged = p / "input" / "staged" / f"{guid}.fna"
        if not (raw_identify / "Complete_json_summary.json").is_file():
            missing_identify.append(gid)
            continue
        arrays = parse_mod.parse_genome(gid, guid, raw_identify, staged)
        contig_seqs = parse_mod.read_staged_fasta(staged) if arrays else {}
        for a in arrays:
            qc_result = qc_mod.classify(a, contig_seqs.get(a.contig_uid), qc_cfg)
            # Attach QC info directly to the ParsedArray for downstream use.
            a.qc_tier = qc_result["qc_tier"]
            a.n_full_repeats = qc_result["n_full_repeats"]
            a.left_flank_bp = qc_result["left_flank_bp"]
            a.right_flank_bp = qc_result["right_flank_bp"]
            a.ambiguous_bases_in_context = qc_result["ambiguous_bases_in_context"]
            a.terminal_repeat_degenerate = qc_result["terminal_repeat_degenerate"]
            # Merge QC reasons with any pre-existing coord-check reasons.
            a.qc_reasons = list(dict.fromkeys(qc_result["qc_reasons"] + a.qc_reasons))
            all_arrays.append(a)

    # arrays.all.tsv — one row per candidate across all genomes (full schema
    # with qc_tier etc.). arrays.hq.tsv / arrays.reserve.tsv are filtered views.
    parse_mod.write_arrays_tsv_full(all_arrays, out_dir / "arrays.all.tsv")
    parse_mod.write_units_tsv(all_arrays, out_dir / "array_units.tsv")
    parse_mod.write_arrays_tsv_full(
        [a for a in all_arrays if getattr(a, "qc_tier", "") == "hq"],
        out_dir / "arrays.hq.tsv",
    )
    parse_mod.write_arrays_tsv_full(
        [a for a in all_arrays if getattr(a, "qc_tier", "") == "reserve"],
        out_dir / "arrays.reserve.tsv",
    )
    n_bed = parse_mod.write_known_array_regions_bed(
        all_arrays, out_dir / "known_array_regions.bed"
    )
    # Per-genome status table
    status_rows = status_mod.build_rows(p, manifest_rows)
    status_mod.write_tsv(status_rows, out_dir / "run_status.tsv")

    # Cas linking
    operons_by_genome, missing_cctyper = cas_mod.load_all_operons(p, manifest_rows)
    max_gap = cas_cfg_dict.get("max_local_gap_bp", 10000)
    links = cas_mod.link_arrays_to_operons(all_arrays, operons_by_genome, max_gap)
    cas_mod.write_cas_links_tsv(links, out_dir / "array_cas_links.tsv")

    # NR100 + families (only HQ arrays feed the search library). To keep family
    # IDs stable across thresholds, the family clustering is run on the UNION
    # of HQ repeats detectable at the broadest threshold we support (0.75). At
    # a stricter threshold (e.g. 0.90) we just filter the catalog down to the
    # HQ-at-this-threshold subset before writing — families remain the same
    # 0.75-induced partition.
    hq = [a for a in all_arrays if a.qc_tier == "hq"]
    reserve = [a for a in all_arrays if a.qc_tier == "reserve"]

    # Build a 0.75 reference catalog from the same arrays set using a 0.75 cfg.
    # Score alone determines HQ vs reserve for a bona-fide candidate at this
    # stage, so we recompute the HQ membership at 0.75 without touching other
    # QC rules.
    qc_cfg_075 = qc_mod.QcConfig(
        min_score=0.75,
        min_full_repeats=qc_cfg.min_full_repeats,
        min_unique_spacers=qc_cfg.min_unique_spacers,
        context_bp=qc_cfg.context_bp,
        max_ambiguous_bases_in_context=qc_cfg.max_ambiguous_bases_in_context,
        min_copy_normalized_edit_similarity=qc_cfg.min_copy_normalized_edit_similarity,
        min_copy_length_ratio=qc_cfg.min_copy_length_ratio,
    )
    hq_075: list = []
    for a in all_arrays:
        if not a.coordinate_check_pass:
            continue
        if a.category_raw != "Bona-fide":
            continue
        # Reuse structural fields already attached during this export.
        contig_seqs_cache = parse_mod.read_staged_fasta(
            p / "input" / "staged" / f"{stage_mod.sanitize(a.genome_id)}.fna"
        ) if False else None  # cached ones are not available at this point
        # Simpler: recompute full tier with the 0.75 config
        qc075 = qc_mod.classify(
            a, parse_mod.read_staged_fasta(
                p / "input" / "staged" / f"{stage_mod.sanitize(a.genome_id)}.fna"
            ).get(a.contig_uid),
            qc_cfg_075,
        )
        if qc075["qc_tier"] == "hq":
            hq_075.append(a)

    catalog_075 = families_mod.build_nr100(hq_075)
    cluster_075 = families_mod.cluster_families(
        catalog_075,
        min_similarity=fam_cfg_dict.get("min_normalized_edit_similarity", 0.9),
        min_length_ratio=fam_cfg_dict.get("min_length_ratio", 0.9),
    )

    # Filter down to the current threshold's HQ subset, re-using the same
    # family IDs (repeat_id → family_id map is injective at both thresholds
    # because the 0.75 catalog is a superset of the 0.90 catalog).
    hq_repeat_ids_here = {parse_mod.repeat_id(a.consensus_forward) for a in hq}
    catalog_here = {rid: entry for rid, entry in catalog_075.items()
                    if rid in hq_repeat_ids_here}
    # keep only families touched by at least one HQ-here repeat
    cluster_here = {
        "families": {},
        "assignments": {rid: info for rid, info in cluster_075["assignments"].items()
                        if rid in hq_repeat_ids_here},
        "centroid_order": [],
    }
    touched = {info["family_id"] for info in cluster_here["assignments"].values()}
    for fid in cluster_075["centroid_order"]:
        if fid not in touched:
            continue
        f = dict(cluster_075["families"][fid])
        f["members"] = [m for m in f["members"] if m in hq_repeat_ids_here]
        cluster_here["families"][fid] = f
        cluster_here["centroid_order"].append(fid)

    fam_summary = families_mod.write_all(
        catalog_here, cluster_here,
        out_dir=out_dir,
        fasta_path=out_dir / "repeats.hq.nr100.fna",
    )

    tier_counts = Counter(a.qc_tier for a in all_arrays)
    category_counts = Counter(a.category for a in all_arrays)
    cas_status_counts = Counter(l["cas_link_status"] for l in links)

    summary = {
        "manifest_rows": len(manifest_rows),
        "genomes_with_identify_output": len(manifest_rows) - len(missing_identify),
        "missing_identify": missing_identify,
        "genomes_with_cctyper_output":
            len(manifest_rows) - len(missing_cctyper),
        "missing_cctyper": missing_cctyper,
        "candidates_total": len(all_arrays),
        "coordinate_check_fail":
            sum(1 for a in all_arrays if not a.coordinate_check_pass),
        "qc_tier_counts": dict(tier_counts),
        "category_counts": dict(category_counts),
        "hq_arrays": len(hq),
        "reserve_arrays": len(reserve),
        "nr100_repeats": fam_summary["nr100_repeats"],
        "families": fam_summary["families"],
        "cas_link_status_counts": dict(cas_status_counts),
        "known_array_regions_bed_rows": n_bed,
        "report_md": str(out_dir / "report.md"),
        "outputs": {
            "arrays.all.tsv": str(out_dir / "arrays.all.tsv"),
            "arrays.hq.tsv": str(out_dir / "arrays.hq.tsv"),
            "arrays.reserve.tsv": str(out_dir / "arrays.reserve.tsv"),
            "array_units.tsv": str(out_dir / "array_units.tsv"),
            "array_cas_links.tsv": str(out_dir / "array_cas_links.tsv"),
            "known_array_regions.bed": str(out_dir / "known_array_regions.bed"),
            "run_status.tsv": str(out_dir / "run_status.tsv"),
            "repeats.hq.nr100.fna": str(out_dir / "repeats.hq.nr100.fna"),
            "repeat_sources.tsv": str(out_dir / "repeat_sources.tsv"),
            "repeat_family_members.tsv": str(out_dir / "repeat_family_members.tsv"),
            "repeat_families.tsv": str(out_dir / "repeat_families.tsv"),
        },
    }
    report_mod.build_report(p, summary, out_dir / "report.md")
    print(json.dumps(summary, indent=2))
    # Exit 0 even if a small number of coord-check failures were isolated —
    # per step1.md §4 we treat upstream-bug arrays as flagged, not fatal.
    return 0


def cmd_run(cfg: dict, limit: int | None, resume: bool) -> int:
    raise SystemExit("run: not implemented yet — see step1.md sections 3-9")


def cmd_validate(cfg: dict) -> int:
    raise SystemExit("validate: not implemented yet — see step1.md section 12")


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="step1.py")
    sub = parser.add_subparsers(dest="cmd", required=True)

    def _add_config(p: argparse.ArgumentParser) -> None:
        p.add_argument("--config", required=True, type=Path)

    _add_config(sub.add_parser("doctor"))
    prep = sub.add_parser("prepare")
    _add_config(prep)
    prep.add_argument("--limit", type=int, default=None,
                      help="stage only the first N manifest rows")
    sm = sub.add_parser("submit")
    _add_config(sm)
    sm.add_argument("--skip-existing", action="store_true",
                    help="skip genomes whose raw output dir already exists")
    pr = sub.add_parser("parse")
    _add_config(pr)
    ex = sub.add_parser("export")
    _add_config(ex)
    ex.add_argument("--min-score", type=float, default=None,
                    help="override seed_qc.min_score (does not touch raw/)")
    ex.add_argument("--out-subdir", type=str, default="step1",
                    help="write under results/<out-subdir>/ (default: step1)")
    ex.add_argument("--identify-source", type=str, default="raw/identify",
                    help="which identify output tree to parse "
                         "(default: raw/identify; patched run uses raw_patched/identify)")
    r = sub.add_parser("run")
    _add_config(r)
    r.add_argument("--limit", type=int, default=None)
    r.add_argument("--resume", action="store_true")
    _add_config(sub.add_parser("validate"))

    args = parser.parse_args(argv)
    cfg = load_config(args.config)

    if args.cmd == "doctor":
        return cmd_doctor(cfg)
    if args.cmd == "prepare":
        return cmd_prepare(cfg, args.limit)
    if args.cmd == "submit":
        return cmd_submit(cfg, args.skip_existing)
    if args.cmd == "parse":
        return cmd_parse(cfg)
    if args.cmd == "export":
        return cmd_export(cfg, args.min_score, args.out_subdir,
                          identify_source=args.identify_source)
    if args.cmd == "run":
        return cmd_run(cfg, args.limit, args.resume)
    if args.cmd == "validate":
        return cmd_validate(cfg)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
