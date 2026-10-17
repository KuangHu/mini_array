#!/usr/bin/env python3
"""Phase L3: pure sequence export from the frozen coordinate catalog.

Writes everything under results/sequence_catalog/. The frozen
results/coordinate_catalog/ is NOT modified.

Non-negotiables encoded in code (per the 2026-10-16 instruction):
  - Only source of sequence: input/staged/<genome_uid>.fna
  - No consensus / no query-sequence substitution
  - All output is genomic-forward DNA; no automatic RC based on strand hypothesis
  - Preserve N and other ambiguous bases exactly as present in staged FASTA
  - Preserve existing count levels: 70,207 / 1,205 / 40 / 11 / 210
  - IDs stable and unique per file; map 1:1 back to the frozen coordinate rows
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _slurm_gate import require_slurm  # noqa: E402

PROJECT = Path("/groups/rubin/projects/kuang/out/CRISPR_array_related")
STAGED_DIR = PROJECT / "input" / "staged"
COORD = PROJECT / "results" / "coordinate_catalog"
OUT = PROJECT / "results" / "sequence_catalog"


# ----------------- fasta IO helpers -----------------

def _read_staged_fasta(path: Path) -> dict[str, str]:
    seqs = {}
    cur = None
    parts: list[str] = []
    with path.open() as fh:
        for ln in fh:
            if ln.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(parts)
                cur = ln[1:].split()[0]
                parts = []
            else:
                parts.append(ln.strip())
        if cur is not None:
            seqs[cur] = "".join(parts)
    return seqs


class _StagedCache:
    def __init__(self):
        self._cache: dict[str, dict[str, str]] = {}

    def get(self, genome_uid: str) -> dict[str, str]:
        if genome_uid not in self._cache:
            self._cache[genome_uid] = _read_staged_fasta(
                STAGED_DIR / f"{genome_uid}.fna"
            )
        return self._cache[genome_uid]


# ----------------- fingerprints -----------------

def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


# ----------------- contig lookup (genome_id, genome_uid) -----------------

def _load_contig_lookup() -> dict[str, dict]:
    out: dict[str, dict] = {}
    with (PROJECT / "results" / "staging" / "contigs.1000.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["contig_uid"]] = {
                "genome_id": r["genome_id"],
                "genome_uid": r["genome_uid"],
                "length": int(r["length"]),
            }
    return out


# ----------------- fasta writer with index entries -----------------

class _FastaWriter:
    def __init__(self, path: Path, cache: _StagedCache, contigs: dict,
                  fasta_label: str):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path
        self.fh = path.open("w")
        self.cache = cache
        self.contigs = contigs
        self.label = fasta_label
        self.n_written = 0
        self.n_oob = 0
        self.index_rows: list[dict] = []

    def write_record(self, *, record_id: str,
                      contig_uid: str, start0: int, end0: int, meta: dict
                      ) -> bool:
        meta_contig = self.contigs.get(contig_uid)
        if meta_contig is None:
            self.n_oob += 1
            return False
        guid = meta_contig["genome_uid"]
        seqs = self.cache.get(guid)
        contig_seq = seqs.get(contig_uid)
        if contig_seq is None or start0 < 0 or end0 > len(contig_seq) or end0 <= start0:
            self.n_oob += 1
            return False
        seq = contig_seq[start0:end0]
        self.fh.write(f">{record_id}\n{seq}\n")
        self.n_written += 1
        has_amb = "true" if any(c not in "ACGT" for c in seq) else "false"
        row = {
            "fasta_file": self.label,
            "record_id": record_id,
            "genome_id": meta_contig["genome_id"],
            "contig_uid": contig_uid,
            "start0": start0,
            "end0": end0,
            "length": end0 - start0,
            "sequence_orientation": "genomic_forward",
            "has_ambiguous_bases": has_amb,
            "sequence_sha256": _sha256_text(seq),
            **meta,
        }
        self.index_rows.append(row)
        return True

    def close(self) -> None:
        self.fh.close()


# ----------------- exporters -----------------

def export_arrays(fname: str, src_tsv: Path, cache: _StagedCache,
                   contigs: dict, index_sink: list[dict]) -> int:
    """Emit <fname> under OUT with one FASTA record per array row in src_tsv.
    IDs are '<source_tool>|<source_id>' so they stay unique in both the all
    and the trusted exports."""
    writer = _FastaWriter(OUT / fname, cache, contigs, fname)
    with src_tsv.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            # Skip rows whose coord_validation flagged them not in bounds
            if r.get("coordinate_in_bounds", "true") != "true":
                continue
            rec_id = f"{r['source_tool']}|{r['source_id']}"
            writer.write_record(
                record_id=rec_id,
                contig_uid=r["contig_uid"],
                start0=int(r["start0"]), end0=int(r["end0"]),
                meta={
                    "call_id": r["source_id"],
                    "physical_locus_id": "",
                    "pair_id": "",
                    "unit_type": "array",
                    "locus_class": r["category_raw"],
                },
            )
    writer.close()
    index_sink.extend(writer.index_rows)
    print(f"[{fname}] wrote {writer.n_written} records, {writer.n_oob} skipped")
    return writer.n_written


def export_mini_physical(fname: str, src_tsv: Path, cache: _StagedCache,
                          contigs: dict, index_sink: list[dict]) -> int:
    writer = _FastaWriter(OUT / fname, cache, contigs, fname)
    with src_tsv.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("coordinate_in_bounds", "true") != "true":
                continue
            writer.write_record(
                record_id=r["physical_locus_id"],
                contig_uid=r["contig_uid"],
                start0=int(r["start0"]), end0=int(r["end0"]),
                meta={
                    "call_id": "",
                    "physical_locus_id": r["physical_locus_id"],
                    "pair_id": "",
                    "unit_type": "physical_locus",
                    "locus_class": r["structural_class"],
                },
            )
    writer.close()
    index_sink.extend(writer.index_rows)
    print(f"[{fname}] wrote {writer.n_written} records, {writer.n_oob} skipped")
    return writer.n_written


def export_units(fname: str, src_tsv: Path, cache: _StagedCache,
                  contigs: dict, index_sink: list[dict]) -> int:
    writer = _FastaWriter(OUT / fname, cache, contigs, fname)
    with src_tsv.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("coordinate_in_bounds", "true") != "true":
                continue
            pair_id = r["pair_id"]
            locus_class = r["locus_class"]
            for utype, s0, e0, extra_meta in (
                ("repeat_1", int(r["repeat_1_start0"]), int(r["repeat_1_end0"]),
                 {"unit_query_id": r.get("repeat_1_query_id", "")}),
                ("spacer", int(r["spacer_start0"]), int(r["spacer_end0"]),
                 {"spacer_len": r.get("spacer_len", "")}),
                ("repeat_2", int(r["repeat_2_start0"]), int(r["repeat_2_end0"]),
                 {"unit_query_id": r.get("repeat_2_query_id", ""),
                  "repeat_2_completeness": r.get("repeat_2_completeness", ""),
                  "repeat_2_pident": r.get("repeat_2_pident", ""),
                  "repeat_2_cov": r.get("repeat_2_cov", "")}),
            ):
                writer.write_record(
                    record_id=f"{pair_id}|{utype}",
                    contig_uid=r["contig_uid"],
                    start0=s0, end0=e0,
                    meta={
                        "call_id": "",
                        "physical_locus_id": r["physical_locus_id"],
                        "pair_id": pair_id,
                        "unit_type": utype,
                        "locus_class": locus_class,
                        **{k: v for k, v in extra_meta.items()},
                    },
                )
    writer.close()
    index_sink.extend(writer.index_rows)
    print(f"[{fname}] wrote {writer.n_written} records, {writer.n_oob} skipped")
    return writer.n_written


# ----------------- array units reference copy -----------------

def export_array_units_reference(index_sink: list[dict]) -> int:
    """Copy / organize the patched identify array_units.tsv into the
    sequence_catalog, and emit a cctyper stub table that labels each
    cctyper array as 'unit annotation unavailable'. We do NOT re-run or
    re-identify units."""
    src = PROJECT / "results" / "step1_patched_min_score_0.75" / "array_units.tsv"
    dst = OUT / "array_units_patched.tsv"
    rows = []
    with src.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(r)
    with dst.open("w", newline="") as fh:
        if rows:
            w = csv.DictWriter(fh,
                               fieldnames=list(rows[0].keys()),
                               delimiter="\t", lineterminator="\n")
            w.writeheader()
            for r in rows:
                w.writerow(r)
    n_rows = len(rows)
    print(f"[array_units_patched.tsv] copied {n_rows} rows from "
           f"step1_patched_min_score_0.75/array_units.tsv")

    # cctyper stub: one row per cctyper array row in crispr_arrays.tsv
    # that says 'unit annotation unavailable'.
    cct_stub = OUT / "array_units_cctyper_unavailable.tsv"
    stub_rows = []
    with (COORD / "crispr_arrays.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["source_tool"] != "crisprcastyper":
                continue
            stub_rows.append({
                "source_tool": "crisprcastyper",
                "source_id": r["source_id"],
                "genome_id": r["genome_id"],
                "contig_uid": r["contig_uid"],
                "start0": r["start0"],
                "end0": r["end0"],
                "unit_annotation": "unit annotation unavailable",
                "note": ("CRISPRCasTyper does not expose per-repeat/spacer unit "
                          "coordinates in crisprs_all.tab; cctyper's spacers/ "
                          "folder has spacer FASTAs but no per-unit genomic "
                          "coordinate table. Not re-identified here."),
            })
    with cct_stub.open("w", newline="") as fh:
        w = csv.DictWriter(fh,
                           fieldnames=["source_tool", "source_id", "genome_id",
                                        "contig_uid", "start0", "end0",
                                        "unit_annotation", "note"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in stub_rows:
            w.writerow(r)
    print(f"[array_units_cctyper_unavailable.tsv] stubbed {len(stub_rows)} rows")
    return n_rows


# ----------------- validation -----------------

def validate_catalog(index_rows: list[dict], expected_counts: dict,
                      contigs: dict, cache: _StagedCache) -> dict:
    """Verify every record_id, length, and sequence content."""
    issues = {"length_mismatch": 0, "seq_mismatch": 0,
               "missing_contig": 0, "unexpected_rows": {}}
    actual_counts: dict[str, int] = defaultdict(int)
    for r in index_rows:
        actual_counts[r["fasta_file"]] += 1
        meta = contigs.get(r["contig_uid"])
        if meta is None:
            issues["missing_contig"] += 1
            continue
        if int(r["length"]) != int(r["end0"]) - int(r["start0"]):
            issues["length_mismatch"] += 1
        # Verify sha256 against a fresh slice
        seqs = cache.get(meta["genome_uid"])
        slice_ = seqs[r["contig_uid"]][int(r["start0"]):int(r["end0"])]
        if hashlib.sha256(slice_.encode()).hexdigest() != r["sequence_sha256"]:
            issues["seq_mismatch"] += 1
    for name, exp in expected_counts.items():
        got = actual_counts.get(name, 0)
        if got != exp:
            issues["unexpected_rows"][name] = {"expected": exp, "actual": got}
    issues["actual_counts"] = dict(actual_counts)
    return issues


# ----------------- fingerprints writer -----------------

def write_fingerprints(coord_dir: Path) -> dict:
    """Record SHA256 of the input coordinate catalog files + a listing of
    staged FASTAs (hash the first few + the total manifest)."""
    inputs = {}
    for name in ("crispr_arrays.tsv", "crispr_arrays_trusted.tsv",
                  "mini_arrays.tsv", "edge_unresolved_candidates.tsv",
                  "mini_array_units.tsv",
                  "coordinate_catalog_summary.md"):
        p = coord_dir / name
        if p.is_file():
            inputs[name] = {"path": str(p), "sha256": _sha256_file(p),
                             "size": p.stat().st_size}
    # Staged FASTA fingerprint: hash concatenation of (filename:size) listing.
    staged_listing = []
    for f in sorted(STAGED_DIR.glob("*.fna")):
        staged_listing.append(f"{f.name}:{f.stat().st_size}")
    listing_text = "\n".join(staged_listing)
    staged_fp = {
        "n_files": len(staged_listing),
        "listing_sha256": _sha256_text(listing_text),
        "note": ("Listing-level fingerprint (filename:size) over all 1000 "
                 "staged FASTAs. Individual per-file hashes omitted for "
                 "speed; each staged FASTA is deterministic from Step 1's "
                 "staging code (see provenance/LOCKED.md).")
    }
    return {"input_coordinate_catalog": inputs, "input_staged_fastas": staged_fp}


# ----------------- README writer -----------------

def write_readme(counts: dict, validation: dict, fingerprints: dict) -> None:
    md = []
    md.append("# Sequence catalog (frozen)\n")
    md.append("Pure genomic-forward DNA export from the frozen coordinate "
               "catalog. No new biology, no new search, no re-identification. "
               "All sequences are literal slices of `input/staged/<genome_uid>."
               "fna` using the 0-based half-open coordinates in "
               "`results/coordinate_catalog/`. Ambiguous bases (N and other "
               "non-ACGT) are preserved as-is.\n")
    md.append("## Scope\n")
    md.append("- This catalog covers ONLY the frozen 1000-genome batch and "
               "the calls that were already recovered under the current "
               "coordinate catalog (allow_full_run: false).")
    md.append("- `crispr_arrays.all_calls.fna` is the union of every call "
               "from CRISPRidentify + CRISPRCasTyper. **It is NOT a set of "
               "independent trusted arrays**; rows from the two tools may "
               "overlap on the same genomic interval.")
    md.append("- `crispr_arrays.trusted_calls.fna` applies the "
               "coordinate_catalog's trusted filter (identify Bona-fide OR "
               "cctyper Trusted) and is still per-call. Downstream work must "
               "dedup if an independent-locus count is needed.")
    md.append("- 210 unit records = 70 pairs × 3 unit_types (repeat_1, "
               "spacer, repeat_2), spanning both `candidate_mini_array` and "
               "`edge_unresolved` loci. The `locus_class` column in "
               "`sequence_index.tsv` partitions them.")
    md.append("- The existing degenerate-repeat coverage limitation from "
               "Step 2's secondary thresholds remains unchanged (see the "
               "coordinate_catalog summary for the note).")
    md.append("")

    md.append("## Expected vs actual record counts (post-write validation)\n")
    md.append("| fasta_file | expected | actual |")
    md.append("|---|---:|---:|")
    for k in ("crispr_arrays.all_calls.fna",
               "crispr_arrays.trusted_calls.fna",
               "mini_arrays.fna",
               "edge_unresolved_candidates.fna",
               "mini_array_units.fna"):
        exp = counts["expected"].get(k, "—")
        act = validation["actual_counts"].get(k, 0)
        md.append(f"| `{k}` | {exp} | {act} |")
    md.append("")

    md.append("### Validation results\n")
    md.append(f"- missing contig lookups: {validation['missing_contig']}")
    md.append(f"- length mismatches (len != end0 - start0): {validation['length_mismatch']}")
    md.append(f"- sequence sha256 mismatches: {validation['seq_mismatch']}")
    if validation["unexpected_rows"]:
        md.append("- COUNT MISMATCHES:")
        for k, v in validation["unexpected_rows"].items():
            md.append(f"  - `{k}`: expected {v['expected']}, actual {v['actual']}")
    else:
        md.append("- all per-file counts match expected values")
    md.append("")

    md.append("## FASTA ID conventions\n")
    md.append("| file | ID format |")
    md.append("|---|---|")
    md.append("| crispr_arrays.all_calls.fna | `<source_tool>|<source_id>` |")
    md.append("| crispr_arrays.trusted_calls.fna | `<source_tool>|<source_id>` |")
    md.append("| mini_arrays.fna | `<physical_locus_id>` |")
    md.append("| edge_unresolved_candidates.fna | `<physical_locus_id>` |")
    md.append("| mini_array_units.fna | `<pair_id>|<unit_type>` where unit_type ∈ {repeat_1, spacer, repeat_2} |")
    md.append("")
    md.append("Every ID is unique within its own FASTA file. `sequence_index.tsv` "
               "maps every ID back to its row in the frozen coordinate catalog.")
    md.append("")

    md.append("## Deliverables\n")
    md.append("```")
    md.append("results/sequence_catalog/")
    md.append("  crispr_arrays.all_calls.fna")
    md.append("  crispr_arrays.trusted_calls.fna")
    md.append("  mini_arrays.fna")
    md.append("  edge_unresolved_candidates.fna")
    md.append("  mini_array_units.fna")
    md.append("  array_units_patched.tsv            # identify repeat/spacer units (copy of patched output)")
    md.append("  array_units_cctyper_unavailable.tsv # cctyper calls: 'unit annotation unavailable'")
    md.append("  sequence_index.tsv")
    md.append("  input_fingerprints.json            # coordinate catalog + staged FASTA fingerprints")
    md.append("  README.md")
    md.append("```")
    md.append("")

    md.append("## Non-negotiables honoured\n")
    md.append("- Sole source of sequence: `input/staged/<genome_uid>.fna`. No "
               "consensus/query substitution. No re-identification of units.")
    md.append("- All sequence output is genomic-forward DNA; no automatic "
               "reverse-complement based on strand hypothesis. Strand "
               "hypotheses stay in `coordinate_catalog/mini_arrays.tsv`.")
    md.append("- Ambiguous bases preserved as-is; `has_ambiguous_bases` column "
               "in `sequence_index.tsv` flags any non-ACGT content.")
    md.append("- Count levels preserved. Identical unit sequences are NOT "
               "deduped — pair-level correspondence must be kept.")
    md.append("- No parameter / threshold / boundary change. No new search. "
               "No CDS / protein / HMM analysis.")
    md.append("- All execution via sbatch; SLURM_JOB_ID gate enforced at "
               "import.")
    md.append("- `allow_full_run: false` preserved.")
    md.append("")

    md.append("## Input fingerprints\n")
    md.append("```json")
    md.append(json.dumps(fingerprints, indent=2))
    md.append("```")
    (OUT / "README.md").write_text("\n".join(md) + "\n")


# ----------------- main -----------------

def main() -> int:
    require_slurm()
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"OUT = {OUT}")

    print("[fingerprint] hashing coordinate catalog + staged listing…")
    fingerprints = write_fingerprints(COORD)
    (OUT / "input_fingerprints.json").write_text(
        json.dumps(fingerprints, indent=2) + "\n"
    )

    print("[load] contigs lookup…")
    contigs = _load_contig_lookup()
    cache = _StagedCache()

    index_rows: list[dict] = []
    counts = {"expected": {
        "crispr_arrays.all_calls.fna": 70207,
        "crispr_arrays.trusted_calls.fna": 1205,
        "mini_arrays.fna": 40,
        "edge_unresolved_candidates.fna": 11,
        "mini_array_units.fna": 210,
    }, "actual": {}}

    counts["actual"]["crispr_arrays.all_calls.fna"] = export_arrays(
        "crispr_arrays.all_calls.fna",
        COORD / "crispr_arrays.tsv",
        cache, contigs, index_rows,
    )
    counts["actual"]["crispr_arrays.trusted_calls.fna"] = export_arrays(
        "crispr_arrays.trusted_calls.fna",
        COORD / "crispr_arrays_trusted.tsv",
        cache, contigs, index_rows,
    )
    counts["actual"]["mini_arrays.fna"] = export_mini_physical(
        "mini_arrays.fna",
        COORD / "mini_arrays.tsv",
        cache, contigs, index_rows,
    )
    counts["actual"]["edge_unresolved_candidates.fna"] = export_mini_physical(
        "edge_unresolved_candidates.fna",
        COORD / "edge_unresolved_candidates.tsv",
        cache, contigs, index_rows,
    )
    counts["actual"]["mini_array_units.fna"] = export_units(
        "mini_array_units.fna",
        COORD / "mini_array_units.tsv",
        cache, contigs, index_rows,
    )

    export_array_units_reference(index_rows)

    # Write sequence_index.tsv
    index_cols = ["fasta_file", "record_id", "genome_id", "contig_uid",
                   "call_id", "physical_locus_id", "pair_id", "unit_type",
                   "locus_class", "start0", "end0", "length",
                   "sequence_orientation", "has_ambiguous_bases",
                   "sequence_sha256",
                   "unit_query_id", "repeat_2_completeness",
                   "repeat_2_pident", "repeat_2_cov", "spacer_len"]
    with (OUT / "sequence_index.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=index_cols, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in index_rows:
            # fill missing fields with "" (not NA — NA could be mis-parsed as
            # a sodium name; empty = clearly not applicable)
            w.writerow({c: r.get(c, "") for c in index_cols})

    print("[validate] re-reading each slice and comparing sha256…")
    validation = validate_catalog(index_rows, counts["expected"],
                                   contigs, cache)
    print(f"[validate] missing_contig={validation['missing_contig']} "
           f"length_mismatch={validation['length_mismatch']} "
           f"seq_mismatch={validation['seq_mismatch']} "
           f"unexpected={validation['unexpected_rows']}")

    write_readme(counts, validation, fingerprints)
    (OUT / "validation.json").write_text(
        json.dumps({"expected": counts["expected"],
                     "actual_counts": validation["actual_counts"],
                     "issues": {k: v for k, v in validation.items()
                                 if k != "actual_counts"}}, indent=2) + "\n"
    )

    dt = time.time() - t0
    print(f"[done] {dt:.1f}s; wrote {len(index_rows)} index rows")
    return 0 if (not validation["unexpected_rows"]
                 and validation["length_mismatch"] == 0
                 and validation["seq_mismatch"] == 0) else 1


if __name__ == "__main__":
    sys.exit(main())
