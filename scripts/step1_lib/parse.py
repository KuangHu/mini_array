"""Parse CRISPRidentify output into audited tables.

Per step1.md §4, the JSON top-level `start/end/length` are NOT reliable BED
coordinates (we have confirmed on real output that for Reversed arrays `end`
reports the smaller coordinate side and `length` can even be negative). The
only authoritative fields are the paired lists:

  list_repeats            — repeat sequences in DISPLAY orientation
  list_repeat_starts      — 0-based genomic-forward start of each repeat
                            (ascending for Forward arrays, descending for
                            Reversed arrays, but the integers are forward coords)
  list_spacers            — spacer sequences in DISPLAY orientation
  consensus_repeat        — repeat consensus in DISPLAY orientation
  strand                  — "Forward (Orientation was not computed)" | "Reversed"
                            | "Forward" | similar

Per `components_output_maker.py` the display orientation equals the forward
genomic strand iff strand starts with "Forward"; otherwise sequences in the
lists are reverse-complements of the forward reference.

For each candidate we rebuild genomic-forward repeat/spacer intervals, verify
every base against the staged FASTA, and only then emit a parsed record.
Failures do not get silently fixed.
"""
from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn-",
                           "TGCAYRWSMKVHDBNtgcayrwsmkvhdbn-")

CATEGORY_MAP = {
    "Bona-fide": "bona_fide",
    "Alternative": "alternative",
    "Possible": "possible",
    "Possible Discarded": "possible_discarded",
    "Low score": "low_score",
}

CATEGORIES = list(CATEGORY_MAP.keys())


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def canonical(seq: str) -> str:
    rc = revcomp(seq)
    return seq if seq <= rc else rc


def repeat_id(forward_consensus: str) -> str:
    canon = canonical(forward_consensus.upper())
    return "rep_" + hashlib.sha256(canon.encode()).hexdigest()[:20]


def array_id(genome_id: str, contig_uid: str, start0: int, end0: int,
             consensus_canonical: str) -> str:
    key = f"{genome_id}|{contig_uid}|{start0}|{end0}|{consensus_canonical}".encode()
    return "arr_" + hashlib.sha256(key).hexdigest()[:20]


def orientation_label(strand_raw: str) -> str:
    s = (strand_raw or "").strip()
    if s.startswith("Forward (Orientation was not computed)"):
        return "unknown"
    if s.startswith("Forward"):
        return "forward"
    if s.startswith("Reversed"):
        return "reversed"
    return "unknown"


def display_is_forward(strand_raw: str) -> bool:
    """True if the display orientation equals the genomic forward strand."""
    return (strand_raw or "").startswith("Forward")


@dataclass
class RepeatUnit:
    unit_index: int
    start0: int
    end0: int
    sequence_forward: str


@dataclass
class SpacerUnit:
    unit_index: int
    start0: int
    end0: int
    sequence_forward: str


@dataclass
class ParsedArray:
    genome_id: str
    contig_uid: str
    locus_id: str          # per-contig region index from upstream
    call_id: str           # category + locus + array_index — unique within genome
    array_id: str          # stable hash
    category_raw: str
    category: str
    score: float
    strand_raw: str
    orientation_prediction: str
    consensus_reported: str      # as output by upstream (display orientation)
    consensus_forward: str       # genomic-forward orientation
    consensus_canonical: str
    repeat_id: str
    start0: int                  # ascending-sorted
    end0: int                    # exclusive
    length: int
    n_repeats: int
    n_spacers: int
    n_unique_spacers: int
    repeats: list[RepeatUnit] = field(default_factory=list)
    spacers: list[SpacerUnit] = field(default_factory=list)
    coordinate_check_pass: bool = False
    qc_reasons: list[str] = field(default_factory=list)


def load_contigs_index(contigs_tsv: Path) -> dict[str, dict]:
    """Return {contig_uid: {...}} from contigs.tsv."""
    out: dict[str, dict] = {}
    with contigs_tsv.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["contig_uid"]] = row
    return out


def read_staged_fasta(path: Path) -> dict[str, str]:
    """Return {contig_uid: sequence (uppercase, single-line)}."""
    out: dict[str, str] = {}
    cur: str | None = None
    parts: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    out[cur] = "".join(parts)
                cur = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if cur is not None:
            out[cur] = "".join(parts)
    return out


def parse_candidate(category_raw: str, locus_key: str, cand: dict,
                    contig_uid: str, contig_seq: str,
                    genome_id: str) -> ParsedArray:
    """Build a ParsedArray from one JSON candidate; verify every base."""
    strand_raw = cand.get("strand", "")
    forward_display = display_is_forward(strand_raw)

    raw_repeats: list[str] = list(cand.get("list_repeats", []))
    raw_starts: list[int] = list(cand.get("list_repeat_starts", []))
    raw_spacers: list[str] = list(cand.get("list_spacers", []))
    consensus_display = cand.get("consensus_repeat", "") or ""
    score = float(cand.get("certainty_score", "nan"))
    array_index = cand.get("array_index", locus_key)

    reasons: list[str] = []
    pairs: list[tuple[int, str]] = []
    for i, (start_fwd, rep_disp) in enumerate(zip(raw_starts, raw_repeats)):
        rep_fwd = rep_disp if forward_display else revcomp(rep_disp)
        pairs.append((int(start_fwd), rep_fwd))

    pairs.sort(key=lambda p: p[0])

    repeats: list[RepeatUnit] = []
    for idx, (s0, rep_fwd) in enumerate(pairs):
        e0 = s0 + len(rep_fwd)
        if s0 < 0 or e0 > len(contig_seq):
            reasons.append(f"repeat {idx} out of contig bounds: [{s0}:{e0})")
        elif contig_seq[s0:e0] != rep_fwd:
            reasons.append(f"repeat {idx} does not match reference at [{s0}:{e0})")
        repeats.append(RepeatUnit(unit_index=idx, start0=s0, end0=e0,
                                  sequence_forward=rep_fwd))

    # Rebuild spacers from reference between successive sorted repeats.
    spacers: list[SpacerUnit] = []
    for idx in range(len(repeats) - 1):
        s0 = repeats[idx].end0
        e0 = repeats[idx + 1].start0
        if e0 < s0:
            reasons.append(f"spacer {idx} has non-positive length ({s0}..{e0})")
            spacer_fwd = ""
        else:
            spacer_fwd = contig_seq[s0:e0]
        spacers.append(SpacerUnit(unit_index=idx, start0=s0, end0=e0,
                                  sequence_forward=spacer_fwd))

    # Cross-check against upstream's spacer list (order differs for Reversed).
    if raw_spacers:
        upstream_fwd = [s if forward_display else revcomp(s) for s in raw_spacers]
        # For reversed display, upstream spacers come top-down in display order,
        # which maps to descending genomic coords, i.e. reverse of sorted order.
        expected = upstream_fwd if forward_display else list(reversed(upstream_fwd))
        got = [sp.sequence_forward for sp in spacers]
        if got != expected:
            reasons.append("spacer sequences do not match upstream list_spacers")

    if repeats:
        start0 = repeats[0].start0
        end0 = repeats[-1].end0
    else:
        start0 = end0 = -1
        reasons.append("no repeats in candidate")

    length = end0 - start0 if end0 > start0 else 0
    consensus_forward = consensus_display if forward_display else revcomp(consensus_display)
    consensus_canon = canonical(consensus_forward.upper())
    rep_id = repeat_id(consensus_forward)

    locus_id = f"{contig_uid}:{locus_key}"
    cat_norm = CATEGORY_MAP.get(category_raw, category_raw.lower().replace(" ", "_"))
    call_id = f"{contig_uid}:{cat_norm}:{array_index}"

    coord_ok = len(reasons) == 0
    n_rep = len(repeats)
    n_spacer = len(spacers)
    n_unique_spacer = len({sp.sequence_forward for sp in spacers if sp.sequence_forward})

    return ParsedArray(
        genome_id=genome_id,
        contig_uid=contig_uid,
        locus_id=locus_id,
        call_id=call_id,
        array_id=array_id(genome_id, contig_uid, start0, end0, consensus_canon),
        category_raw=category_raw,
        category=cat_norm,
        score=score,
        strand_raw=strand_raw,
        orientation_prediction=orientation_label(strand_raw),
        consensus_reported=consensus_display,
        consensus_forward=consensus_forward,
        consensus_canonical=consensus_canon,
        repeat_id=rep_id,
        start0=start0,
        end0=end0,
        length=length,
        n_repeats=n_rep,
        n_spacers=n_spacer,
        n_unique_spacers=n_unique_spacer,
        repeats=repeats,
        spacers=spacers,
        coordinate_check_pass=coord_ok,
        qc_reasons=reasons,
    )


def parse_genome(genome_id: str, genome_uid: str, raw_identify_dir: Path,
                 staged_fasta: Path) -> list[ParsedArray]:
    json_summary = raw_identify_dir / "Complete_json_summary.json"
    if not json_summary.is_file():
        return []

    with json_summary.open() as fh:
        summary = json.load(fh)

    seqs = read_staged_fasta(staged_fasta)
    arrays: list[ParsedArray] = []
    for contig_uid, blob in summary.items():
        contig_seq = seqs.get(contig_uid)
        if contig_seq is None:
            # upstream can emit entries for contigs we didn't stage; skip
            continue
        buckets = blob.get("Arrays", {})
        for category_raw in CATEGORIES:
            cat_blob = buckets.get(category_raw, {}) or {}
            for locus_key, cand in cat_blob.items():
                arrays.append(parse_candidate(
                    category_raw=category_raw,
                    locus_key=str(locus_key),
                    cand=cand,
                    contig_uid=contig_uid,
                    contig_seq=contig_seq,
                    genome_id=genome_id,
                ))
    return arrays


ARRAYS_COLS = [
    "array_id", "call_id", "locus_id", "genome_id", "contig_uid",
    "start0", "end0", "length",
    "category_raw", "category", "score",
    "n_repeats", "n_spacers", "n_unique_spacers",
    "consensus_reported", "consensus_forward", "consensus_canonical",
    "repeat_id", "strand_raw", "orientation_prediction",
    "coordinate_check_pass", "qc_reasons",
]

UNITS_COLS = [
    "array_id", "contig_uid", "unit_kind", "unit_index",
    "start0", "end0", "sequence_forward",
]

ARRAYS_COLS_FULL = ARRAYS_COLS + [
    "qc_tier", "n_full_repeats", "left_flank_bp", "right_flank_bp",
    "ambiguous_bases_in_context", "terminal_repeat_degenerate",
]


def write_arrays_tsv_full(arrays: list[ParsedArray], out_path: Path) -> None:
    """Variant of write_arrays_tsv that also serializes QC columns attached as
    attributes on each ParsedArray (qc_tier, n_full_repeats, ...)."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(".tsv.tmp")
    with tmp.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(ARRAYS_COLS_FULL)
        for a in arrays:
            w.writerow([
                a.array_id, a.call_id, a.locus_id, a.genome_id, a.contig_uid,
                a.start0, a.end0, a.length,
                a.category_raw, a.category, f"{a.score:.6f}",
                a.n_repeats, a.n_spacers, a.n_unique_spacers,
                a.consensus_reported, a.consensus_forward, a.consensus_canonical,
                a.repeat_id, a.strand_raw, a.orientation_prediction,
                "true" if a.coordinate_check_pass else "false",
                json.dumps(a.qc_reasons, separators=(",", ":")),
                getattr(a, "qc_tier", ""),
                getattr(a, "n_full_repeats", ""),
                getattr(a, "left_flank_bp", ""),
                getattr(a, "right_flank_bp", ""),
                getattr(a, "ambiguous_bases_in_context", ""),
                "true" if getattr(a, "terminal_repeat_degenerate", False) else "false",
            ])
    tmp.rename(out_path)


BED_COLS = ["contig_uid", "start0", "end0", "name", "score", "strand"]


def write_known_array_regions_bed(arrays: list[ParsedArray], out_path: Path) -> int:
    """BED6 of bona-fide + reserve arrays per step1.md §9.
    `name` is array_id, `score` is UCSC 0-1000 scaled from certainty,
    `strand` reflects orientation_prediction (+/-/?)."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(".bed.tmp")
    n = 0
    with tmp.open("w") as fh:
        for a in sorted(arrays, key=lambda x: (x.contig_uid, x.start0, x.end0)):
            if a.category_raw != "Bona-fide" or not a.coordinate_check_pass:
                continue
            if a.start0 < 0 or a.end0 <= a.start0:
                continue
            strand = {"forward": "+", "reversed": "-", "unknown": "?"}.get(
                a.orientation_prediction, "?"
            )
            score = max(0, min(1000, int(round(a.score * 1000))))
            fh.write("\t".join([
                a.contig_uid, str(a.start0), str(a.end0),
                a.array_id, str(score), strand,
            ]) + "\n")
            n += 1
    tmp.rename(out_path)
    return n


def write_arrays_tsv(arrays: list[ParsedArray], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(".tsv.tmp")
    with tmp.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(ARRAYS_COLS)
        for a in arrays:
            w.writerow([
                a.array_id, a.call_id, a.locus_id, a.genome_id, a.contig_uid,
                a.start0, a.end0, a.length,
                a.category_raw, a.category, f"{a.score:.6f}",
                a.n_repeats, a.n_spacers, a.n_unique_spacers,
                a.consensus_reported, a.consensus_forward, a.consensus_canonical,
                a.repeat_id, a.strand_raw, a.orientation_prediction,
                "true" if a.coordinate_check_pass else "false",
                json.dumps(a.qc_reasons, separators=(",", ":")),
            ])
    tmp.rename(out_path)


def write_units_tsv(arrays: list[ParsedArray], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(".tsv.tmp")
    with tmp.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(UNITS_COLS)
        for a in arrays:
            for r in a.repeats:
                w.writerow([a.array_id, a.contig_uid, "repeat",
                            r.unit_index, r.start0, r.end0, r.sequence_forward])
            for s in a.spacers:
                w.writerow([a.array_id, a.contig_uid, "spacer",
                            s.unit_index, s.start0, s.end0, s.sequence_forward])
    tmp.rename(out_path)


def parse_all(project_root: Path, manifest_rows: list[dict],
              arrays_tsv: Path, units_tsv: Path) -> dict:
    all_arrays: list[ParsedArray] = []
    missing: list[str] = []
    for row in manifest_rows:
        gid = row["genome_id"]
        guid = row.get("genome_uid") or _sanitize(gid)
        raw_dir = project_root / "raw" / "identify" / guid
        staged = project_root / "input" / "staged" / f"{guid}.fna"
        if not (raw_dir / "Complete_json_summary.json").is_file():
            missing.append(gid)
            continue
        all_arrays.extend(parse_genome(gid, guid, raw_dir, staged))

    write_arrays_tsv(all_arrays, arrays_tsv)
    write_units_tsv(all_arrays, units_tsv)
    verified = sum(1 for a in all_arrays if a.coordinate_check_pass)
    return {
        "genomes_parsed": len(manifest_rows) - len(missing),
        "genomes_missing_output": len(missing),
        "missing_genomes": missing,
        "candidates_total": len(all_arrays),
        "coordinate_check_pass": verified,
        "coordinate_check_fail": len(all_arrays) - verified,
        "arrays_tsv": str(arrays_tsv),
        "units_tsv": str(units_tsv),
    }


def _sanitize(name: str) -> str:
    # Local re-import to avoid a cycle with stage.py at import time.
    from . import stage as _stage
    return _stage.sanitize(name)
