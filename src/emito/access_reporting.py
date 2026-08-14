#!/usr/bin/env python3
"""Replay eMito-access and report how many probes each filter removes."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from . import pipeline


MODE_ROOTS = {
    "taxa": Path("03_eMito_taxa_generate/probe_sets"),
    "group": Path("04_eMito_group_generate/probe_sets"),
    "node": Path("05_eMito_node_generate/probe_sets"),
}


class StatisticsError(RuntimeError):
    pass


@dataclass(frozen=True)
class Parameters:
    gc_min: float
    gc_max: float
    complexity_min: float
    complexity_max: float
    dimer_k: int
    dimer_threshold: float


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Replay one mode's eMito-access step and summarize sequential GC, "
            "complexity, dimer, and exact-sequence deduplication counts."
        ),
    )
    parser.add_argument("--output-root", type=Path, default=Path("emito_output"))
    parser.add_argument("--mode", choices=tuple(MODE_ROOTS), default="taxa")
    parser.add_argument(
        "--report-dir",
        type=Path,
        help=(
            "Output directory; default is "
            "<output-root>/08_probe_statistics/<mode>_access_filters."
        ),
    )
    return parser.parse_args(argv)


def load_config(path: Path) -> Mapping[str, object]:
    if not path.is_file():
        raise StatisticsError(f"Pipeline configuration does not exist: {path}")
    with path.open("rt", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise StatisticsError(f"Pipeline configuration is not a JSON object: {path}")
    return value


def load_parameters(config: Mapping[str, object]) -> Parameters:
    if "dimer_threshold" not in config:
        raise StatisticsError(
            "This run predates the eProbe-compatible dimer_threshold algorithm; "
            "rerun it with eMito 0.1.1 or later before using access-summary."
        )
    return Parameters(
        gc_min=float(config.get("gc_min", 35.0)),
        gc_max=float(config.get("gc_max", 65.0)),
        complexity_min=float(config.get("complexity_min", 0.0)),
        complexity_max=float(config.get("complexity_max", 2.0)),
        dimer_k=int(config.get("dimer_k", 11)),
        dimer_threshold=float(config["dimer_threshold"]),
    )


def load_manifest(path: Path) -> List[Dict[str, str]]:
    if not path.is_file():
        raise StatisticsError(f"Generation manifest does not exist: {path}")
    with path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"target_id", "rank", "name", "taxid", "output_fasta"}
    if not rows:
        raise StatisticsError(f"Generation manifest is empty: {path}")
    missing = required - set(rows[0])
    if missing:
        raise StatisticsError(f"Manifest lacks columns {sorted(missing)}: {path}")
    return rows


def resolve_generation_fasta(root: Path, mode: str, raw_path: str) -> Path:
    candidate = Path(raw_path)
    if candidate.is_file():
        return candidate.resolve()
    matches = sorted((root / MODE_ROOTS[mode]).rglob(candidate.name), key=str)
    if len(matches) == 1:
        return matches[0].resolve()
    if not matches:
        raise StatisticsError(f"Generation FASTA does not exist: {raw_path}")
    raise StatisticsError(
        f"Generation FASTA basename is ambiguous: {candidate.name}"
    )


def count_fasta(path: Path) -> int:
    return sum(1 for _ in pipeline.iter_fasta(path))


def analyze_target(
    root: Path,
    mode: str,
    row: Mapping[str, str],
    parameters: Parameters,
) -> Dict[str, object]:
    input_path = resolve_generation_fasta(root, mode, row["output_fasta"])
    all_records = list(pipeline.iter_fasta(input_path))
    atcg_records = [
        record for record in all_records if pipeline.DNA_RE.fullmatch(record[1])
    ]

    stage_one: List[Tuple[str, str, float, float]] = []
    removed_gc = 0
    removed_complexity = 0
    gc_fail_marginal = 0
    complexity_fail_marginal = 0
    failed_both = 0
    for header, sequence in atcg_records:
        gc = pipeline.gc_percent(sequence)
        complexity = pipeline.dust_complexity(sequence)
        gc_pass = parameters.gc_min <= gc <= parameters.gc_max
        complexity_pass = (
            parameters.complexity_min <= complexity <= parameters.complexity_max
        )
        if not gc_pass:
            removed_gc += 1
            gc_fail_marginal += 1
        if not complexity_pass:
            complexity_fail_marginal += 1
            if gc_pass:
                removed_complexity += 1
        if not gc_pass and not complexity_pass:
            failed_both += 1
        if gc_pass and complexity_pass:
            stage_one.append((header, sequence, gc, complexity))

    if len(stage_one) > 1 and parameters.dimer_threshold > 0:
        scores = pipeline.dimer_scores(
            [(header, sequence) for header, sequence, _, _ in stage_one],
            parameters.dimer_k,
        )
        cutoff = pipeline.eprobe_dimer_cutoff(scores, parameters.dimer_threshold)
    else:
        scores = [0.0] * len(stage_one)
        cutoff = math.inf
    after_dimer = [
        (*record, score)
        for record, score in zip(stage_one, scores)
        if score <= cutoff
    ]
    after_dimer.sort(key=lambda item: (item[4], item[3], item[2]))
    deduplicated = []
    seen = set()
    for record in after_dimer:
        if record[1] in seen:
            continue
        seen.add(record[1])
        deduplicated.append(record)

    assessed_name = input_path.name.replace(".fasta", ".assessed.fasta")
    assessed_path = root / "06_optional_modes" / mode / "access" / assessed_name
    if not assessed_path.is_file():
        raise StatisticsError(f"Assessed FASTA does not exist: {assessed_path}")
    actual_assessed = count_fasta(assessed_path)
    after_atcg = len(atcg_records)
    after_gc = after_atcg - removed_gc
    return {
        "target_id": row["target_id"],
        "rank": row["rank"],
        "name": row["name"],
        "taxid": row["taxid"],
        "input_records": len(all_records),
        "removed_non_atcg": len(all_records) - after_atcg,
        "after_atcg": after_atcg,
        "removed_by_gc": removed_gc,
        "after_gc": after_gc,
        "removed_by_complexity_after_gc": removed_complexity,
        "after_gc_and_complexity": len(stage_one),
        "removed_by_dimer": len(stage_one) - len(after_dimer),
        "after_dimer": len(after_dimer),
        "removed_by_exact_sequence_dedup": len(after_dimer) - len(deduplicated),
        "recomputed_assessed": len(deduplicated),
        "actual_assessed": actual_assessed,
        "recomputed_matches_actual": (
            "YES" if len(deduplicated) == actual_assessed else "NO"
        ),
        "gc_fail_regardless_of_complexity": gc_fail_marginal,
        "complexity_fail_regardless_of_gc": complexity_fail_marginal,
        "failed_both_gc_and_complexity": failed_both,
        "input_fasta": str(input_path),
        "assessed_fasta": str(assessed_path),
    }


def total(rows: Sequence[Mapping[str, object]], key: str) -> int:
    return sum(int(row[key]) for row in rows)


def write_per_target(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_summary(
    path: Path,
    mode: str,
    rows: Sequence[Mapping[str, object]],
    parameters: Parameters,
) -> None:
    metrics = [
        ("mode", mode, ""),
        ("target_files", len(rows), "Number of independently assessed targets"),
        ("gc_min_percent", parameters.gc_min, "Inclusive lower GC cutoff"),
        ("gc_max_percent", parameters.gc_max, "Inclusive upper GC cutoff"),
        ("complexity_min", parameters.complexity_min, "Inclusive DUST lower cutoff"),
        ("complexity_max", parameters.complexity_max, "Inclusive DUST upper cutoff"),
        ("dimer_k", parameters.dimer_k, "Reverse-complement k-mer length"),
        (
            "dimer_threshold",
            parameters.dimer_threshold,
            "0<x<1 quantile; x>=1 absolute cutoff; x<=0 disabled",
        ),
        ("input_records", total(rows, "input_records"), "Before eMito-access"),
        ("removed_non_atcg", total(rows, "removed_non_atcg"), "Before metric filters"),
        ("after_atcg", total(rows, "after_atcg"), ""),
        ("removed_by_gc", total(rows, "removed_by_gc"), "Sequential attribution"),
        ("after_gc", total(rows, "after_gc"), ""),
        (
            "removed_by_complexity_after_gc",
            total(rows, "removed_by_complexity_after_gc"),
            "Sequential attribution among GC-passing probes",
        ),
        ("after_gc_and_complexity", total(rows, "after_gc_and_complexity"), ""),
        ("removed_by_dimer", total(rows, "removed_by_dimer"), ""),
        ("after_dimer", total(rows, "after_dimer"), ""),
        (
            "removed_by_exact_sequence_dedup",
            total(rows, "removed_by_exact_sequence_dedup"),
            "Within each target",
        ),
        ("recomputed_assessed", total(rows, "recomputed_assessed"), ""),
        ("actual_assessed", total(rows, "actual_assessed"), "Existing assessed FASTAs"),
        (
            "targets_matching_existing_assessed_count",
            sum(row["recomputed_matches_actual"] == "YES" for row in rows),
            f"of {len(rows)} targets",
        ),
        (
            "gc_fail_regardless_of_complexity",
            total(rows, "gc_fail_regardless_of_complexity"),
            "Marginal count; may overlap complexity failures",
        ),
        (
            "complexity_fail_regardless_of_gc",
            total(rows, "complexity_fail_regardless_of_gc"),
            "Marginal count; may overlap GC failures",
        ),
        (
            "failed_both_gc_and_complexity",
            total(rows, "failed_both_gc_and_complexity"),
            "Overlap of the two marginal counts",
        ),
    ]
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value", "note"))
        writer.writerows(metrics)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        root = args.output_root.resolve()
        if not root.is_dir():
            raise StatisticsError(f"Pipeline output root does not exist: {root}")
        config = load_config(root / "00_config.json")
        if not bool(config.get(f"{args.mode}_access", False)):
            raise StatisticsError(
                f"eMito-access was not enabled for mode {args.mode!r} in this run"
            )
        parameters = load_parameters(config)
        rows = load_manifest(root / MODE_ROOTS[args.mode] / "manifest.tsv")
        results = []
        for index, row in enumerate(rows, start=1):
            results.append(analyze_target(root, args.mode, row, parameters))
            if index % 25 == 0 or index == len(rows):
                print(f"Analyzed {index:,}/{len(rows):,} {args.mode} targets")

        report_dir = (
            args.report_dir.resolve()
            if args.report_dir is not None
            else root / "08_probe_statistics" / f"{args.mode}_access_filters"
        )
        report_dir.mkdir(parents=True, exist_ok=True)
        per_target = report_dir / f"{args.mode}_access_filter_per_target.tsv"
        summary = report_dir / f"{args.mode}_access_filter_summary.tsv"
        write_per_target(per_target, results)
        write_summary(summary, args.mode, results, parameters)
        print(f"Summary report: {summary}")
        print(f"Per-target report: {per_target}")
        return 0
    except (
        StatisticsError,
        pipeline.PipelineError,
        OSError,
        UnicodeError,
        ValueError,
        json.JSONDecodeError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
