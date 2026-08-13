#!/usr/bin/env python3
"""Summarize eMito probe counts by mode and processing stage.

The report deliberately contains no per-species, per-subgroup, per-node or
per-file rows.  "Merged" means logical FASTA concatenation (record count), and
"deduplicated" means exact deduplication by uppercase ATCG sequence.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple


DEFAULT_OUTPUT_ROOT = Path("emito_output")
MODE_ORDER = ("taxa", "group", "node")
DNA_TO_BITS = {"A": 0, "C": 1, "G": 2, "T": 3}


class CountError(RuntimeError):
    pass


@dataclass(frozen=True)
class SequenceCount:
    files: int
    merged: int
    unique: int
    invalid: int


@dataclass(frozen=True)
class ModeSummary:
    mode: str
    generation: SequenceCount
    access_enabled: bool
    access: Optional[SequenceCount]
    collapse_enabled: bool
    collapse: Optional[SequenceCount]
    terminal_stage: str
    terminal: SequenceCount


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Show generation/access/collapse and final merged/deduplicated probe "
            "counts by mode, without per-taxon details."
        ),
    )
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--report",
        type=Path,
        help=(
            "Mode TSV output; default is "
            "<output-root>/08_probe_statistics/mode_processing_summary.tsv."
        ),
    )
    parser.add_argument(
        "--final-report",
        type=Path,
        help=(
            "Final TSV output; default is "
            "<output-root>/08_probe_statistics/final_merge_summary.tsv."
        ),
    )
    return parser.parse_args(argv)


def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header: Optional[str] = None
    parts: List[str] = []
    with path.open("rt", encoding="ascii", errors="strict") as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.strip()
            if not text:
                continue
            if text.startswith(">"):
                if header is not None:
                    if not parts:
                        raise CountError(f"Empty FASTA record in {path}")
                    yield header, "".join(parts).upper()
                header = text[1:].strip()
                parts = []
            else:
                if header is None:
                    raise CountError(
                        f"Sequence before first FASTA header: {path}:{line_number}"
                    )
                parts.append(text)
    if header is not None:
        if not parts:
            raise CountError(f"Empty FASTA record in {path}")
        yield header, "".join(parts).upper()


def encode_atcg(sequence: str) -> Optional[int]:
    """Encode an exact sequence compactly, retaining sequence-length identity."""
    encoded = 1
    for base in sequence:
        value = DNA_TO_BITS.get(base)
        if value is None:
            return None
        encoded = (encoded << 2) | value
    return encoded


def fasta_files(root: Path) -> Tuple[Path, ...]:
    if not root.is_dir():
        return ()
    return tuple(
        sorted(
            (path.resolve() for path in root.rglob("*.fasta") if path.is_file()),
            key=str,
        )
    )


def count_paths(paths: Iterable[Path]) -> SequenceCount:
    path_list = tuple(sorted(set(paths), key=str))
    merged = 0
    invalid = 0
    unique: Set[int] = set()
    for path in path_list:
        for _, sequence in iter_fasta(path):
            merged += 1
            encoded = encode_atcg(sequence)
            if encoded is None:
                invalid += 1
            else:
                unique.add(encoded)
    return SequenceCount(len(path_list), merged, len(unique), invalid)


def load_config(output_root: Path) -> Mapping[str, object]:
    path = output_root / "00_config.json"
    if not path.is_file():
        raise CountError(f"Pipeline configuration does not exist: {path}")
    with path.open("rt", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise CountError(f"Pipeline configuration is not a JSON object: {path}")
    return value


def config_boolean(config: Mapping[str, object], key: str) -> bool:
    if key not in config:
        raise CountError(f"00_config.json lacks required option: {key}")
    value = config[key]
    if not isinstance(value, bool):
        raise CountError(f"00_config.json option {key} is not boolean: {value!r}")
    return value


def generation_files(output_root: Path) -> Dict[str, Tuple[Path, ...]]:
    return {
        "taxa": fasta_files(output_root / "03_eMito_taxa_generate" / "probe_sets"),
        "group": fasta_files(output_root / "04_eMito_group_generate" / "probe_sets"),
        "node": fasta_files(output_root / "05_eMito_node_generate" / "probe_sets"),
    }


def resolve_manifest_path(raw: str, output_root: Path) -> Path:
    path = Path(raw)
    if path.is_file():
        return path.resolve()
    matches = sorted(
        (match.resolve() for match in output_root.rglob(path.name)), key=str
    )
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise CountError(f"Final-input manifest FASTA does not exist: {path}")
    raise CountError(
        f"Final-input manifest basename is ambiguous after path remapping: {path.name}"
    )


def terminal_files(output_root: Path) -> Dict[str, Tuple[Path, ...]]:
    manifest = output_root / "07_final_probe_set" / "input_probe_files.tsv"
    if not manifest.is_file():
        raise CountError(f"Final-input manifest does not exist: {manifest}")
    result: Dict[str, List[Path]] = {mode: [] for mode in MODE_ORDER}
    with manifest.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"mode", "input_probe_file"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise CountError(f"Final-input manifest lacks columns {sorted(required)}")
        for row in reader:
            mode = row["mode"].strip()
            if mode not in result:
                raise CountError(f"Unknown mode in final-input manifest: {mode!r}")
            result[mode].append(
                resolve_manifest_path(row["input_probe_file"].strip(), output_root)
            )
    return {
        mode: tuple(sorted(set(paths), key=str)) for mode, paths in result.items()
    }


def expected_terminal_stage(access_enabled: bool, collapse_enabled: bool) -> str:
    if collapse_enabled:
        return "collapse"
    if access_enabled:
        return "access"
    return "generation"


def build_summaries(
    output_root: Path, config: Mapping[str, object]
) -> Tuple[List[ModeSummary], SequenceCount, SequenceCount, bool]:
    generated = generation_files(output_root)
    terminal = terminal_files(output_root)
    summaries: List[ModeSummary] = []
    all_terminal_paths: List[Path] = []
    for mode in MODE_ORDER:
        access_enabled = config_boolean(config, f"{mode}_access")
        collapse_enabled = config_boolean(config, f"{mode}_collapse")
        access_paths = fasta_files(output_root / "06_optional_modes" / mode / "access")
        collapse_paths = fasta_files(
            output_root / "06_optional_modes" / mode / "collapse"
        )
        if access_enabled and not access_paths:
            raise CountError(f"{mode} access is enabled but no access FASTA files exist")
        if collapse_enabled and not collapse_paths:
            raise CountError(f"{mode} collapse is enabled but no collapse FASTA files exist")
        expected_stage = expected_terminal_stage(access_enabled, collapse_enabled)
        expected_paths = {
            "generation": generated[mode],
            "access": access_paths,
            "collapse": collapse_paths,
        }[expected_stage]
        if set(terminal[mode]) != set(expected_paths):
            raise CountError(
                f"{mode} final-input manifest does not match configured terminal "
                f"stage {expected_stage}: manifest_files={len(terminal[mode])}, "
                f"expected_files={len(expected_paths)}"
            )
        all_terminal_paths.extend(terminal[mode])
        summaries.append(
            ModeSummary(
                mode=mode,
                generation=count_paths(generated[mode]),
                access_enabled=access_enabled,
                access=count_paths(access_paths) if access_enabled else None,
                collapse_enabled=collapse_enabled,
                collapse=count_paths(collapse_paths) if collapse_enabled else None,
                terminal_stage=expected_stage,
                terminal=count_paths(terminal[mode]),
            )
        )

    terminal_total = count_paths(all_terminal_paths)
    final_path = output_root / "07_final_probe_set" / "final_probe_set.fasta"
    if not final_path.is_file():
        raise CountError(f"Final probe FASTA does not exist: {final_path}")
    final_fasta = count_paths((final_path,))
    final_dedup = config_boolean(config, "final_dedup")
    expected_final_records = (
        terminal_total.unique if final_dedup else terminal_total.merged - terminal_total.invalid
    )
    if final_fasta.merged != expected_final_records:
        raise CountError(
            "Final FASTA record count differs from the configured final merge "
            f"behavior (final_dedup={final_dedup}): final={final_fasta.merged}, "
            f"expected={expected_final_records}"
        )
    return summaries, terminal_total, final_fasta, final_dedup


def optional_value(value: Optional[SequenceCount], field: str) -> str:
    if value is None:
        return "NA"
    return str(getattr(value, field))


def write_mode_report(path: Path, summaries: Sequence[ModeSummary]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "mode",
                "generation_input_files",
                "generation_merged_sequences",
                "generation_unique_after_dedup",
                "access_executed",
                "access_merged_sequences",
                "access_unique_after_dedup",
                "collapse_executed",
                "collapse_merged_sequences",
                "collapse_unique_after_dedup",
                "terminal_stage_used_for_final",
                "terminal_merged_sequences",
                "terminal_unique_after_dedup",
            )
        )
        for row in summaries:
            writer.writerow(
                (
                    row.mode,
                    row.generation.files,
                    row.generation.merged,
                    row.generation.unique,
                    "YES" if row.access_enabled else "NO",
                    optional_value(row.access, "merged"),
                    optional_value(row.access, "unique"),
                    "YES" if row.collapse_enabled else "NO",
                    optional_value(row.collapse, "merged"),
                    optional_value(row.collapse, "unique"),
                    row.terminal_stage,
                    row.terminal.merged,
                    row.terminal.unique,
                )
            )


def write_final_report(
    path: Path,
    terminal_total: SequenceCount,
    final_fasta: SequenceCount,
    final_dedup: bool,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "final_input_files",
                "sequences_after_merging_terminal_taxa_group_node",
                "final_dedup_executed",
                "unique_sequences_in_merged_inputs",
                "duplicate_sequences_removed_from_final_fasta",
                "records_in_final_probe_set_fasta",
            )
        )
        writer.writerow(
            (
                terminal_total.files,
                terminal_total.merged,
                "YES" if final_dedup else "NO",
                terminal_total.unique,
                (
                    terminal_total.merged
                    - terminal_total.invalid
                    - terminal_total.unique
                    if final_dedup
                    else 0
                ),
                final_fasta.merged,
            )
        )


def print_summary(
    summaries: Sequence[ModeSummary],
    terminal_total: SequenceCount,
    final_fasta: SequenceCount,
    final_dedup: bool,
) -> None:
    print(
        "mode\tgeneration_merged\tgeneration_unique\taccess\t"
        "access_merged\taccess_unique\tcollapse\tcollapse_merged\t"
        "collapse_unique\tfinal_uses"
    )
    for row in summaries:
        print(
            f"{row.mode}\t{row.generation.merged:,}\t{row.generation.unique:,}\t"
            f"{'YES' if row.access_enabled else 'NO'}\t"
            f"{optional_value(row.access, 'merged')}\t"
            f"{optional_value(row.access, 'unique')}\t"
            f"{'YES' if row.collapse_enabled else 'NO'}\t"
            f"{optional_value(row.collapse, 'merged')}\t"
            f"{optional_value(row.collapse, 'unique')}\t{row.terminal_stage}"
        )
    print("\nFinal merge:")
    print(
        "terminal_files\tsequences_after_merge\tfinal_dedup\tunique_in_inputs\t"
        "duplicates_removed\tfinal_fasta_records"
    )
    print(
        f"{terminal_total.files:,}\t{terminal_total.merged:,}\t"
        f"{'YES' if final_dedup else 'NO'}\t{terminal_total.unique:,}\t"
        f"{terminal_total.merged - terminal_total.invalid - terminal_total.unique if final_dedup else 0:,}\t"
        f"{final_fasta.merged:,}"
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        output_root = args.output_root.resolve()
        if not output_root.is_dir():
            raise CountError(f"Pipeline output root does not exist: {output_root}")
        report_dir = output_root / "08_probe_statistics"
        mode_report = (
            args.report.resolve()
            if args.report is not None
            else report_dir / "mode_processing_summary.tsv"
        )
        final_report = (
            args.final_report.resolve()
            if args.final_report is not None
            else report_dir / "final_merge_summary.tsv"
        )
        config = load_config(output_root)
        summaries, terminal_total, final_fasta, final_dedup = build_summaries(
            output_root, config
        )
        write_mode_report(mode_report, summaries)
        write_final_report(final_report, terminal_total, final_fasta, final_dedup)
        print_summary(summaries, terminal_total, final_fasta, final_dedup)
        print(f"\nMode report: {mode_report}")
        print(f"Final report: {final_report}")
        return 0
    except (CountError, OSError, UnicodeError, ValueError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
