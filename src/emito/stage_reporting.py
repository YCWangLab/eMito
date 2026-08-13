#!/usr/bin/env python3
"""Count probe records and sequence-deduplicated probes at every eMito stage.

The script is read-only with respect to an existing pipeline result. It writes
TSV reports under <output-root>/08_probe_statistics by default. Deduplication is
based only on the uppercase ATCG sequence, regardless of FASTA header, accession
or coordinate.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple


DEFAULT_OUTPUT_ROOT = Path("emito_output")
DNA_TO_BITS = {"A": 0, "C": 1, "G": 2, "T": 3}
MODE_ORDER = ("taxa", "group", "node")


class StatisticsError(RuntimeError):
    pass


@dataclass(frozen=True)
class Stage:
    stage_id: str
    description: str
    files: Tuple[Path, ...]


@dataclass(frozen=True)
class FileCount:
    records: int
    atcg_records: int
    invalid_records: int
    unique_sequences: int


@dataclass(frozen=True)
class StageCount:
    stage_id: str
    description: str
    file_count: int
    empty_file_count: int
    records: int
    atcg_records: int
    invalid_records: int
    sum_unique_per_file: int
    unique_after_stage_merge: int

    @property
    def duplicates_within_files(self) -> int:
        return self.atcg_records - self.sum_unique_per_file

    @property
    def duplicates_across_files(self) -> int:
        return self.sum_unique_per_file - self.unique_after_stage_merge

    @property
    def duplicates_removed_by_stage_merge(self) -> int:
        return self.atcg_records - self.unique_after_stage_merge

    @property
    def retained_percent(self) -> float:
        if self.atcg_records == 0:
            return 0.0
        return 100.0 * self.unique_after_stage_merge / self.atcg_records


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Count FASTA probe records and sequence-deduplicated probes for "
            "every stage of an eMito pipeline result."
        ),
    )
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--report-dir",
        type=Path,
        help="Report directory; default is <output-root>/08_probe_statistics.",
    )
    parser.add_argument(
        "--skip-kmers",
        action="store_true",
        help="Count only probe stages and skip the larger k-mer reports.",
    )
    return parser.parse_args(argv)


def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header: Optional[str] = None
    sequence_parts: List[str] = []
    with path.open("rt", encoding="ascii", errors="strict") as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.strip()
            if not text:
                continue
            if text.startswith(">"):
                if header is not None:
                    if not sequence_parts:
                        raise StatisticsError(f"Empty FASTA record in {path}")
                    yield header, "".join(sequence_parts).upper()
                header = text[1:].strip()
                sequence_parts = []
            else:
                if header is None:
                    raise StatisticsError(
                        f"Sequence occurs before the first header: {path}:{line_number}"
                    )
                sequence_parts.append(text)
    if header is not None:
        if not sequence_parts:
            raise StatisticsError(f"Empty FASTA record in {path}")
        yield header, "".join(sequence_parts).upper()


def encode_atcg(sequence: str) -> Optional[int]:
    """Encode ATCG plus a leading sentinel as one compact arbitrary-size int."""
    encoded = 1
    for base in sequence:
        value = DNA_TO_BITS.get(base)
        if value is None:
            return None
        encoded = (encoded << 2) | value
    return encoded


def count_file(path: Path, stage_sequences: Set[int]) -> FileCount:
    records = 0
    atcg_records = 0
    invalid_records = 0
    file_sequences: Set[int] = set()
    for _, sequence in iter_fasta(path):
        records += 1
        encoded = encode_atcg(sequence)
        if encoded is None:
            invalid_records += 1
            continue
        atcg_records += 1
        file_sequences.add(encoded)
        stage_sequences.add(encoded)
    return FileCount(
        records=records,
        atcg_records=atcg_records,
        invalid_records=invalid_records,
        unique_sequences=len(file_sequences),
    )


def count_stage(
    stage: Stage,
    output_root: Path,
) -> Tuple[StageCount, List[Tuple[str, str, FileCount]]]:
    stage_sequences: Set[int] = set()
    records = 0
    atcg_records = 0
    invalid_records = 0
    sum_unique_per_file = 0
    empty_file_count = 0
    file_rows: List[Tuple[str, str, FileCount]] = []

    for index, path in enumerate(stage.files, start=1):
        result = count_file(path, stage_sequences)
        records += result.records
        atcg_records += result.atcg_records
        invalid_records += result.invalid_records
        sum_unique_per_file += result.unique_sequences
        if result.records == 0:
            empty_file_count += 1
        try:
            display_path = str(path.relative_to(output_root))
        except ValueError:
            display_path = str(path)
        file_rows.append((stage.stage_id, display_path, result))
        if index % 250 == 0 or index == len(stage.files):
            log(f"{stage.stage_id}: counted {index:,}/{len(stage.files):,} FASTA files")

    return (
        StageCount(
            stage_id=stage.stage_id,
            description=stage.description,
            file_count=len(stage.files),
            empty_file_count=empty_file_count,
            records=records,
            atcg_records=atcg_records,
            invalid_records=invalid_records,
            sum_unique_per_file=sum_unique_per_file,
            unique_after_stage_merge=len(stage_sequences),
        ),
        file_rows,
    )


def fasta_files(root: Path, pattern: str = "*.fasta") -> Tuple[Path, ...]:
    if not root.is_dir():
        return ()
    return tuple(sorted((path for path in root.rglob(pattern) if path.is_file()), key=str))


def unique_paths(paths: Iterable[Path]) -> Tuple[Path, ...]:
    return tuple(sorted(set(paths), key=str))


def resolve_manifest_path(raw_path: str, output_root: Path) -> Path:
    path = Path(raw_path)
    if path.is_file():
        return path
    candidates = list(output_root.rglob(path.name))
    if len(candidates) == 1:
        log(f"WARNING: remapped missing manifest path {path} to {candidates[0]}")
        return candidates[0]
    if not candidates:
        raise StatisticsError(f"Manifest input FASTA does not exist: {path}")
    raise StatisticsError(
        f"Manifest input FASTA is missing and basename is ambiguous: {path}; "
        f"matches={len(candidates)}"
    )


def load_terminal_files(output_root: Path) -> Dict[str, Tuple[Path, ...]]:
    manifest = output_root / "07_final_probe_set" / "input_probe_files.tsv"
    result: Dict[str, List[Path]] = {mode: [] for mode in MODE_ORDER}
    if not manifest.is_file():
        log(f"WARNING: terminal-input manifest is absent: {manifest}")
        return {mode: () for mode in MODE_ORDER}
    with manifest.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"mode", "input_probe_file"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise StatisticsError(
                f"Terminal manifest lacks columns {sorted(required)}: {manifest}"
            )
        for row in reader:
            mode = row["mode"].strip()
            if mode not in result:
                raise StatisticsError(f"Unknown mode {mode!r} in {manifest}")
            result[mode].append(resolve_manifest_path(row["input_probe_file"], output_root))
    return {mode: unique_paths(paths) for mode, paths in result.items()}


def build_stages(output_root: Path) -> Tuple[List[Stage], Mapping[str, Tuple[Path, ...]]]:
    taxa_sets = fasta_files(output_root / "03_eMito_taxa_generate" / "probe_sets")
    group_sets = fasta_files(output_root / "04_eMito_group_generate" / "probe_sets")
    node_sets = fasta_files(output_root / "05_eMito_node_generate" / "probe_sets")
    terminal = load_terminal_files(output_root)

    stages = [
        Stage(
            "01_initial_probe_windows",
            "All step-sized probe windows generated from individual genomes",
            fasta_files(output_root / "02_windows", "*.probe.fasta"),
        ),
        Stage(
            "02_taxa_individual_intersections",
            "Per-accession intersections with combined species/genus-specific k-mers",
            fasta_files(output_root / "03_eMito_taxa_generate" / "intersections"),
        ),
        Stage(
            "03_taxa_specific_probe_sets",
            "One merged and deduplicated Taxa-specific probe set per species",
            taxa_sets,
        ),
        Stage(
            "04_group_individual_intersections",
            "Per-accession intersections with subgroup-specific k-mers",
            fasta_files(output_root / "04_eMito_group_generate" / "intersections"),
        ),
        Stage(
            "05_group_specific_probe_sets",
            "One merged and deduplicated Group-specific probe set per subgroup",
            group_sets,
        ),
        Stage(
            "06_node_tiling_probe_sets",
            "Node-tiling probe sets before optional processing",
            node_sets,
        ),
        Stage(
            "07_generation_outputs_all_modes",
            "Taxa, group and node generation outputs combined",
            unique_paths((*taxa_sets, *group_sets, *node_sets)),
        ),
    ]

    for mode in MODE_ORDER:
        stages.append(
            Stage(
                f"08_{mode}_access_outputs",
                f"All {mode} eMito-access assessed FASTA outputs",
                fasta_files(output_root / "06_optional_modes" / mode / "access"),
            )
        )
        stages.append(
            Stage(
                f"09_{mode}_collapse_outputs",
                f"All {mode} eMito-collapse FASTA outputs",
                fasta_files(output_root / "06_optional_modes" / mode / "collapse"),
            )
        )
        stages.append(
            Stage(
                f"10_{mode}_terminal_outputs",
                f"{mode} files actually selected for the final merge",
                terminal[mode],
            )
        )

    terminal_all = unique_paths(path for mode in MODE_ORDER for path in terminal[mode])
    final_path = output_root / "07_final_probe_set" / "final_probe_set.fasta"
    stages.extend(
        [
            Stage(
                "11_terminal_outputs_all_modes",
                "All terminal mode outputs selected for final merge",
                terminal_all,
            ),
            Stage(
                "12_final_probe_set",
                "Pipeline final sequence-deduplicated probe set",
                (final_path,) if final_path.is_file() else (),
            ),
        ]
    )
    return stages, terminal


def build_kmer_stages(output_root: Path) -> List[Stage]:
    return [
        Stage(
            "K01_individual_forward_kmers",
            "Forward k-mer windows from every individual genome",
            fasta_files(output_root / "02_windows", "*.kmer.forward.fasta"),
        ),
        Stage(
            "K02_individual_reverse_complement_kmers",
            "Reverse-complement k-mer windows from every individual genome",
            fasta_files(
                output_root / "02_windows", "*.kmer.reverse_complement.fasta"
            ),
        ),
        Stage(
            "K03_individual_merged_kmers",
            "Forward plus reverse-complement k-mers for every individual",
            fasta_files(output_root / "02_windows", "*.kmer.merged.fasta"),
        ),
        Stage(
            "K04_taxa_species_representative_kmers",
            "Representative k-mers for each eligible species",
            fasta_files(output_root / "03_eMito_taxa_generate" / "representative_kmers"),
        ),
        Stage(
            "K05_taxa_species_specific_kmers",
            "Species-specific k-mers before combination with genus-specific k-mers",
            fasta_files(
                output_root / "03_eMito_taxa_generate" / "specific_kmers" / "species"
            ),
        ),
        Stage(
            "K06_taxa_genus_specific_kmers",
            "Genus-specific k-mers retained as filtering intermediates",
            fasta_files(
                output_root / "03_eMito_taxa_generate" / "specific_kmers" / "genus"
            ),
        ),
        Stage(
            "K06b_taxa_genus_specific_kmers_by_species",
            "Genus-specific k-mers restricted to each species' representative set",
            fasta_files(
                output_root
                / "03_eMito_taxa_generate"
                / "specific_kmers"
                / "genus_by_species"
            ),
        ),
        Stage(
            "K07_taxa_combined_kmers_by_species",
            "Species-specific plus its genus-specific k-mers, one file per species",
            fasta_files(
                output_root
                / "03_eMito_taxa_generate"
                / "specific_kmers"
                / "combined_by_species"
            ),
        ),
        Stage(
            "K08_group_representative_kmers",
            "Representative k-mers for each subgroup",
            fasta_files(output_root / "04_eMito_group_generate" / "representative_kmers"),
        ),
        Stage(
            "K09_group_specific_kmers",
            "Subgroup-specific k-mers",
            fasta_files(output_root / "04_eMito_group_generate" / "specific_kmers"),
        ),
    ]


def read_sequence_set(paths: Iterable[Path]) -> Set[int]:
    result: Set[int] = set()
    for path in paths:
        for _, sequence in iter_fasta(path):
            encoded = encode_atcg(sequence)
            if encoded is not None:
                result.add(encoded)
    return result


def write_terminal_contributions(
    path: Path,
    terminal: Mapping[str, Tuple[Path, ...]],
) -> None:
    mode_sets = {mode: read_sequence_set(terminal[mode]) for mode in MODE_ORDER}
    union: Set[int] = set()
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "mode",
                "unique_probe_sequences",
                "new_sequences_added_in_taxa_group_node_order",
                "sequences_already_present_in_earlier_modes",
                "cumulative_unique_sequences",
            )
        )
        for mode in MODE_ORDER:
            sequences = mode_sets[mode]
            new_sequences = sequences - union
            already_present = len(sequences) - len(new_sequences)
            union.update(sequences)
            writer.writerow(
                (mode, len(sequences), len(new_sequences), already_present, len(union))
            )

        writer.writerow(())
        writer.writerow(("overlap", "unique_probe_sequences"))
        writer.writerow(("taxa_and_group", len(mode_sets["taxa"] & mode_sets["group"])))
        writer.writerow(("taxa_and_node", len(mode_sets["taxa"] & mode_sets["node"])))
        writer.writerow(("group_and_node", len(mode_sets["group"] & mode_sets["node"])))
        writer.writerow(
            (
                "taxa_and_group_and_node",
                len(mode_sets["taxa"] & mode_sets["group"] & mode_sets["node"]),
            )
        )


def write_summary(path: Path, counts: Sequence[StageCount]) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "stage",
                "description",
                "fasta_files",
                "empty_fasta_files",
                "fasta_records",
                "atcg_records",
                "invalid_records",
                "sum_unique_within_each_file",
                "duplicates_within_files",
                "unique_after_merging_all_stage_files",
                "duplicates_across_files",
                "total_duplicates_removed_by_stage_merge",
                "unique_retained_percent",
            )
        )
        for result in counts:
            writer.writerow(
                (
                    result.stage_id,
                    result.description,
                    result.file_count,
                    result.empty_file_count,
                    result.records,
                    result.atcg_records,
                    result.invalid_records,
                    result.sum_unique_per_file,
                    result.duplicates_within_files,
                    result.unique_after_stage_merge,
                    result.duplicates_across_files,
                    result.duplicates_removed_by_stage_merge,
                    f"{result.retained_percent:.4f}",
                )
            )


def write_file_details(
    path: Path,
    rows: Sequence[Tuple[str, str, FileCount]],
) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "stage",
                "fasta_file",
                "fasta_records",
                "atcg_records",
                "invalid_records",
                "unique_sequences_within_file",
                "duplicate_records_within_file",
            )
        )
        for stage_id, display_path, result in rows:
            writer.writerow(
                (
                    stage_id,
                    display_path,
                    result.records,
                    result.atcg_records,
                    result.invalid_records,
                    result.unique_sequences,
                    result.atcg_records - result.unique_sequences,
                )
            )


def print_summary(counts: Sequence[StageCount]) -> None:
    print(
        "stage\tfiles\trecords\tsum_unique_per_file\t"
        "unique_after_stage_merge\tduplicates_removed"
    )
    for result in counts:
        print(
            f"{result.stage_id}\t{result.file_count:,}\t{result.atcg_records:,}\t"
            f"{result.sum_unique_per_file:,}\t{result.unique_after_stage_merge:,}\t"
            f"{result.duplicates_removed_by_stage_merge:,}"
        )


def validate_taxa_outputs(output_root: Path) -> None:
    probe_root = output_root / "03_eMito_taxa_generate" / "probe_sets"
    unexpected = [
        path
        for path in fasta_files(probe_root)
        if "__genus_taxid_" in path.name or "genus" in path.relative_to(probe_root).parts
    ]
    if unexpected:
        log(
            "WARNING: found genus-level taxa probe files. This result may contain "
            "outputs from the old pipeline version:"
        )
        for path in unexpected[:20]:
            log(f"  {path}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        output_root = args.output_root.resolve()
        if not output_root.is_dir():
            raise StatisticsError(f"Pipeline output root does not exist: {output_root}")
        report_dir = (
            args.report_dir.resolve()
            if args.report_dir is not None
            else output_root / "08_probe_statistics"
        )
        report_dir.mkdir(parents=True, exist_ok=True)

        validate_taxa_outputs(output_root)
        probe_stages, terminal = build_stages(output_root)
        probe_counts: List[StageCount] = []
        probe_file_rows: List[Tuple[str, str, FileCount]] = []
        for stage in probe_stages:
            log(f"Counting {stage.stage_id}: {stage.description}")
            result, details = count_stage(stage, output_root)
            probe_counts.append(result)
            probe_file_rows.extend(details)

        kmer_counts: List[StageCount] = []
        kmer_file_rows: List[Tuple[str, str, FileCount]] = []
        if not args.skip_kmers:
            for stage in build_kmer_stages(output_root):
                log(f"Counting {stage.stage_id}: {stage.description}")
                result, details = count_stage(stage, output_root)
                kmer_counts.append(result)
                kmer_file_rows.extend(details)

        summary_path = report_dir / "probe_stage_summary.tsv"
        details_path = report_dir / "probe_file_counts.tsv"
        kmer_summary_path = report_dir / "kmer_stage_summary.tsv"
        kmer_details_path = report_dir / "kmer_file_counts.tsv"
        contribution_path = report_dir / "probe_terminal_mode_contributions.tsv"
        write_summary(summary_path, probe_counts)
        write_file_details(details_path, probe_file_rows)
        if not args.skip_kmers:
            write_summary(kmer_summary_path, kmer_counts)
            write_file_details(kmer_details_path, kmer_file_rows)
        write_terminal_contributions(contribution_path, terminal)

        by_stage = {result.stage_id: result for result in probe_counts}
        terminal_count = by_stage["11_terminal_outputs_all_modes"].unique_after_stage_merge
        final_count = by_stage["12_final_probe_set"].unique_after_stage_merge
        if terminal_count != final_count:
            log(
                "WARNING: terminal merged unique count does not equal final FASTA count: "
                f"{terminal_count:,} vs {final_count:,}"
            )

        print("Probe stages:")
        print_summary(probe_counts)
        if not args.skip_kmers:
            print("\nK-mer stages:")
            print_summary(kmer_counts)
        print(f"\nSummary report: {summary_path}")
        print(f"Per-file report: {details_path}")
        if not args.skip_kmers:
            print(f"K-mer summary report: {kmer_summary_path}")
            print(f"K-mer per-file report: {kmer_details_path}")
        print(f"Terminal-mode contribution report: {contribution_path}")
        return 0
    except (StatisticsError, OSError, UnicodeError, ValueError) as exc:
        log(f"ERROR: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
