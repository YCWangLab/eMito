#!/usr/bin/env python3
"""Export compact publication source tables from one completed eMito run.

This script is read-only with respect to the pipeline and plotting directories.
It recomputes core aggregate reports, replays taxa eMito-access to attribute
filter losses, and copies the three TSVs underlying the manuscript probe plots.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Mapping, Optional, Sequence


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from emito import access_reporting  # noqa: E402
from emito import reporting  # noqa: E402
from emito import stage_reporting  # noqa: E402
from emito import taxonomy_reporting  # noqa: E402


PLOT_TABLES = {
    "panel_A_probe_counts.tsv": "eMito_v3_2_panel_A_probe_counts.tsv",
    "panel_B_focal_probe_counts.tsv": "eMito_v3_2_panel_B_focal_probe_counts.tsv",
    "genus_probe_counts.tsv": "genus_probe_counts.v3_2.tsv",
}

PARAMETER_KEYS = (
    "pipeline_version",
    "modes",
    "window_length",
    "probe_step",
    "kmer_step",
    "min_genome_length",
    "representative_fraction",
    "small_group_all_max",
    "node_rank",
    "node_step",
    "random_seed",
    "taxa_access",
    "taxa_collapse",
    "group_access",
    "group_collapse",
    "node_access",
    "node_collapse",
    "gc_min",
    "gc_max",
    "complexity_min",
    "complexity_max",
    "dimer_k",
    "dimer_threshold",
    "final_dedup",
)


class ExportError(RuntimeError):
    pass


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Export eMito manuscript metadata and figure source tables.",
    )
    parser.add_argument("--pipeline-root", type=Path, required=True)
    parser.add_argument(
        "--plot-dir",
        type=Path,
        required=True,
        help="Directory containing the TSV files written by the two plotting scripts.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def require_success(exit_code: int, label: str) -> None:
    if exit_code != 0:
        raise ExportError(f"{label} failed with exit code {exit_code}")


def load_config(path: Path) -> Mapping[str, object]:
    if not path.is_file():
        raise ExportError(f"Pipeline configuration does not exist: {path}")
    with path.open("rt", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ExportError(f"Pipeline configuration is not a JSON object: {path}")
    return value


def write_parameter_table(path: Path, config: Mapping[str, object]) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("parameter", "value"))
        for key in PARAMETER_KEYS:
            if key not in config:
                continue
            value = config[key]
            if isinstance(value, bool):
                value = "true" if value else "false"
            writer.writerow((key, value))


def make_pipeline_paths_portable(
    source: Path,
    destination: Path,
    pipeline_root: Path,
) -> None:
    """Copy a TSV while making pipeline-output paths root-relative."""
    with source.open("rt", encoding="utf-8", newline="") as input_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ExportError(f"Plot source table has no header: {source}")
        path_columns = [
            name
            for name in ("Output_fasta", "Terminal_fasta")
            if name in reader.fieldnames
        ]
        with destination.open("wt", encoding="utf-8", newline="") as output_handle:
            writer = csv.DictWriter(
                output_handle,
                fieldnames=reader.fieldnames,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for row in reader:
                for column in path_columns:
                    value = row.get(column, "")
                    if not value:
                        continue
                    candidate = Path(value)
                    if not candidate.is_absolute():
                        continue
                    try:
                        row[column] = candidate.relative_to(pipeline_root).as_posix()
                    except ValueError as exc:
                        raise ExportError(
                            f"{column} path is outside --pipeline-root: {candidate}"
                        ) from exc
                writer.writerow(row)


def copy_plot_tables(
    plot_dir: Path,
    output_dir: Path,
    pipeline_root: Path,
) -> None:
    missing = [
        plot_dir / source_name
        for source_name in PLOT_TABLES.values()
        if not (plot_dir / source_name).is_file()
    ]
    if missing:
        formatted = "\n".join(f"  {path}" for path in missing)
        raise ExportError(
            "Required plot source table(s) are missing. Run the final A/B and "
            f"genus-tree plotting scripts first:\n{formatted}"
        )
    for destination_name, source_name in PLOT_TABLES.items():
        make_pipeline_paths_portable(
            plot_dir / source_name,
            output_dir / destination_name,
            pipeline_root,
        )


def write_manifest(path: Path) -> None:
    descriptions = {
        "pipeline_parameters.tsv": "Resolved parameters that materially affect probe generation and filtering",
        "mode_routing.tsv": "Per-species taxa/group input decisions and accession membership after length QC",
        "input_taxonomy_summary.tsv": "Input and retained genome/species/genus/family counts",
        "mode_processing_summary.tsv": "Generation, access, collapse, and terminal counts by mode",
        "final_merge_summary.tsv": "Final cross-mode merge and sequence-deduplication counts",
        "probe_stage_summary.tsv": "Aggregate FASTA counts at all probe-processing stages",
        "probe_terminal_mode_contributions.tsv": "Unique and overlapping terminal contributions from taxa/group/node",
        "taxa_access_filter_summary.tsv": "Sequential GC, complexity, dimer, and exact-dedup counts for taxa access",
        "panel_A_probe_counts.tsv": "Absolute counts and plotted percentages underlying probe-filtering panel A",
        "panel_B_focal_probe_counts.tsv": "Focal taxon counts underlying panel B",
        "genus_probe_counts.tsv": "Final probe counts grouped by genus for the circular tree",
    }
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("file", "description"))
        for name, description in descriptions.items():
            writer.writerow((name, description))


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        pipeline_root = args.pipeline_root.resolve()
        plot_dir = args.plot_dir.resolve()
        output_dir = args.output_dir.resolve()
        if not pipeline_root.is_dir():
            raise ExportError(f"Pipeline root does not exist: {pipeline_root}")
        if not plot_dir.is_dir():
            raise ExportError(f"Plot directory does not exist: {plot_dir}")
        output_dir.mkdir(parents=True, exist_ok=True)

        config = load_config(pipeline_root / "00_config.json")
        write_parameter_table(output_dir / "pipeline_parameters.tsv", config)
        routing_source = pipeline_root / "00_input_qc" / "mode_routing.tsv"
        if not routing_source.is_file():
            raise ExportError(
                "Mode-routing audit does not exist; rerun this dataset with "
                f"eMito 0.1.4 or later: {routing_source}"
            )
        shutil.copy2(routing_source, output_dir / "mode_routing.tsv")

        require_success(
            taxonomy_reporting.main(
                [
                    "--output-root",
                    str(pipeline_root),
                    "--report",
                    str(output_dir / "input_taxonomy_summary.tsv"),
                ]
            ),
            "taxonomy summary",
        )
        require_success(
            reporting.main(
                [
                    "--output-root",
                    str(pipeline_root),
                    "--report",
                    str(output_dir / "mode_processing_summary.tsv"),
                    "--final-report",
                    str(output_dir / "final_merge_summary.tsv"),
                ]
            ),
            "mode/final summary",
        )

        with tempfile.TemporaryDirectory(prefix="emito_publication_") as temporary:
            temporary_path = Path(temporary)
            stage_dir = temporary_path / "stages"
            require_success(
                stage_reporting.main(
                    [
                        "--output-root",
                        str(pipeline_root),
                        "--report-dir",
                        str(stage_dir),
                        "--skip-kmers",
                    ]
                ),
                "probe stage summary",
            )
            shutil.copy2(
                stage_dir / "probe_stage_summary.tsv",
                output_dir / "probe_stage_summary.tsv",
            )
            shutil.copy2(
                stage_dir / "probe_terminal_mode_contributions.tsv",
                output_dir / "probe_terminal_mode_contributions.tsv",
            )

            access_dir = temporary_path / "access"
            require_success(
                access_reporting.main(
                    [
                        "--output-root",
                        str(pipeline_root),
                        "--mode",
                        "taxa",
                        "--report-dir",
                        str(access_dir),
                    ]
                ),
                "taxa access summary",
            )
            shutil.copy2(
                access_dir / "taxa_access_filter_summary.tsv",
                output_dir / "taxa_access_filter_summary.tsv",
            )

        copy_plot_tables(plot_dir, output_dir, pipeline_root)
        write_manifest(output_dir / "EXPORT_MANIFEST.tsv")
        print(f"Publication metadata export complete: {output_dir}")
        for path in sorted(output_dir.glob("*.tsv"), key=lambda item: item.name):
            print(f"  {path.name}")
        return 0
    except (
        ExportError,
        OSError,
        UnicodeError,
        ValueError,
        json.JSONDecodeError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
