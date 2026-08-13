#!/usr/bin/env python3
"""Count genomes and distinct species/genus/family nodes entering eMito."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Dict, Optional, Sequence, Set

from . import pipeline


DEFAULT_OUTPUT_ROOT = Path("emito_output")


class CountError(RuntimeError):
    pass


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Count distinct species, genera and families among genomes retained "
            "after input length QC."
        ),
    )
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--nodes-dmp",
        type=Path,
        help="Override nodes.dmp; otherwise read its path from 00_config.json.",
    )
    parser.add_argument(
        "--report",
        type=Path,
        help=(
            "Output TSV; default is "
            "<output-root>/08_probe_statistics/input_taxonomy_summary.tsv."
        ),
    )
    return parser.parse_args(argv)


def read_config(output_root: Path) -> Dict[str, object]:
    path = output_root / "00_config.json"
    if not path.is_file():
        raise CountError(f"Pipeline configuration does not exist: {path}")
    with path.open("rt", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise CountError(f"Configuration is not a JSON object: {path}")
    return value


def read_retained_manifest(path: Path) -> tuple[int, Set[int], Set[int], Set[int]]:
    if not path.is_file():
        raise CountError(f"Retained-genome manifest does not exist: {path}")
    species: Set[int] = set()
    genera: Set[int] = set()
    actual_taxids: Set[int] = set()
    genomes = 0
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"accession_id", "species_taxid", "genus_taxid", "actual_taxid"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise CountError(f"Manifest lacks columns {sorted(required)}: {path}")
        seen_accessions: Set[str] = set()
        for row in reader:
            accession = row["accession_id"].strip()
            if accession in seen_accessions:
                raise CountError(f"Duplicate accession in retained manifest: {accession}")
            seen_accessions.add(accession)
            genomes += 1
            species.add(int(row["species_taxid"]))
            genera.add(int(row["genus_taxid"]))
            actual_taxids.add(int(row["actual_taxid"]))
    return genomes, species, genera, actual_taxids


def count_qc_rows(path: Path) -> tuple[int, int, Set[int], Set[int]]:
    if not path.is_file():
        return 0, 0, set(), set()
    total = 0
    excluded = 0
    species: Set[int] = set()
    actual_taxids: Set[int] = set()
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"status", "species_taxid", "actual_taxid"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise CountError(f"Input-QC table lacks columns {sorted(required)}: {path}")
        for row in reader:
            total += 1
            species.add(int(row["species_taxid"]))
            actual_taxids.add(int(row["actual_taxid"]))
            if row["status"].strip().lower() == "excluded":
                excluded += 1
    return total, excluded, species, actual_taxids


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        output_root = args.output_root.resolve()
        if not output_root.is_dir():
            raise CountError(f"Pipeline output root does not exist: {output_root}")
        config = read_config(output_root)
        if args.nodes_dmp is not None:
            nodes_path = args.nodes_dmp
        else:
            raw_nodes = config.get("nodes_dmp")
            if not isinstance(raw_nodes, str) or not raw_nodes:
                raise CountError("00_config.json lacks a usable nodes_dmp path")
            nodes_path = Path(raw_nodes)

        genomes, species, genera, actual_taxids = read_retained_manifest(
            output_root / "01_organized_genomes" / "manifest.tsv"
        )
        nodes = pipeline.load_nodes(nodes_path)
        families = {
            pipeline.ancestor_at_rank(taxid, "family", nodes)
            for taxid in actual_taxids
        }
        qc_total, qc_excluded, raw_species, raw_actual_taxids = count_qc_rows(
            output_root / "00_input_qc" / "genome_length_qc.tsv"
        )
        if not raw_species:
            raw_species = set(species)
            raw_actual_taxids = set(actual_taxids)
        raw_genera = {
            pipeline.ancestor_at_rank(taxid, "genus", nodes)
            for taxid in raw_actual_taxids
        }
        raw_families = {
            pipeline.ancestor_at_rank(taxid, "family", nodes)
            for taxid in raw_actual_taxids
        }

        report = (
            args.report.resolve()
            if args.report is not None
            else output_root / "08_probe_statistics" / "input_taxonomy_summary.tsv"
        )
        report.parent.mkdir(parents=True, exist_ok=True)
        with report.open("wt", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                (
                    "raw_input_genomes",
                    "raw_species",
                    "raw_genera",
                    "raw_families",
                    "excluded_by_length_qc",
                    "retained_genomes",
                    "retained_species",
                    "retained_genera",
                    "retained_families",
                )
            )
            writer.writerow(
                (
                    qc_total or genomes,
                    len(raw_species),
                    len(raw_genera),
                    len(raw_families),
                    qc_excluded,
                    genomes,
                    len(species),
                    len(genera),
                    len(families),
                )
            )

        print(
            "raw_input_genomes\traw_species\traw_genera\traw_families\t"
            "excluded_by_length_qc\tretained_genomes\tretained_species\t"
            "retained_genera\tretained_families"
        )
        print(
            f"{qc_total or genomes:,}\t{len(raw_species):,}\t{len(raw_genera):,}\t"
            f"{len(raw_families):,}\t{qc_excluded:,}\t{genomes:,}\t"
            f"{len(species):,}\t{len(genera):,}\t{len(families):,}"
        )
        print(f"\nReport: {report}")
        return 0
    except (
        CountError,
        pipeline.PipelineError,
        KeyError,
        OSError,
        UnicodeError,
        ValueError,
        json.JSONDecodeError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
