#!/usr/bin/env python3
"""Reproducible mitochondrial probe-generation pipeline.

Stages
------
1. Organize accession FASTA files as genus/species[/subgroup].
2. Before organizing or alignment, count each input genome's sequence length and
   exclude records shorter than the configurable minimum. Within each retained
   species (species-directory genomes) or subgroup, normalize strand and circular
   origin, run MAFFT, and generate probes/k-mers from shared alignment-coordinate
   starts. This prevents partial records, indels or different GenBank origins
   from emptying representative sets or creating permanent step-phase shifts.
   Windows containing non-ATCG bases are removed.
3. eMito-taxa-generate: representative species k-mers plus species- and
   genus-specific k-mers. For each species, its species-specific k-mers and its
   genus-specific k-mers are combined, intersected only with that species'
   probes, and merged into one Taxa-specific probe set. Only genomes placed
   directly under a species (empty subgroup label) participate; subgroup
   genomes are reserved for eMito-group-generate.
4. eMito-group-generate: representative subgroup k-mers, subgroup-specific
   k-mers, probe intersection, and within-subgroup sequence deduplication.
5. eMito-node-generate: choose one genome per requested taxonomy node (NC_
   preferred, otherwise seeded random selection) and tile probes.
6. Optionally run eMito-access and/or eMito-collapse independently for each of
   the three generation modes.
7. Merge the terminal output selected for each mode and, by default, deduplicate
   by sequence. Final deduplication can be disabled while retaining the merged
   record order and duplicate records.

Coordinates in generated headers are 1-based and inclusive. Reverse-complement
k-mer coordinates refer to the reverse-complement sequence itself, as requested.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
from array import array
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple


PIPELINE_VERSION = "0.1.1"
DEFAULT_METADATA = Path("metadata.tsv")
DEFAULT_FASTA_DIR = Path("fasta")
DEFAULT_OUTPUT_ROOT = Path("emito_output")
DEFAULT_NAMES_DMP = Path("taxdump/names.dmp")
DEFAULT_NODES_DMP = Path("taxdump/nodes.dmp")

ACCESSION_RE = re.compile(r"^[A-Za-z]{1,6}_?\d+(?:\.\d+)?$")
COORDINATE_RE = re.compile(
    r"^(?P<accession>[A-Za-z]{1,6}_?\d+(?:\.\d+)?)"
    r"\|start=(?P<start>\d+)\|end=(?P<end>\d+)"
    r"(?:\|strand=(?P<strand>[^\s]+))?"
)
DNA_RE = re.compile(r"^[ACGT]+$")
DNA_TO_BITS = {"A": 0, "C": 1, "G": 2, "T": 3}
RC_TRANSLATION = str.maketrans("ACGTacgt", "TGCAtgca")
NORMALIZED_COORDINATE_RE = re.compile(
    r"\|normalized_start=(?P<start>\d+)\|normalized_end=(?P<end>\d+)"
)


class PipelineError(RuntimeError):
    pass


class TaxonomyNodes:
    """Memory-efficient TaxID -> parent/rank storage for a full NCBI taxdump."""

    def __init__(self) -> None:
        self.parents = array("I", [0])
        self.rank_codes = bytearray(1)
        self.rank_to_code: Dict[str, int] = {"": 0}
        self.code_to_rank: List[str] = [""]
        self.node_count = 0

    def _ensure(self, taxid: int) -> None:
        if taxid < len(self.parents):
            return
        new_length = max(taxid + 1, len(self.parents) * 2)
        growth = new_length - len(self.parents)
        self.parents.extend(array("I", [0]) * growth)
        self.rank_codes.extend(b"\x00" * growth)

    def add(self, taxid: int, parent: int, rank: str) -> None:
        if taxid < 1 or parent < 1 or taxid > 100_000_000 or parent > 100_000_000:
            raise PipelineError(f"Unexpected TaxID in nodes.dmp: {taxid}, {parent}")
        self._ensure(taxid)
        if self.rank_codes[taxid] != 0:
            raise PipelineError(f"Duplicate TaxID {taxid} in nodes.dmp")
        code = self.rank_to_code.get(rank)
        if code is None:
            code = len(self.code_to_rank)
            if code > 255:
                raise PipelineError("More than 255 taxonomy ranks are not supported")
            self.rank_to_code[rank] = code
            self.code_to_rank.append(rank)
        self.parents[taxid] = parent
        self.rank_codes[taxid] = code
        self.node_count += 1

    def __contains__(self, taxid: object) -> bool:
        return (
            isinstance(taxid, int)
            and 0 < taxid < len(self.rank_codes)
            and self.rank_codes[taxid] != 0
        )

    def __getitem__(self, taxid: int) -> Tuple[int, str]:
        if taxid not in self:
            raise KeyError(taxid)
        return self.parents[taxid], self.code_to_rank[self.rank_codes[taxid]]

    def __len__(self) -> int:
        return self.node_count


@dataclass(frozen=True)
class MetadataRow:
    accession: str
    species_name: str
    taxid: int
    subgroup_label: str


@dataclass(frozen=True)
class Genome:
    accession: str
    species_name: str
    species_taxid: int
    taxid: int
    subgroup_label: str
    genus_name: str
    genus_taxid: int
    source_fasta: Path
    organized_fasta: Path
    sequence_length: int


@dataclass(frozen=True)
class WindowFiles:
    accession: str
    probe: Path
    kmer_forward: Path
    kmer_reverse: Path
    kmer_merged: Path


@dataclass(frozen=True)
class Target:
    target_id: str
    rank: str
    name: str
    taxid: Optional[int]
    genus_taxid: Optional[int] = None
    species_taxid: Optional[int] = None


@dataclass(frozen=True)
class WindowTask:
    """One alignment group's normalization, MAFFT and window-generation task."""

    group_id: str
    members: Tuple["WindowMember", ...]
    reference_accession: str
    normalized_input: str
    alignment_output: str
    normalization_manifest: str
    mafft_executable: str
    circular_anchor_length: int
    circular_min_anchor_length: int
    circular_min_anchor_hits: int
    window_length: int
    probe_step: int
    kmer_step: int


@dataclass(frozen=True)
class WindowMember:
    accession: str
    input_fasta: str
    probe_output: str
    kmer_forward_output: str
    kmer_reverse_output: str
    kmer_merged_output: str


@dataclass(frozen=True)
class GenomeLengthQC:
    accession: str
    species_name: str
    species_taxid: int
    taxid: int
    subgroup_label: str
    source_fasta: str
    sequence_length: int
    minimum_length: int
    retained: bool
    reason: str


@dataclass(frozen=True)
class AccessParameters:
    gc_min: float
    gc_max: float
    complexity_min: float
    complexity_max: float
    dimer_k: int
    dimer_threshold: float


@dataclass(frozen=True)
class AccessTask:
    input_fasta: str
    output_fasta: str
    output_tsv: str
    parameters: AccessParameters


@dataclass(frozen=True)
class CollapseTask:
    input_fasta: str
    output_fasta: str


@dataclass(frozen=True)
class NormalizationResult:
    sequence: str
    orientation: str
    rotation: int
    anchor_length: int
    anchor_hits: int
    modal_anchor_support: int
    attempted_hits: Tuple[Tuple[int, int], ...]


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise PipelineError(f"{label} is not a readable file: {path}")


def require_directory(path: Path, label: str) -> None:
    if not path.is_dir():
        raise PipelineError(f"{label} is not a readable directory: {path}")


def resolve_executable(command: str, label: str) -> str:
    candidate = Path(command)
    if candidate.parent != Path("."):
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate.resolve())
        raise PipelineError(f"{label} is not executable: {candidate}")
    resolved = shutil.which(command)
    if resolved is not None:
        return str(Path(resolved).resolve())
    environment_candidate = Path(sys.executable).resolve().parent / command
    if environment_candidate.is_file() and os.access(environment_candidate, os.X_OK):
        return str(environment_candidate)
    raise PipelineError(
        f"Cannot find executable {command!r} for {label}. Install MAFFT in the "
        "active environment or pass --mafft /absolute/path/to/mafft."
    )


def add_boolean_option(
    parser: argparse.ArgumentParser,
    name: str,
    default: bool,
    help_text: str,
) -> None:
    destination = name.replace("-", "_")
    group = parser.add_mutually_exclusive_group()
    group.add_argument(f"--{name}", dest=destination, action="store_true", help=help_text)
    group.add_argument(
        f"--no-{name}",
        dest=destination,
        action="store_false",
        help=f"Disable {help_text.lower()}",
    )
    parser.set_defaults(**{destination: default})


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Run all configurable eMito probe-generation modes.",
    )
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--fasta-dir", type=Path, default=DEFAULT_FASTA_DIR)
    parser.add_argument("--names-dmp", type=Path, default=DEFAULT_NAMES_DMP)
    parser.add_argument("--nodes-dmp", type=Path, default=DEFAULT_NODES_DMP)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--modes",
        default="taxa,group,node",
        help="Comma-separated generation modes chosen from taxa,group,node.",
    )
    parser.add_argument("--window-length", type=int, default=52)
    parser.add_argument("--probe-step", type=int, default=5)
    parser.add_argument("--kmer-step", type=int, default=1)
    parser.add_argument(
        "--min-genome-length",
        type=int,
        default=10_000,
        help=(
            "Exclude input mitochondrial FASTA records shorter than this many "
            "sequence characters before organizing, alignment and all modes."
        ),
    )
    parser.add_argument(
        "--mafft",
        default="mafft",
        help="MAFFT executable used for within-species/subgroup alignment.",
    )
    parser.add_argument(
        "--circular-anchor-length",
        type=int,
        default=31,
        help="Initial exact reference-anchor length used to normalize circular origins.",
    )
    parser.add_argument(
        "--circular-min-anchor-length",
        type=int,
        default=11,
        help="Shortest exact anchor allowed during automatic fallback.",
    )
    parser.add_argument(
        "--circular-min-anchor-hits",
        type=int,
        default=3,
        help="Minimum consistent reference anchors required for circular normalization.",
    )
    parser.add_argument(
        "--representative-fraction",
        type=float,
        default=0.75,
        help="Minimum individual-presence fraction when a group has >3 genomes.",
    )
    parser.add_argument(
        "--small-group-all-max",
        type=int,
        default=3,
        help="At or below this genome count, representative k-mers must occur in all.",
    )
    parser.add_argument(
        "--node-rank",
        default="family",
        help="NCBI taxonomy rank sampled by eMito-node-generate.",
    )
    parser.add_argument("--node-step", type=int, default=5)
    parser.add_argument("--random-seed", type=int, default=20250812)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--organize-method",
        choices=("copy", "symlink", "hardlink"),
        default="copy",
    )

    add_boolean_option(parser, "taxa-access", True, "Run eMito-access for taxa outputs.")
    add_boolean_option(parser, "taxa-collapse", False, "Run eMito-collapse for taxa outputs.")
    add_boolean_option(parser, "group-access", False, "Run eMito-access for group outputs.")
    add_boolean_option(parser, "group-collapse", False, "Run eMito-collapse for group outputs.")
    add_boolean_option(parser, "node-access", False, "Run eMito-access for node outputs.")
    add_boolean_option(parser, "node-collapse", False, "Run eMito-collapse for node outputs.")
    add_boolean_option(
        parser,
        "final-dedup",
        True,
        "Deduplicate the final merged probe FASTA by exact uppercase ATCG sequence.",
    )

    parser.add_argument(
        "--gc-min",
        type=float,
        default=35.0,
        help="Minimum GC percentage for eMito-access (eProbe default: 35).",
    )
    parser.add_argument(
        "--gc-max",
        type=float,
        default=65.0,
        help="Maximum GC percentage for eMito-access (eProbe default: 65).",
    )
    parser.add_argument("--complexity-min", type=float, default=0.0)
    parser.add_argument("--complexity-max", type=float, default=2.0)
    parser.add_argument(
        "--dimer-k",
        type=int,
        default=11,
        help="Reverse-complement k-mer size for the eProbe-compatible dimer score.",
    )
    parser.add_argument(
        "--dimer",
        "--dimer-threshold",
        dest="dimer_threshold",
        type=float,
        default=0.15,
        help=(
            "eProbe-compatible dimer threshold: 0<x<1 uses the x quantile "
            "(default 0.15), x>=1 is an absolute score, and x<=0 disables "
            "dimer filtering."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Remove and recreate the exact --output-root if it already exists.",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate metadata, FASTA and taxonomy inputs without generating output.",
    )
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> Set[str]:
    modes = {item.strip().lower() for item in args.modes.split(",") if item.strip()}
    unknown = modes - {"taxa", "group", "node"}
    if unknown or not modes:
        raise PipelineError(f"Invalid --modes value; unknown={sorted(unknown)}")
    for label in (
        "window_length",
        "probe_step",
        "kmer_step",
        "node_step",
        "threads",
        "circular_anchor_length",
        "circular_min_anchor_length",
        "circular_min_anchor_hits",
        "min_genome_length",
    ):
        if getattr(args, label) < 1:
            raise PipelineError(f"--{label.replace('_', '-')} must be >= 1")
    if args.circular_min_anchor_length > args.circular_anchor_length:
        raise PipelineError(
            "--circular-min-anchor-length cannot exceed --circular-anchor-length"
        )
    if not 0 < args.representative_fraction <= 1:
        raise PipelineError("--representative-fraction must be in (0,1]")
    if args.small_group_all_max < 0:
        raise PipelineError("--small-group-all-max must be >= 0")
    if not 0 <= args.gc_min <= args.gc_max <= 100:
        raise PipelineError("GC bounds must satisfy 0 <= min <= max <= 100")
    if args.complexity_min > args.complexity_max:
        raise PipelineError("Complexity minimum exceeds maximum")
    if args.dimer_k < 1:
        raise PipelineError("--dimer-k must be >= 1")
    return modes


def load_nodes(path: Path) -> TaxonomyNodes:
    require_file(path, "nodes.dmp")
    nodes = TaxonomyNodes()
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split("|")
            if len(fields) < 3:
                raise PipelineError(f"Malformed nodes.dmp line {line_number}")
            try:
                nodes.add(int(fields[0].strip()), int(fields[1].strip()), fields[2].strip())
            except ValueError as exc:
                raise PipelineError(f"Invalid TaxID at nodes.dmp line {line_number}") from exc
    log(f"Loaded {len(nodes):,} taxonomy nodes")
    return nodes


def ancestor_at_rank(taxid: int, rank: str, nodes: TaxonomyNodes) -> int:
    current = taxid
    visited: Set[int] = set()
    while current not in visited:
        visited.add(current)
        try:
            parent, current_rank = nodes[current]
        except KeyError as exc:
            raise PipelineError(f"TaxID {current} is absent from nodes.dmp") from exc
        if current_rank == rank:
            return current
        if parent == current:
            break
        current = parent
    raise PipelineError(f"TaxID {taxid} has no ancestor at rank {rank!r}")


def load_names(path: Path, wanted: Set[int]) -> Dict[int, str]:
    require_file(path, "names.dmp")
    names: Dict[int, str] = {}
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split("|")
            if len(fields) < 4:
                raise PipelineError(f"Malformed names.dmp line {line_number}")
            try:
                taxid = int(fields[0].strip())
            except ValueError as exc:
                raise PipelineError(f"Invalid TaxID at names.dmp line {line_number}") from exc
            if taxid in wanted and fields[3].strip() == "scientific name":
                names[taxid] = fields[1].strip()
    missing = sorted(wanted - set(names))
    if missing:
        raise PipelineError(f"Scientific names missing for TaxID(s): {missing[:20]}")
    return names


def load_metadata(path: Path) -> List[MetadataRow]:
    require_file(path, "metadata")
    rows: List[MetadataRow] = []
    seen: Set[str] = set()
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"accession_id", "species_name", "taxid", "subgroup_label"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise PipelineError(
                f"Metadata must contain columns {sorted(required)}; got {reader.fieldnames}"
            )
        for line_number, raw in enumerate(reader, start=2):
            accession = raw["accession_id"].strip().upper()
            if not ACCESSION_RE.fullmatch(accession):
                raise PipelineError(f"Invalid accession at metadata line {line_number}: {accession}")
            if accession in seen:
                raise PipelineError(f"Duplicate accession in metadata: {accession}")
            seen.add(accession)
            try:
                taxid = int(raw["taxid"].strip())
            except ValueError as exc:
                raise PipelineError(f"Invalid TaxID at metadata line {line_number}") from exc
            species_name = raw["species_name"].strip()
            if not species_name:
                raise PipelineError(f"Empty species name at metadata line {line_number}")
            rows.append(
                MetadataRow(
                    accession=accession,
                    species_name=species_name,
                    taxid=taxid,
                    subgroup_label=raw["subgroup_label"].strip(),
                )
            )
    if not rows:
        raise PipelineError("Metadata contains no genomes")
    log(f"Loaded {len(rows):,} metadata rows")
    return rows


def safe_component(name: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    value = re.sub(r"_+", "_", value).strip("._-")
    if not value:
        raise PipelineError(f"Taxonomy name cannot form a safe directory name: {name!r}")
    return value


def read_single_fasta(path: Path) -> Tuple[str, str]:
    require_file(path, "genome FASTA")
    header: Optional[str] = None
    parts: List[str] = []
    with path.open("rt", encoding="ascii", errors="strict") as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.strip()
            if not text:
                continue
            if text.startswith(">"):
                if header is not None:
                    raise PipelineError(f"Multiple records found in genome FASTA: {path}")
                header = text[1:].strip()
            else:
                if header is None:
                    raise PipelineError(f"Sequence before header at {path}:{line_number}")
                parts.append(re.sub(r"\s+", "", text).upper())
    sequence = "".join(parts)
    if header is None or not sequence:
        raise PipelineError(f"Empty FASTA: {path}")
    if not re.fullmatch(r"[A-Z*.-]+", sequence):
        raise PipelineError(f"Invalid sequence characters in genome FASTA: {path}")
    return header, sequence


def biological_sequence_length(sequence: str) -> int:
    """Count sequence letters while ignoring any pre-existing alignment gaps."""
    return sum(1 for base in sequence if base.isalpha())


def build_genomes(
    rows: Sequence[MetadataRow],
    fasta_dir: Path,
    nodes: TaxonomyNodes,
    names_dmp: Path,
    organized_root: Path,
) -> List[Genome]:
    require_directory(fasta_dir, "input FASTA directory")
    lineage: Dict[int, Tuple[int, int]] = {}
    wanted: Set[int] = set()
    for row in rows:
        species_taxid = ancestor_at_rank(row.taxid, "species", nodes)
        genus_taxid = ancestor_at_rank(row.taxid, "genus", nodes)
        lineage[row.taxid] = (species_taxid, genus_taxid)
        wanted.update((species_taxid, genus_taxid))
    names = load_names(names_dmp, wanted)

    genomes: List[Genome] = []
    organized_seen: Dict[Path, str] = {}
    for row in rows:
        species_taxid, genus_taxid = lineage[row.taxid]
        taxonomy_species_name = names[species_taxid]
        if taxonomy_species_name != row.species_name:
            raise PipelineError(
                f"Metadata species mismatch for {row.accession}: metadata={row.species_name!r}, "
                f"taxonomy={taxonomy_species_name!r} (species TaxID {species_taxid})"
            )
        source = fasta_dir / f"{row.accession}.fasta"
        _, source_sequence = read_single_fasta(source)
        sequence_length = biological_sequence_length(source_sequence)
        parts = [safe_component(names[genus_taxid]), safe_component(row.species_name)]
        if row.subgroup_label:
            parts.append(safe_component(row.subgroup_label))
        destination = organized_root.joinpath(*parts, f"{row.accession}.fasta")
        previous = organized_seen.get(destination)
        if previous is not None and previous != row.accession:
            raise PipelineError(f"Organized path collision: {destination}")
        organized_seen[destination] = row.accession
        genomes.append(
            Genome(
                accession=row.accession,
                species_name=row.species_name,
                species_taxid=species_taxid,
                taxid=row.taxid,
                subgroup_label=row.subgroup_label,
                genus_name=names[genus_taxid],
                genus_taxid=genus_taxid,
                source_fasta=source,
                organized_fasta=destination,
                sequence_length=sequence_length,
            )
        )
    return sorted(genomes, key=lambda g: (g.genus_name, g.species_name, g.subgroup_label, g.accession))


def filter_genomes_by_length(
    genomes: Sequence[Genome], minimum_length: int
) -> Tuple[List[Genome], List[GenomeLengthQC]]:
    retained: List[Genome] = []
    rows: List[GenomeLengthQC] = []
    for genome in genomes:
        keep = genome.sequence_length >= minimum_length
        reason = "retained" if keep else f"sequence_length_below_{minimum_length}"
        rows.append(
            GenomeLengthQC(
                accession=genome.accession,
                species_name=genome.species_name,
                species_taxid=genome.species_taxid,
                taxid=genome.taxid,
                subgroup_label=genome.subgroup_label,
                source_fasta=str(genome.source_fasta),
                sequence_length=genome.sequence_length,
                minimum_length=minimum_length,
                retained=keep,
                reason=reason,
            )
        )
        if keep:
            retained.append(genome)
    return retained, rows


def write_genome_length_qc(
    output_root: Path, rows: Sequence[GenomeLengthQC]
) -> None:
    qc_root = output_root / "00_input_qc"
    qc_root.mkdir(parents=True, exist_ok=True)
    header = (
        "accession_id",
        "species_name",
        "species_taxid",
        "actual_taxid",
        "subgroup_label",
        "source_fasta",
        "sequence_length",
        "minimum_length",
        "status",
        "reason",
    )

    def write_rows(path: Path, selected: Iterable[GenomeLengthQC]) -> None:
        with path.open("wt", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            for row in selected:
                writer.writerow(
                    (
                        row.accession,
                        row.species_name,
                        row.species_taxid,
                        row.taxid,
                        row.subgroup_label,
                        row.source_fasta,
                        row.sequence_length,
                        row.minimum_length,
                        "retained" if row.retained else "excluded",
                        row.reason,
                    )
                )

    ordered = sorted(rows, key=lambda row: (row.species_name, row.subgroup_label, row.accession))
    write_rows(qc_root / "genome_length_qc.tsv", ordered)
    write_rows(qc_root / "excluded_short_genomes.tsv", (row for row in ordered if not row.retained))


def organize_genomes(genomes: Sequence[Genome], method: str) -> None:
    for genome in genomes:
        destination = genome.organized_fasta
        destination.parent.mkdir(parents=True, exist_ok=True)
        if method == "copy":
            shutil.copy2(genome.source_fasta, destination)
        elif method == "symlink":
            destination.symlink_to(genome.source_fasta.resolve())
        elif method == "hardlink":
            os.link(genome.source_fasta, destination)
        else:
            raise AssertionError(method)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(RC_TRANSLATION)[::-1]


def iter_windows(sequence: str, length: int, step: int) -> Iterator[Tuple[int, int, str]]:
    if len(sequence) < length:
        return
    for offset in range(0, len(sequence) - length + 1, step):
        window = sequence[offset : offset + length]
        if DNA_RE.fullmatch(window):
            yield offset + 1, offset + length, window


def forward_header(accession: str, start: int, end: int, kind: str) -> str:
    return f"{accession}|start={start}|end={end}|strand=Forward|type={kind}"


def reverse_header(accession: str, start: int, end: int) -> str:
    return (
        f"{accession}|start={start}|end={end}|strand=ReverseComplement|type=kmer "
        "[Reverse Complement]"
    )


def circular_unique_kmer_positions(sequence: str, k: int) -> Dict[str, Optional[int]]:
    if len(sequence) < k:
        return {}
    doubled = sequence + sequence[: k - 1]
    result: Dict[str, Optional[int]] = {}
    for offset in range(len(sequence)):
        kmer = doubled[offset : offset + k]
        if not DNA_RE.fullmatch(kmer):
            continue
        if kmer in result:
            result[kmer] = None
        else:
            result[kmer] = offset
    return result


def normalization_at_anchor_length(
    reference: str,
    sequence: str,
    anchor_length: int,
) -> NormalizationResult:
    if len(reference) < anchor_length or len(sequence) < anchor_length:
        raise PipelineError(
            f"Genome is shorter than circular anchor length {anchor_length}"
        )
    stride = max(1, anchor_length // 2)
    reference_doubled = reference + reference[: anchor_length - 1]
    anchors = [
        (offset, reference_doubled[offset : offset + anchor_length])
        for offset in range(0, len(reference), stride)
        if DNA_RE.fullmatch(reference_doubled[offset : offset + anchor_length])
    ]
    if not anchors:
        raise PipelineError("Reference has no ATCG-only circular-normalization anchors")

    candidates = (("Forward", sequence), ("ReverseComplement", reverse_complement(sequence)))
    scored: List[Tuple[int, int, int, int, int, str, str]] = []
    for orientation, oriented in candidates:
        positions = circular_unique_kmer_positions(oriented, anchor_length)
        origin_estimates: Counter[int] = Counter()
        matched: List[Tuple[int, int, int]] = []
        hits = 0
        for reference_offset, anchor in anchors:
            candidate_offset = positions.get(anchor)
            if candidate_offset is None:
                continue
            hits += 1
            signed_reference_offset = (
                reference_offset
                if reference_offset <= len(reference) // 2
                else reference_offset - len(reference)
            )
            estimate = (candidate_offset - signed_reference_offset) % len(oriented)
            distance_to_origin = abs(signed_reference_offset)
            origin_estimates[estimate] += 1
            matched.append((distance_to_origin, reference_offset, estimate))
        if matched:
            _, _, rotation = min(matched)
            modal_support = max(origin_estimates.values())
            nearest_distance = min(item[0] for item in matched)
        else:
            rotation, modal_support, nearest_distance = 0, 0, len(reference)
        scored.append(
            (
                hits,
                modal_support,
                int(orientation == "Forward"),
                -nearest_distance,
                -rotation,
                orientation,
                oriented,
            )
        )

    hits, support, _, _, negative_rotation, orientation, oriented = max(scored)
    rotation = -negative_rotation
    normalized = oriented[rotation:] + oriented[:rotation]
    return NormalizationResult(
        sequence=normalized,
        orientation=orientation,
        rotation=rotation,
        anchor_length=anchor_length,
        anchor_hits=hits,
        modal_anchor_support=support,
        attempted_hits=((anchor_length, hits),),
    )


def circular_anchor_lengths(initial: int, minimum: int) -> List[int]:
    lengths: List[int] = []
    current = initial
    while current > minimum:
        lengths.append(current)
        current -= 4
    lengths.append(minimum)
    return list(dict.fromkeys(lengths))


def circular_normalize_to_reference(
    reference: str,
    sequence: str,
    anchor_length: int,
    minimum_anchor_length: int,
    minimum_anchor_hits: int,
) -> NormalizationResult:
    """Normalize with exact anchors, automatically falling back in length."""
    attempts: List[Tuple[int, int]] = []
    for candidate_length in circular_anchor_lengths(anchor_length, minimum_anchor_length):
        if len(reference) < candidate_length or len(sequence) < candidate_length:
            attempts.append((candidate_length, 0))
            continue
        result = normalization_at_anchor_length(reference, sequence, candidate_length)
        attempts.append((candidate_length, result.anchor_hits))
        if result.anchor_hits >= minimum_anchor_hits:
            return NormalizationResult(
                sequence=result.sequence,
                orientation=result.orientation,
                rotation=result.rotation,
                anchor_length=result.anchor_length,
                anchor_hits=result.anchor_hits,
                modal_anchor_support=result.modal_anchor_support,
                attempted_hits=tuple(attempts),
            )
    attempt_text = ", ".join(f"{length}mer={hits}" for length, hits in attempts)
    raise PipelineError(
        "Cannot normalize circular origin/orientation after adaptive exact-anchor "
        f"fallback ({attempt_text}); minimum required hits={minimum_anchor_hits}"
    )


def original_coordinate(
    normalized_offset: int,
    sequence_length: int,
    orientation: str,
    rotation: int,
) -> int:
    oriented_offset = (rotation + normalized_offset) % sequence_length
    if orientation == "Forward":
        return oriented_offset + 1
    if orientation == "ReverseComplement":
        return sequence_length - oriented_offset
    raise AssertionError(orientation)


def read_alignment(path: Path) -> Dict[str, str]:
    result: Dict[str, str] = {}
    length: Optional[int] = None
    for header, sequence in iter_fasta(path):
        accession = header.split()[0]
        if accession in result:
            raise PipelineError(f"Duplicate accession in alignment: {accession}")
        if not re.fullmatch(r"[A-Z*?.-]+", sequence):
            raise PipelineError(f"Unexpected alignment characters for {accession}: {path}")
        if length is None:
            length = len(sequence)
        elif len(sequence) != length:
            raise PipelineError(f"Unequal aligned sequence lengths in {path}")
        result[accession] = sequence
    if not result:
        raise PipelineError(f"MAFFT produced an empty alignment: {path}")
    return result


def aligned_reference_columns(reference_alignment: str) -> List[int]:
    return [index for index, base in enumerate(reference_alignment) if base != "-"]


def alignment_coordinate_maps(
    aligned_sequence: str,
) -> Tuple[str, List[int], List[int]]:
    ungapped: List[str] = []
    column_to_offset: List[int] = []
    base_columns: List[int] = []
    offset = 0
    for column, base in enumerate(aligned_sequence):
        column_to_offset.append(offset)
        if base != "-":
            ungapped.append(base)
            base_columns.append(column)
            offset += 1
    return "".join(ungapped), column_to_offset, base_columns


def extract_aligned_window(
    aligned_sequence: str,
    ungapped_sequence: str,
    column_to_offset: Sequence[int],
    base_columns: Sequence[int],
    alignment_start: int,
    length: int,
) -> Optional[Tuple[str, int, int, int, int]]:
    if aligned_sequence[alignment_start] == "-":
        return None
    normalized_start = column_to_offset[alignment_start]
    normalized_end = normalized_start + length - 1
    if normalized_end >= len(ungapped_sequence):
        return None
    sequence = ungapped_sequence[normalized_start : normalized_end + 1]
    if not DNA_RE.fullmatch(sequence):
        return None
    return (
        sequence,
        normalized_start,
        normalized_end,
        alignment_start,
        base_columns[normalized_end],
    )


def aligned_forward_header(
    accession: str,
    original_start: int,
    original_end: int,
    normalized_start: int,
    normalized_end: int,
    alignment_start: int,
    alignment_end: int,
    kind: str,
    source_orientation: str,
) -> str:
    return (
        f"{accession}|start={original_start}|end={original_end}|strand=Forward"
        f"|normalized_start={normalized_start}|normalized_end={normalized_end}"
        f"|alignment_start={alignment_start}|alignment_end={alignment_end}"
        f"|source_orientation={source_orientation}|type={kind}"
    )


def aligned_reverse_header(
    accession: str,
    reverse_start: int,
    reverse_end: int,
    normalized_start: int,
    normalized_end: int,
    alignment_start: int,
    alignment_end: int,
    source_orientation: str,
) -> str:
    return (
        f"{accession}|start={reverse_start}|end={reverse_end}"
        "|strand=ReverseComplement"
        f"|normalized_start={normalized_start}|normalized_end={normalized_end}"
        f"|alignment_start={alignment_start}|alignment_end={alignment_end}"
        f"|source_orientation={source_orientation}|type=kmer [Reverse Complement]"
    )


def run_window_task(task: WindowTask) -> Dict[str, int]:
    members = {member.accession: member for member in task.members}
    raw_sequences = {
        accession: read_single_fasta(Path(member.input_fasta))[1]
        for accession, member in members.items()
    }
    reference = raw_sequences[task.reference_accession]
    normalized_input = Path(task.normalized_input)
    alignment_output = Path(task.alignment_output)
    normalization_manifest = Path(task.normalization_manifest)
    normalization_failure = normalization_manifest.with_name("normalization.failed.tsv")
    mafft_log = alignment_output.with_suffix(".mafft.log")
    for path in (
        normalized_input,
        alignment_output,
        normalization_manifest,
        normalization_failure,
        mafft_log,
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
    normalized: Dict[str, str] = {}
    normalization: Dict[str, NormalizationResult] = {}
    for accession in sorted(members):
        if accession == task.reference_accession:
            normalized[accession] = reference
            normalization[accession] = NormalizationResult(
                sequence=reference,
                orientation="Forward",
                rotation=0,
                anchor_length=task.circular_anchor_length,
                anchor_hits=0,
                modal_anchor_support=0,
                attempted_hits=(),
            )
        else:
            try:
                result = circular_normalize_to_reference(
                    reference,
                    raw_sequences[accession],
                    task.circular_anchor_length,
                    task.circular_min_anchor_length,
                    task.circular_min_anchor_hits,
                )
            except PipelineError as exc:
                with normalization_failure.open(
                    "wt", encoding="utf-8", newline=""
                ) as handle:
                    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(
                        (
                            "group_id",
                            "reference_accession",
                            "failed_accession",
                            "reference_length",
                            "accession_length",
                            "error",
                        )
                    )
                    writer.writerow(
                        (
                            task.group_id,
                            task.reference_accession,
                            accession,
                            len(reference),
                            len(raw_sequences[accession]),
                            str(exc),
                        )
                    )
                raise PipelineError(
                    f"Circular normalization failed: group={task.group_id}, "
                    f"reference={task.reference_accession}, accession={accession}, "
                    f"reference_length={len(reference)}, "
                    f"accession_length={len(raw_sequences[accession])}: {exc}"
                ) from exc
            normalized[accession] = result.sequence
            normalization[accession] = result

    with normalized_input.open("wt", encoding="ascii") as handle:
        for accession in sorted(normalized):
            handle.write(f">{accession}\n{normalized[accession]}\n")
    with normalization_manifest.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "group_id",
                "accession_id",
                "reference_accession",
                "source_length",
                "source_orientation",
                "circular_rotation_0_based",
                "anchor_length_used",
                "anchor_hits",
                "modal_anchor_support",
                "adaptive_anchor_attempts",
            )
        )
        for accession in sorted(normalized):
            result = normalization[accession]
            writer.writerow(
                (
                    task.group_id,
                    accession,
                    task.reference_accession,
                    len(raw_sequences[accession]),
                    result.orientation,
                    result.rotation,
                    result.anchor_length,
                    result.anchor_hits,
                    result.modal_anchor_support,
                    ",".join(
                        f"{length}:{hits}" for length, hits in result.attempted_hits
                    ),
                )
            )

    if len(members) == 1:
        shutil.copyfile(normalized_input, alignment_output)
        mafft_log.write_text("Singleton group: MAFFT was not required.\n", encoding="utf-8")
    else:
        with alignment_output.open("wt", encoding="ascii") as stdout, mafft_log.open(
            "wt", encoding="utf-8"
        ) as stderr:
            completed = subprocess.run(
                [task.mafft_executable, "--auto", "--thread", "1", str(normalized_input)],
                stdout=stdout,
                stderr=stderr,
                text=True,
                check=False,
            )
        if completed.returncode != 0:
            raise PipelineError(
                f"MAFFT failed for {task.group_id} with exit code "
                f"{completed.returncode}; see {mafft_log}"
            )

    aligned = read_alignment(alignment_output)
    if set(aligned) != set(members):
        raise PipelineError(
            f"MAFFT accession mismatch for {task.group_id}: "
            f"expected={sorted(members)}, observed={sorted(aligned)}"
        )
    for accession, aligned_sequence in aligned.items():
        if aligned_sequence.replace("-", "") != normalized[accession]:
            raise PipelineError(
                f"MAFFT changed or mismatched the ungapped sequence for {accession} "
                f"in {task.group_id}"
            )
    reference_columns = aligned_reference_columns(aligned[task.reference_accession])
    if len(reference_columns) != len(reference):
        raise PipelineError(f"Reference length changed in MAFFT output: {task.group_id}")

    coordinate_maps = {
        accession: alignment_coordinate_maps(aligned_sequence)
        for accession, aligned_sequence in aligned.items()
    }
    total_probe = 0
    total_forward = 0
    total_reverse = 0
    for accession in sorted(members):
        member = members[accession]
        outputs = [
            Path(member.probe_output),
            Path(member.kmer_forward_output),
            Path(member.kmer_reverse_output),
            Path(member.kmer_merged_output),
        ]
        for output in outputs:
            output.parent.mkdir(parents=True, exist_ok=True)
        normalization_result = normalization[accession]
        orientation = normalization_result.orientation
        rotation = normalization_result.rotation
        source_length = len(raw_sequences[accession])
        ungapped, column_to_offset, base_columns = coordinate_maps[accession]

        with outputs[0].open("wt", encoding="ascii") as probe_handle:
            for reference_offset in range(
                0, len(reference_columns) - task.window_length + 1, task.probe_step
            ):
                window = extract_aligned_window(
                    aligned[accession],
                    ungapped,
                    column_to_offset,
                    base_columns,
                    reference_columns[reference_offset],
                    task.window_length,
                )
                if window is None:
                    continue
                sequence, normalized_start, normalized_end, aligned_start, aligned_end = window
                original_start = original_coordinate(
                    normalized_start, source_length, orientation, rotation
                )
                original_end = original_coordinate(
                    normalized_end, source_length, orientation, rotation
                )
                header = aligned_forward_header(
                    accession,
                    original_start,
                    original_end,
                    normalized_start + 1,
                    normalized_end + 1,
                    aligned_start + 1,
                    aligned_end + 1,
                    "probe",
                    orientation,
                )
                probe_handle.write(f">{header}\n{sequence}\n")
                total_probe += 1

        with outputs[1].open("wt", encoding="ascii") as forward_handle, outputs[2].open(
            "wt", encoding="ascii"
        ) as reverse_handle, outputs[3].open("wt", encoding="ascii") as merged_handle:
            for reference_offset in range(
                0, len(reference_columns) - task.window_length + 1, task.kmer_step
            ):
                window = extract_aligned_window(
                    aligned[accession],
                    ungapped,
                    column_to_offset,
                    base_columns,
                    reference_columns[reference_offset],
                    task.window_length,
                )
                if window is None:
                    continue
                sequence, normalized_start, normalized_end, aligned_start, aligned_end = window
                original_start = original_coordinate(
                    normalized_start, source_length, orientation, rotation
                )
                original_end = original_coordinate(
                    normalized_end, source_length, orientation, rotation
                )
                forward = aligned_forward_header(
                    accession,
                    original_start,
                    original_end,
                    normalized_start + 1,
                    normalized_end + 1,
                    aligned_start + 1,
                    aligned_end + 1,
                    "kmer",
                    orientation,
                )
                forward_text = f">{forward}\n{sequence}\n"
                forward_handle.write(forward_text)
                merged_handle.write(forward_text)
                total_forward += 1

                reverse_start = source_length - normalized_end
                reverse_end = source_length - normalized_start
                reverse = aligned_reverse_header(
                    accession,
                    reverse_start,
                    reverse_end,
                    normalized_start + 1,
                    normalized_end + 1,
                    aligned_start + 1,
                    aligned_end + 1,
                    orientation,
                )
                reverse_sequence = reverse_complement(sequence)
                reverse_text = f">{reverse}\n{reverse_sequence}\n"
                reverse_handle.write(reverse_text)
                merged_handle.write(reverse_text)
                total_reverse += 1
    return {
        "groups": 1,
        "genomes": len(members),
        "probe": total_probe,
        "kmer_forward": total_forward,
        "kmer_reverse": total_reverse,
    }


def build_window_tasks(
    genomes: Sequence[Genome],
    output_root: Path,
    length: int,
    probe_step: int,
    kmer_step: int,
    mafft_executable: str,
    circular_anchor_length: int,
    circular_min_anchor_length: int,
    circular_min_anchor_hits: int,
) -> Tuple[List[WindowTask], Dict[str, WindowFiles]]:
    grouped: Dict[Tuple[str, int, int, str], List[Genome]] = defaultdict(list)
    for genome in genomes:
        if genome.subgroup_label:
            key = ("subgroup", genome.species_taxid, genome.taxid, genome.subgroup_label)
        else:
            key = ("species", genome.species_taxid, genome.species_taxid, "")
        grouped[key].append(genome)

    tasks: List[WindowTask] = []
    files: Dict[str, WindowFiles] = {}
    for key in sorted(grouped):
        group = sorted(grouped[key], key=lambda item: item.accession)
        preferred = [genome for genome in group if genome.accession.startswith("NC_")]
        reference = (preferred or group)[0]
        label = key[3] if key[0] == "subgroup" else reference.species_name
        digest = hashlib.sha1(f"{key}".encode("utf-8")).hexdigest()[:10]
        group_id = f"{key[0]}_{safe_component(label)}_{key[2]}_{digest}"
        alignment_root = output_root / "02_alignments" / group_id
        members: List[WindowMember] = []
        for genome in group:
            relative_parent = genome.organized_fasta.parent.relative_to(
                output_root / "01_organized_genomes"
            )
            base_dir = output_root / "02_windows" / relative_parent
            window_files = WindowFiles(
                accession=genome.accession,
                probe=base_dir / f"{genome.accession}.probe.fasta",
                kmer_forward=base_dir / f"{genome.accession}.kmer.forward.fasta",
                kmer_reverse=base_dir / f"{genome.accession}.kmer.reverse_complement.fasta",
                kmer_merged=base_dir / f"{genome.accession}.kmer.merged.fasta",
            )
            files[genome.accession] = window_files
            members.append(
                WindowMember(
                    accession=genome.accession,
                    input_fasta=str(genome.organized_fasta),
                    probe_output=str(window_files.probe),
                    kmer_forward_output=str(window_files.kmer_forward),
                    kmer_reverse_output=str(window_files.kmer_reverse),
                    kmer_merged_output=str(window_files.kmer_merged),
                )
            )
        tasks.append(
            WindowTask(
                group_id=group_id,
                members=tuple(members),
                reference_accession=reference.accession,
                normalized_input=str(alignment_root / "normalized_input.fasta"),
                alignment_output=str(alignment_root / "alignment.fasta"),
                normalization_manifest=str(alignment_root / "normalization.tsv"),
                mafft_executable=mafft_executable,
                circular_anchor_length=circular_anchor_length,
                circular_min_anchor_length=circular_min_anchor_length,
                circular_min_anchor_hits=circular_min_anchor_hits,
                window_length=length,
                probe_step=probe_step,
                kmer_step=kmer_step,
            )
        )
    return tasks, files


def write_alignment_group_manifest(path: Path, tasks: Sequence[WindowTask]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "group_id",
                "reference_accession",
                "genome_count",
                "member_accessions",
                "normalized_input_fasta",
                "alignment_fasta",
                "normalization_manifest",
            )
        )
        for task in sorted(tasks, key=lambda item: item.group_id):
            writer.writerow(
                (
                    task.group_id,
                    task.reference_accession,
                    len(task.members),
                    ",".join(member.accession for member in task.members),
                    task.normalized_input,
                    task.alignment_output,
                    task.normalization_manifest,
                )
            )


def execute_process_tasks(function, tasks: Sequence[object], threads: int, label: str) -> List[object]:
    if not tasks:
        return []
    results: List[object] = []
    if threads == 1:
        for index, task in enumerate(tasks, start=1):
            results.append(function(task))
            if index % 100 == 0 or index == len(tasks):
                log(f"{label}: {index:,}/{len(tasks):,}")
        return results
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(function, task): task for task in tasks}
        completed = 0
        for future in as_completed(futures):
            results.append(future.result())
            completed += 1
            if completed % 100 == 0 or completed == len(tasks):
                log(f"{label}: {completed:,}/{len(tasks):,}")
    return results


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
                    sequence = "".join(parts).upper()
                    if not sequence:
                        raise PipelineError(f"Empty record in {path}")
                    yield header, sequence
                header = text[1:].strip()
                parts = []
            else:
                if header is None:
                    raise PipelineError(f"Sequence before header at {path}:{line_number}")
                parts.append(text)
    if header is not None:
        sequence = "".join(parts).upper()
        if not sequence:
            raise PipelineError(f"Empty record in {path}")
        yield header, sequence


def write_fasta(path: Path, records: Iterable[Tuple[str, str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("wt", encoding="ascii") as handle:
        for header, sequence in records:
            if not DNA_RE.fullmatch(sequence):
                continue
            handle.write(f">{header}\n{sequence}\n")
            count += 1
    return count


def sequence_dictionary(path: Path) -> Dict[str, str]:
    result: Dict[str, str] = {}
    for header, sequence in iter_fasta(path):
        if DNA_RE.fullmatch(sequence) and sequence not in result:
            result[sequence] = header
    return result


def representative_threshold(
    individual_count: int, fraction: float, small_group_all_max: int
) -> int:
    if individual_count < 1:
        raise PipelineError("Representative k-mer group has zero individuals")
    if individual_count <= small_group_all_max:
        return individual_count
    return max(1, math.ceil(fraction * individual_count))


def representative_kmers(
    kmer_files: Sequence[Path], fraction: float, small_group_all_max: int
) -> Tuple[Dict[str, str], int]:
    threshold = representative_threshold(len(kmer_files), fraction, small_group_all_max)
    counts: Dict[str, int] = {}
    first_header: Dict[str, str] = {}
    for path in sorted(kmer_files, key=lambda p: str(p)):
        individual = sequence_dictionary(path)
        for sequence, header in individual.items():
            counts[sequence] = counts.get(sequence, 0) + 1
            first_header.setdefault(sequence, header)
    representative = {
        sequence: first_header[sequence]
        for sequence, count in counts.items()
        if count >= threshold
    }
    return representative, threshold


def intersect_probe_file(input_path: Path, kmers: Set[str], output: Path) -> int:
    return write_fasta(
        output,
        (
            (header, sequence)
            for header, sequence in iter_fasta(input_path)
            if sequence in kmers and DNA_RE.fullmatch(sequence)
        ),
    )


def merge_probe_files_deduplicated(probe_files: Sequence[Path], output: Path) -> int:
    seen: Set[str] = set()

    def records() -> Iterator[Tuple[str, str]]:
        for path in sorted(probe_files, key=lambda p: str(p)):
            for header, sequence in iter_fasta(path):
                if sequence not in seen and DNA_RE.fullmatch(sequence):
                    seen.add(sequence)
                    yield header, sequence

    return write_fasta(output, records())


def intersect_individuals_then_merge(
    genomes: Sequence[Genome],
    windows: Mapping[str, WindowFiles],
    kmers: Set[str],
    intersection_root: Path,
    merged_output: Path,
) -> int:
    individual_outputs: List[Path] = []
    for genome in sorted(genomes, key=lambda item: item.accession):
        output = intersection_root / f"{genome.accession}.intersection.probe.fasta"
        intersect_probe_file(windows[genome.accession].probe, kmers, output)
        individual_outputs.append(output)
    return merge_probe_files_deduplicated(individual_outputs, merged_output)


def make_species_target(genome: Genome) -> Target:
    return Target(
        target_id=f"species_taxid_{genome.species_taxid}",
        rank="species",
        name=genome.species_name,
        taxid=genome.species_taxid,
        genus_taxid=genome.genus_taxid,
        species_taxid=genome.species_taxid,
    )


def make_genus_target(genome: Genome) -> Target:
    return Target(
        target_id=f"genus_taxid_{genome.genus_taxid}",
        rank="genus",
        name=genome.genus_name,
        taxid=genome.genus_taxid,
        genus_taxid=genome.genus_taxid,
    )


def make_group_target(genome: Genome) -> Target:
    digest = hashlib.sha1(genome.subgroup_label.encode("utf-8")).hexdigest()[:10]
    return Target(
        target_id=(
            f"group_species_{genome.species_taxid}_taxid_{genome.taxid}_{digest}"
        ),
        rank="subgroup",
        name=genome.subgroup_label,
        taxid=genome.taxid,
        genus_taxid=genome.genus_taxid,
        species_taxid=genome.species_taxid,
    )


def target_file_component(target: Target) -> str:
    return f"{safe_component(target.name)}__{target.target_id}"


def write_target_manifest(path: Path, rows: Sequence[Tuple[Target, Path, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("target_id", "rank", "name", "taxid", "output_fasta", "sequence_count"))
        for target, output, count in rows:
            writer.writerow(
                (target.target_id, target.rank, target.name, target.taxid or "", output, count)
            )


def run_taxa_mode(
    genomes: Sequence[Genome],
    windows: Mapping[str, WindowFiles],
    mode_root: Path,
    fraction: float,
    small_group_all_max: int,
) -> List[Tuple[Target, Path]]:
    # Species-level taxa mode deliberately uses only FASTAs located directly in
    # the species directory. Subgroup-labelled FASTAs belong exclusively to the
    # group mode even though metadata maps them to the same species.
    selected = [g for g in genomes if not g.subgroup_label]
    excluded_count = len(genomes) - len(selected)
    if excluded_count:
        log(
            f"taxa mode: excluded {excluded_count:,} subgroup-labelled genomes; "
            "using species-directory genomes only"
        )
    if not selected:
        log("WARNING: taxa mode has no eligible genomes; producing no taxa probe sets")
        return []
    species_groups: Dict[Target, List[Genome]] = defaultdict(list)
    genus_groups: Dict[Target, List[Genome]] = defaultdict(list)
    for genome in selected:
        species_groups[make_species_target(genome)].append(genome)
        genus_groups[make_genus_target(genome)].append(genome)

    representative_manifest: List[Tuple[Target, Path, int]] = []
    first_header: Dict[str, str] = {}
    species_owner: Dict[str, Optional[Target]] = {}
    genus_owner: Dict[str, Optional[Target]] = {}

    for target in sorted(species_groups, key=lambda t: t.target_id):
        group = species_groups[target]
        representative, threshold = representative_kmers(
            [windows[g.accession].kmer_merged for g in group],
            fraction,
            small_group_all_max,
        )
        output = (
            mode_root
            / "representative_kmers"
            / f"{target_file_component(target)}.representative_kmer.fasta"
        )
        count = write_fasta(output, ((header, sequence) for sequence, header in representative.items()))
        representative_manifest.append((target, output, count))
        log(
            f"taxa representative: {target.name}: {len(group)} genomes, "
            f"threshold={threshold}, k-mers={count:,}"
        )
        for sequence, header in representative.items():
            first_header.setdefault(sequence, header)
            previous_species = species_owner.get(sequence)
            if sequence not in species_owner:
                species_owner[sequence] = target
            elif previous_species != target:
                species_owner[sequence] = None
            genus_target = make_genus_target(group[0])
            previous_genus = genus_owner.get(sequence)
            if sequence not in genus_owner:
                genus_owner[sequence] = genus_target
            elif previous_genus != genus_target:
                genus_owner[sequence] = None

    write_target_manifest(mode_root / "representative_kmers" / "manifest.tsv", representative_manifest)

    species_specific: Dict[Target, Dict[str, str]] = defaultdict(dict)
    genus_specific: Dict[Target, Dict[str, str]] = defaultdict(dict)
    for sequence, owner in species_owner.items():
        if owner is not None:
            species_specific[owner][sequence] = first_header[sequence]
    for sequence, owner in genus_owner.items():
        if owner is not None:
            genus_specific[owner][sequence] = first_header[sequence]

    # Retain genus-specific k-mers as inspectable intermediates. They are a
    # specificity criterion, not a separate genus-level probe output target.
    for target in sorted(genus_groups, key=lambda t: t.target_id):
        kmers = genus_specific.get(target, {})
        kmer_path = (
            mode_root
            / "specific_kmers"
            / "genus"
            / f"{target_file_component(target)}.genus_specific_kmer.fasta"
        )
        write_fasta(kmer_path, ((header, seq) for seq, header in kmers.items()))

    outputs: List[Tuple[Target, Path]] = []
    probe_manifest: List[Tuple[Target, Path, int]] = []
    for target in sorted(species_groups, key=lambda t: t.target_id):
        representative_sequences = set(
            sequence_dictionary(
                mode_root
                / "representative_kmers"
                / f"{target_file_component(target)}.representative_kmer.fasta"
            )
        )
        species_kmers = species_specific.get(target, {})
        species_kmer_path = (
            mode_root
            / "specific_kmers"
            / "species"
            / f"{target_file_component(target)}.species_specific_kmer.fasta"
        )
        write_fasta(
            species_kmer_path,
            ((header, seq) for seq, header in species_kmers.items()),
        )

        genus_target = make_genus_target(species_groups[target][0])
        genus_kmers = {
            sequence: header
            for sequence, header in genus_specific.get(genus_target, {}).items()
            if sequence in representative_sequences
        }
        genus_by_species_path = (
            mode_root
            / "specific_kmers"
            / "genus_by_species"
            / f"{target_file_component(target)}.genus_specific_kmer.fasta"
        )
        write_fasta(
            genus_by_species_path,
            ((header, seq) for seq, header in genus_kmers.items()),
        )
        combined_kmers = dict(genus_kmers)
        combined_kmers.update(species_kmers)
        combined_kmer_path = (
            mode_root
            / "specific_kmers"
            / "combined_by_species"
            / f"{target_file_component(target)}.taxa_specific_kmer.fasta"
        )
        write_fasta(
            combined_kmer_path,
            ((header, seq) for seq, header in combined_kmers.items()),
        )

        probe_path = (
            mode_root
            / "probe_sets"
            / "species"
            / f"{target_file_component(target)}.taxa_specific.probe.fasta"
        )
        count = intersect_individuals_then_merge(
            species_groups[target],
            windows,
            set(combined_kmers),
            mode_root / "intersections" / "species" / target_file_component(target),
            probe_path,
        )
        outputs.append((target, probe_path))
        probe_manifest.append((target, probe_path, count))

    write_target_manifest(mode_root / "probe_sets" / "manifest.tsv", probe_manifest)
    log(
        f"eMito-taxa-generate produced {len(outputs):,} species Taxa-specific "
        "probe files (no genus-level probe files)"
    )
    return outputs


def run_group_mode(
    genomes: Sequence[Genome],
    windows: Mapping[str, WindowFiles],
    mode_root: Path,
    fraction: float,
    small_group_all_max: int,
) -> List[Tuple[Target, Path]]:
    selected = [g for g in genomes if g.subgroup_label]
    if not selected:
        log("WARNING: group mode has no subgroup-labelled genomes")
        return []
    groups: Dict[Target, List[Genome]] = defaultdict(list)
    for genome in selected:
        groups[make_group_target(genome)].append(genome)

    owner: Dict[str, Optional[Target]] = {}
    first_header: Dict[str, str] = {}
    representative_manifest: List[Tuple[Target, Path, int]] = []
    for target in sorted(groups, key=lambda t: t.target_id):
        group = groups[target]
        representative, threshold = representative_kmers(
            [windows[g.accession].kmer_merged for g in group],
            fraction,
            small_group_all_max,
        )
        output = (
            mode_root
            / "representative_kmers"
            / safe_component(group[0].species_name)
            / f"{target_file_component(target)}.representative_kmer.fasta"
        )
        count = write_fasta(output, ((header, sequence) for sequence, header in representative.items()))
        representative_manifest.append((target, output, count))
        log(
            f"group representative: {target.name}: {len(group)} genomes, "
            f"threshold={threshold}, k-mers={count:,}"
        )
        for sequence, header in representative.items():
            first_header.setdefault(sequence, header)
            previous = owner.get(sequence)
            if sequence not in owner:
                owner[sequence] = target
            elif previous != target:
                owner[sequence] = None

    write_target_manifest(mode_root / "representative_kmers" / "manifest.tsv", representative_manifest)
    specific: Dict[Target, Dict[str, str]] = defaultdict(dict)
    for sequence, target in owner.items():
        if target is not None:
            specific[target][sequence] = first_header[sequence]

    outputs: List[Tuple[Target, Path]] = []
    probe_manifest: List[Tuple[Target, Path, int]] = []
    for target in sorted(groups, key=lambda t: t.target_id):
        species_component = safe_component(groups[target][0].species_name)
        kmers = specific.get(target, {})
        kmer_path = (
            mode_root
            / "specific_kmers"
            / species_component
            / f"{target_file_component(target)}.group_specific_kmer.fasta"
        )
        write_fasta(kmer_path, ((header, seq) for seq, header in kmers.items()))
        probe_path = (
            mode_root
            / "probe_sets"
            / species_component
            / f"{target_file_component(target)}.group_specific.probe.fasta"
        )
        count = intersect_individuals_then_merge(
            groups[target],
            windows,
            set(kmers),
            mode_root
            / "intersections"
            / species_component
            / target_file_component(target),
            probe_path,
        )
        outputs.append((target, probe_path))
        probe_manifest.append((target, probe_path, count))
    write_target_manifest(mode_root / "probe_sets" / "manifest.tsv", probe_manifest)
    log(f"eMito-group-generate produced {len(outputs):,} target probe files")
    return outputs


def stable_node_choice(genomes: Sequence[Genome], node_taxid: int, seed: int) -> Genome:
    preferred = [genome for genome in genomes if genome.accession.startswith("NC_")]
    pool = sorted(preferred or list(genomes), key=lambda genome: genome.accession)
    chooser = random.Random(f"{seed}:{node_taxid}")
    return pool[chooser.randrange(len(pool))]


def run_node_mode(
    genomes: Sequence[Genome],
    nodes: TaxonomyNodes,
    names_dmp: Path,
    mode_root: Path,
    rank: str,
    window_length: int,
    node_step: int,
    seed: int,
) -> List[Tuple[Target, Path]]:
    by_node: Dict[int, List[Genome]] = defaultdict(list)
    for genome in genomes:
        node_taxid = ancestor_at_rank(genome.taxid, rank, nodes)
        by_node[node_taxid].append(genome)
    names = load_names(names_dmp, set(by_node))
    outputs: List[Tuple[Target, Path]] = []
    manifest_rows: List[Tuple[Target, Path, int]] = []
    selection_path = mode_root / "selected_genomes.tsv"
    selection_path.parent.mkdir(parents=True, exist_ok=True)
    with selection_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("node_rank", "node_taxid", "node_name", "selected_accession", "nc_preferred"))
        for node_taxid in sorted(by_node):
            chosen = stable_node_choice(by_node[node_taxid], node_taxid, seed)
            target = Target(
                target_id=f"{safe_component(rank)}_taxid_{node_taxid}",
                rank=rank,
                name=names[node_taxid],
                taxid=node_taxid,
            )
            output = (
                mode_root
                / "probe_sets"
                / f"{target_file_component(target)}.node_tiling.probe.fasta"
            )
            _, sequence = read_single_fasta(chosen.organized_fasta)
            count = write_fasta(
                output,
                (
                    (forward_header(chosen.accession, start, end, "node_probe"), window)
                    for start, end, window in iter_windows(sequence, window_length, node_step)
                ),
            )
            writer.writerow(
                (rank, node_taxid, names[node_taxid], chosen.accession, chosen.accession.startswith("NC_"))
            )
            outputs.append((target, output))
            manifest_rows.append((target, output, count))
    write_target_manifest(mode_root / "probe_sets" / "manifest.tsv", manifest_rows)
    log(f"eMito-node-generate selected {len(outputs):,} {rank} node genomes")
    return outputs


def gc_percent(sequence: str) -> float:
    if not sequence:
        return 0.0
    return round(100.0 * (sequence.count("G") + sequence.count("C")) / len(sequence), 4)


def dust_complexity(sequence: str) -> float:
    k = 3
    if len(sequence) <= k:
        return 0.0
    counts = Counter(sequence[index : index + k] for index in range(len(sequence) - k + 1))
    pairs = sum(count * (count - 1) / 2 for count in counts.values())
    return round(pairs / (len(sequence) - 3), 4)


def dimer_scores(records: Sequence[Tuple[str, str]], k: int) -> List[float]:
    """Calculate the eProbe-compatible ``DimerCalculatorFast`` score.

    Reverse-complement k-mers from the complete post-GC/complexity probe pool
    are indexed. Each probe's forward k-mers query that index, its own reverse-
    complement contribution is subtracted, and the remaining inter-probe signal
    is normalized by both k-mer count and pool size.
    """
    if not records:
        return []

    sequences = [sequence.upper() for _, sequence in records]
    reverse_kmer_frequency: Counter[str] = Counter()
    for sequence in sequences:
        rc = reverse_complement(sequence)
        if len(rc) >= k:
            reverse_kmer_frequency.update(
                rc[index : index + k] for index in range(len(rc) - k + 1)
            )

    probe_count = len(sequences)
    scores: List[float] = []
    for sequence in sequences:
        if len(sequence) < k:
            scores.append(0.0)
            continue

        rc = reverse_complement(sequence)
        self_reverse_frequency = Counter(
            rc[index : index + k] for index in range(len(rc) - k + 1)
        )
        total_frequency = 0
        kmer_count = len(sequence) - k + 1
        for index in range(kmer_count):
            kmer = sequence[index : index + k]
            other_frequency = max(
                reverse_kmer_frequency.get(kmer, 0)
                - self_reverse_frequency.get(kmer, 0),
                0,
            )
            total_frequency += other_frequency
        scores.append(round(total_frequency / (kmer_count * probe_count) * 10000, 2))
    return scores


def eprobe_dimer_cutoff(scores: Sequence[float], threshold: float) -> float:
    """Resolve eProbe's quantile/absolute dual-mode dimer threshold."""
    if threshold <= 0 or not scores:
        return math.inf
    if threshold < 1.0:
        ordered = sorted(scores)
        index = min(int(threshold * len(ordered)), len(ordered) - 1)
        return ordered[index]
    return threshold


def run_access_task(task: AccessTask) -> Dict[str, int]:
    input_path = Path(task.input_fasta)
    records = [
        (header, sequence)
        for header, sequence in iter_fasta(input_path)
        if DNA_RE.fullmatch(sequence)
    ]
    # Match eProbe's order: absolute GC/DUST filters first, then construct the
    # dimer pool only from probes that passed those two filters.
    stage_one: List[Tuple[str, str, float, float]] = []
    for header, sequence in records:
        gc = gc_percent(sequence)
        complexity = dust_complexity(sequence)
        if (
            task.parameters.gc_min <= gc <= task.parameters.gc_max
            and task.parameters.complexity_min <= complexity <= task.parameters.complexity_max
        ):
            stage_one.append((header, sequence, gc, complexity))

    if len(stage_one) > 1 and task.parameters.dimer_threshold > 0:
        dimer_input = [(header, sequence) for header, sequence, _, _ in stage_one]
        scores = dimer_scores(dimer_input, task.parameters.dimer_k)
        cutoff = eprobe_dimer_cutoff(scores, task.parameters.dimer_threshold)
    else:
        scores = [0.0] * len(stage_one)
        cutoff = math.inf

    assessed = [
        (header, sequence, gc, complexity, dimer_score)
        for (header, sequence, gc, complexity), dimer_score in zip(stage_one, scores)
        if dimer_score <= cutoff
    ]
    assessed.sort(key=lambda row: (row[4], row[3], row[2]))
    deduplicated: List[Tuple[str, str, float, float, float]] = []
    seen: Set[str] = set()
    for row in assessed:
        if row[1] not in seen:
            seen.add(row[1])
            deduplicated.append(row)

    output_tsv = Path(task.output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("source_file", "id", "full_header", "gc_content", "Sequence", "Complexity", "Dimer_Score")
        )
        for header, sequence, gc, complexity, dimer_score in deduplicated:
            writer.writerow(
                (input_path.name, header.split()[0], header, gc, sequence, complexity, dimer_score)
            )
    write_fasta(
        Path(task.output_fasta), ((header, sequence) for header, sequence, *_ in deduplicated)
    )
    return {"input": len(records), "filtered": len(deduplicated)}


def run_collapse_task(task: CollapseTask) -> Dict[str, int]:
    records = list(iter_fasta(Path(task.input_fasta)))
    grouped: Dict[str, List[Tuple[int, int, int, str, str]]] = defaultdict(list)
    for order, (header, sequence) in enumerate(records):
        match = COORDINATE_RE.match(header)
        if match is None:
            raise PipelineError(f"Cannot parse accession/coordinates from header: {header}")
        normalized_match = NORMALIZED_COORDINATE_RE.search(header)
        if normalized_match is not None:
            start = int(normalized_match.group("start"))
            end = int(normalized_match.group("end"))
        else:
            start = int(match.group("start"))
            end = int(match.group("end"))
        grouped[match.group("accession")].append(
            (start, end, order, header, sequence)
        )
    retained: List[Tuple[str, str]] = []
    for accession in sorted(grouped):
        intervals = sorted(grouped[accession], key=lambda row: (row[0], row[1], row[2]))
        last_retained_end = -1
        for start, end, order, header, sequence in intervals:
            if start > last_retained_end:
                retained.append((header, sequence))
                last_retained_end = end
    write_fasta(Path(task.output_fasta), retained)
    return {"input": len(records), "retained": len(retained)}


def optional_mode_pipeline(
    mode: str,
    inputs: Sequence[Tuple[Target, Path]],
    optional_root: Path,
    run_access: bool,
    run_collapse: bool,
    access_parameters: AccessParameters,
    threads: int,
) -> List[Tuple[Target, Path]]:
    current = list(inputs)
    if run_access:
        tasks: List[AccessTask] = []
        output_map: Dict[Path, Tuple[Target, Path]] = {}
        for target, input_path in current:
            relative = input_path.name
            output = optional_root / mode / "access" / relative.replace(".fasta", ".assessed.fasta")
            table = output.with_suffix(".tsv")
            tasks.append(
                AccessTask(
                    input_fasta=str(input_path),
                    output_fasta=str(output),
                    output_tsv=str(table),
                    parameters=access_parameters,
                )
            )
            output_map[input_path] = (target, output)
        execute_process_tasks(run_access_task, tasks, threads, f"eMito-access {mode}")
        current = [output_map[input_path] for _, input_path in current]
    if run_collapse:
        tasks = []
        output_map = {}
        for target, input_path in current:
            output = (
                optional_root
                / mode
                / "collapse"
                / input_path.name.replace(".fasta", ".collapsed.fasta")
            )
            tasks.append(CollapseTask(input_fasta=str(input_path), output_fasta=str(output)))
            output_map[input_path] = (target, output)
        execute_process_tasks(run_collapse_task, tasks, threads, f"eMito-collapse {mode}")
        current = [output_map[input_path] for _, input_path in current]
    return current


def merge_final_outputs(
    mode_outputs: Mapping[str, Sequence[Tuple[Target, Path]]],
    final_root: Path,
    deduplicate: bool = True,
) -> Tuple[Path, int]:
    final_root.mkdir(parents=True, exist_ok=True)
    final_path = final_root / "final_probe_set.fasta"
    manifest_path = final_root / "input_probe_files.tsv"
    seen: Set[str] = set()
    count = 0
    with final_path.open("wt", encoding="ascii") as output, manifest_path.open(
        "wt", encoding="utf-8", newline=""
    ) as manifest:
        writer = csv.writer(manifest, delimiter="\t", lineterminator="\n")
        writer.writerow(("mode", "target_id", "rank", "target_name", "input_probe_file"))
        for mode in ("taxa", "group", "node"):
            for target, path in mode_outputs.get(mode, []):
                writer.writerow((mode, target.target_id, target.rank, target.name, path))
                for header, sequence in iter_fasta(path):
                    if not DNA_RE.fullmatch(sequence):
                        continue
                    if deduplicate:
                        if sequence in seen:
                            continue
                        seen.add(sequence)
                    output.write(f">{header}\n{sequence}\n")
                    count += 1
    return final_path, count


def write_genome_manifest(path: Path, genomes: Sequence[Genome]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "accession_id",
                "genus_name",
                "genus_taxid",
                "species_name",
                "species_taxid",
                "actual_taxid",
                "subgroup_label",
                "sequence_length",
                "organized_fasta",
            )
        )
        for genome in genomes:
            writer.writerow(
                (
                    genome.accession,
                    genome.genus_name,
                    genome.genus_taxid,
                    genome.species_name,
                    genome.species_taxid,
                    genome.taxid,
                    genome.subgroup_label,
                    genome.sequence_length,
                    genome.organized_fasta,
                )
            )


def prepare_output_root(path: Path, force: bool) -> None:
    if path.exists():
        if not force:
            raise PipelineError(
                f"Output root already exists: {path}. Use a new --output-root or --force."
            )
        resolved = path.resolve()
        if str(resolved) in {"/", str(Path.home().resolve())} or len(resolved.parts) < 4:
            raise PipelineError(f"Refusing to remove unsafe output path: {resolved}")
        shutil.rmtree(resolved)
    path.mkdir(parents=True)


def configuration_dict(args: argparse.Namespace, modes: Set[str]) -> Dict[str, object]:
    result: Dict[str, object] = {}
    for key, value in vars(args).items():
        result[key] = str(value) if isinstance(value, Path) else value
    result["resolved_modes"] = sorted(modes)
    result["pipeline_version"] = PIPELINE_VERSION
    return result


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        log(f"eMito pipeline version: {PIPELINE_VERSION}")
        modes = validate_args(args)
        mafft_executable = resolve_executable(args.mafft, "within-group alignment")
        log(f"MAFFT executable: {mafft_executable}")
        require_file(args.metadata, "metadata")
        require_directory(args.fasta_dir, "input FASTA directory")
        require_file(args.names_dmp, "names.dmp")
        require_file(args.nodes_dmp, "nodes.dmp")
        rows = load_metadata(args.metadata)
        nodes = load_nodes(args.nodes_dmp)

        organized_root = args.output_root / "01_organized_genomes"
        all_genomes = build_genomes(
            rows, args.fasta_dir, nodes, args.names_dmp, organized_root
        )
        log(
            f"Validated {len(all_genomes):,} genomes across "
            f"{len({g.species_taxid for g in all_genomes}):,} species and "
            f"{len({g.genus_taxid for g in all_genomes}):,} genera"
        )
        genomes, length_qc = filter_genomes_by_length(
            all_genomes, args.min_genome_length
        )
        excluded_short = len(all_genomes) - len(genomes)
        log(
            f"Input length QC: retained {len(genomes):,}; excluded "
            f"{excluded_short:,} genomes shorter than {args.min_genome_length:,} bp"
        )
        if not genomes:
            raise PipelineError(
                f"All input genomes are shorter than --min-genome-length "
                f"{args.min_genome_length:,}"
            )
        if args.validate_only:
            log("Validation completed; no output was generated")
            return 0

        prepare_output_root(args.output_root, args.force)
        (args.output_root / "00_config.json").write_text(
            json.dumps(configuration_dict(args, modes), indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        write_genome_length_qc(args.output_root, length_qc)
        organize_genomes(genomes, args.organize_method)
        write_genome_manifest(args.output_root / "01_organized_genomes" / "manifest.tsv", genomes)
        log(f"Organized {len(genomes):,} genome FASTA files")

        tasks, window_files = build_window_tasks(
            genomes,
            args.output_root,
            args.window_length,
            args.probe_step,
            args.kmer_step,
            mafft_executable,
            args.circular_anchor_length,
            args.circular_min_anchor_length,
            args.circular_min_anchor_hits,
        )
        write_alignment_group_manifest(
            args.output_root / "02_alignments" / "manifest.tsv", tasks
        )
        execute_process_tasks(
            run_window_task,
            tasks,
            args.threads,
            "group alignment and aligned-window generation",
        )

        generated: Dict[str, List[Tuple[Target, Path]]] = {}
        if "taxa" in modes:
            generated["taxa"] = run_taxa_mode(
                genomes,
                window_files,
                args.output_root / "03_eMito_taxa_generate",
                args.representative_fraction,
                args.small_group_all_max,
            )
        if "group" in modes:
            generated["group"] = run_group_mode(
                genomes,
                window_files,
                args.output_root / "04_eMito_group_generate",
                args.representative_fraction,
                args.small_group_all_max,
            )
        if "node" in modes:
            generated["node"] = run_node_mode(
                genomes,
                nodes,
                args.names_dmp,
                args.output_root / "05_eMito_node_generate",
                args.node_rank,
                args.window_length,
                args.node_step,
                args.random_seed,
            )

        access_parameters = AccessParameters(
            gc_min=args.gc_min,
            gc_max=args.gc_max,
            complexity_min=args.complexity_min,
            complexity_max=args.complexity_max,
            dimer_k=args.dimer_k,
            dimer_threshold=args.dimer_threshold,
        )
        terminal: Dict[str, List[Tuple[Target, Path]]] = {}
        mode_options = {
            "taxa": (args.taxa_access, args.taxa_collapse),
            "group": (args.group_access, args.group_collapse),
            "node": (args.node_access, args.node_collapse),
        }
        for mode, outputs in generated.items():
            run_access, run_collapse = mode_options[mode]
            terminal[mode] = optional_mode_pipeline(
                mode,
                outputs,
                args.output_root / "06_optional_modes",
                run_access,
                run_collapse,
                access_parameters,
                args.threads,
            )

        final_path, count = merge_final_outputs(
            terminal,
            args.output_root / "07_final_probe_set",
            deduplicate=args.final_dedup,
        )
        if args.final_dedup:
            log(f"Final merged and sequence-deduplicated probe set: {final_path}")
        else:
            log(f"Final merged probe set (deduplication disabled): {final_path}")
        log(f"Final probe count: {count:,}")
        return 0
    except (PipelineError, OSError, UnicodeError, ValueError) as exc:
        log(f"ERROR: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
