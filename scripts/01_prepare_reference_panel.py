#!/usr/bin/env python3
"""
01_prepare_reference_panel.py

Prepare a clean mitochondrial reference panel for eMito probe design.

Typical use
-----------
python scripts/design/01_prepare_reference_panel.py \
  --metadata metadata/reference_panel_raw.tsv \
  --fasta-dir data/reference_mitogenomes/ \
  --out-fasta results/design/01_reference_panel/reference_panel.clean.fasta \
  --out-metadata results/design/01_reference_panel/reference_panel.clean.tsv

Expected metadata columns
-------------------------
Recommended columns:
  sequence_id, accession, fasta_path, family, genus, species,
  group, haplogroup, population, is_representative,
  used_for_kmer_screening, used_for_probe_tiling

Only sequence_id/accession and taxonomic columns are essential. If fasta_path is
present, it is used directly. Otherwise the script searches --fasta-dir for a
file whose basename starts with sequence_id/accession.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple


def open_text(path: str | Path, mode: str = "rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, newline="" if "t" in mode else None)


def sanitize_id(value: str) -> str:
    value = str(value).strip()
    value = re.sub(r"\s+", "_", value)
    value = re.sub(r"[^A-Za-z0-9_.|:+-]+", "_", value)
    return value.strip("_") or "unknown"


def read_table(path: str | Path) -> List[dict]:
    with open_text(path, "rt") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delimiter = "\t" if sample.count("\t") >= sample.count(",") else ","
        reader = csv.DictReader(fh, delimiter=delimiter)
        rows = []
        for row in reader:
            clean = {str(k).strip(): ("" if v is None else str(v).strip()) for k, v in row.items()}
            if any(clean.values()):
                rows.append(clean)
        return rows


def write_table(path: str | Path, rows: List[dict], fieldnames: List[str]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def parse_fasta(path: str | Path) -> Iterator[Tuple[str, str, str]]:
    name = None
    desc = ""
    seq_chunks: List[str] = []
    with open_text(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, desc, "".join(seq_chunks)
                desc = line[1:].strip()
                name = desc.split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if name is not None:
            yield name, desc, "".join(seq_chunks)


def write_fasta_record(handle, name: str, seq: str, width: int = 80) -> None:
    handle.write(f">{name}\n")
    for i in range(0, len(seq), width):
        handle.write(seq[i:i + width] + "\n")


def normalize_sequence(seq: str, replace_non_acgtn: bool = True) -> str:
    seq = re.sub(r"\s+", "", seq).upper().replace("U", "T")
    if replace_non_acgtn:
        seq = re.sub(r"[^ACGTN]", "N", seq)
    return seq


def first_existing_column(row: dict, candidates: Iterable[str]) -> Optional[str]:
    lower = {k.lower(): k for k in row.keys()}
    for c in candidates:
        if c.lower() in lower and row.get(lower[c.lower()], "") != "":
            return lower[c.lower()]
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    return None


def auto_sequence_id(row: dict, id_column: Optional[str]) -> str:
    if id_column and id_column in row and row[id_column]:
        return sanitize_id(row[id_column])
    col = first_existing_column(row, ["sequence_id", "seq_id", "accession", "id", "name"])
    if col and row.get(col):
        return sanitize_id(row[col])
    raise ValueError(f"Could not infer sequence id from row: {row}")


def build_fasta_index(fasta_paths: List[Path]) -> Dict[str, Tuple[Path, str]]:
    """Map first-token FASTA id to (path, full_header)."""
    idx: Dict[str, Tuple[Path, str]] = {}
    for fp in fasta_paths:
        for rec_id, desc, _seq in parse_fasta(fp):
            idx.setdefault(rec_id, (fp, desc))
    return idx


def discover_fasta_files(fasta_dir: str | Path, extensions: List[str]) -> List[Path]:
    fasta_dir = Path(fasta_dir)
    files: List[Path] = []
    for ext in extensions:
        files.extend(fasta_dir.rglob(f"*{ext}"))
    return sorted(set(files))


def find_fasta_for_row(row: dict, seq_id: str, fasta_dir: Optional[Path], extensions: List[str], path_column: str) -> Optional[Path]:
    if path_column in row and row[path_column]:
        p = Path(row[path_column])
        if not p.is_absolute() and fasta_dir is not None:
            p2 = fasta_dir / p
            if p2.exists():
                return p2
        if p.exists():
            return p

    if fasta_dir is None:
        return None

    patterns = []
    for key in ["sequence_id", "accession", "id", "name"]:
        if key in row and row[key]:
            patterns.append(sanitize_id(row[key]))
            patterns.append(row[key].strip())
    patterns.append(seq_id)

    for pat in patterns:
        for ext in extensions:
            hits = list(fasta_dir.rglob(f"{pat}*{ext}"))
            if hits:
                return sorted(hits)[0]
    return None


def metadata_value(row: dict, names: List[str]) -> str:
    col = first_existing_column(row, names)
    return row.get(col, "") if col else ""


def make_output_header(seq_id: str, row: dict) -> str:
    fields = {
        "family": metadata_value(row, ["family"]),
        "genus": metadata_value(row, ["genus"]),
        "species": metadata_value(row, ["species"]),
        "group": metadata_value(row, ["group", "taxon_group"]),
        "haplogroup": metadata_value(row, ["haplogroup", "lineage"]),
        "population": metadata_value(row, ["population", "pop"]),
    }
    tag = "|".join([seq_id] + [f"{k}={sanitize_id(v)}" for k, v in fields.items() if v])
    return tag


def main() -> None:
    p = argparse.ArgumentParser(
        description="Prepare a clean combined mitochondrial reference FASTA and metadata table."
    )
    p.add_argument("--metadata", required=True, help="Reference metadata table, TSV or CSV.")
    p.add_argument("--fasta-dir", default=None, help="Directory containing individual FASTA files.")
    p.add_argument("--combined-fasta", nargs="*", default=None, help="Optional already-combined FASTA file(s).")
    p.add_argument("--out-fasta", required=True, help="Output combined clean FASTA.")
    p.add_argument("--out-metadata", required=True, help="Output clean metadata TSV.")
    p.add_argument("--id-column", default=None, help="Metadata column to use as sequence_id. Auto-detected if omitted.")
    p.add_argument("--path-column", default="fasta_path", help="Metadata column containing FASTA path.")
    p.add_argument("--extensions", default=".fa,.fasta,.fna,.fas,.fa.gz,.fasta.gz,.fna.gz",
                   help="Comma-separated FASTA extensions for --fasta-dir search.")
    p.add_argument("--min-length", type=int, default=1000, help="Discard sequences shorter than this length.")
    p.add_argument("--max-n-fraction", type=float, default=1.0, help="Discard sequences with N fraction above this value.")
    p.add_argument("--replace-non-acgtn-with-n", action="store_true", default=True,
                   help="Replace non-ACGTN characters with N. Enabled by default.")
    args = p.parse_args()

    rows = read_table(args.metadata)
    if not rows:
        sys.exit("[ERROR] metadata table is empty")

    extensions = [x.strip() for x in args.extensions.split(",") if x.strip()]
    fasta_dir = Path(args.fasta_dir) if args.fasta_dir else None

    combined_records: Dict[str, Tuple[str, str]] = {}
    if args.combined_fasta:
        for fp in args.combined_fasta:
            for rec_id, desc, seq in parse_fasta(fp):
                combined_records.setdefault(rec_id, (desc, seq))

    fasta_files = discover_fasta_files(fasta_dir, extensions) if fasta_dir else []

    out_rows: List[dict] = []
    kept = 0
    skipped = 0

    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.out_fasta, "wt") as out_fa:
        for row in rows:
            try:
                seq_id = auto_sequence_id(row, args.id_column)
            except ValueError as exc:
                print(f"[WARN] {exc}", file=sys.stderr)
                skipped += 1
                continue

            seq = None

            if seq_id in combined_records:
                seq = combined_records[seq_id][1]
            else:
                # Try accession or raw id against combined FASTA ids.
                for key in ["accession", "sequence_id", "id", "name"]:
                    val = row.get(key, "")
                    if val and val in combined_records:
                        seq = combined_records[val][1]
                        break

            if seq is None:
                fp = find_fasta_for_row(row, seq_id, fasta_dir, extensions, args.path_column)
                if fp is None:
                    print(f"[WARN] no FASTA found for {seq_id}", file=sys.stderr)
                    skipped += 1
                    continue
                records = list(parse_fasta(fp))
                if not records:
                    print(f"[WARN] empty FASTA: {fp}", file=sys.stderr)
                    skipped += 1
                    continue
                if len(records) > 1:
                    # Prefer a record whose id matches the metadata; otherwise use the longest.
                    match = [r for r in records if r[0] == seq_id or r[0] == row.get("accession", "")]
                    rec = match[0] if match else max(records, key=lambda x: len(x[2]))
                else:
                    rec = records[0]
                seq = rec[2]

            seq = normalize_sequence(seq, args.replace_non_acgtn_with_n)
            n_fraction = seq.count("N") / max(1, len(seq))
            if len(seq) < args.min_length or n_fraction > args.max_n_fraction:
                print(f"[WARN] skipped {seq_id}: length={len(seq)}, N_fraction={n_fraction:.4f}", file=sys.stderr)
                skipped += 1
                continue

            header = make_output_header(seq_id, row)
            write_fasta_record(out_fa, header, seq)

            clean = {
                "sequence_id": seq_id,
                "accession": metadata_value(row, ["accession"]),
                "family": metadata_value(row, ["family"]),
                "genus": metadata_value(row, ["genus"]),
                "species": metadata_value(row, ["species"]),
                "group": metadata_value(row, ["group", "taxon_group"]),
                "haplogroup": metadata_value(row, ["haplogroup", "lineage"]),
                "population": metadata_value(row, ["population", "pop"]),
                "is_representative": metadata_value(row, ["is_representative", "representative"]),
                "used_for_kmer_screening": metadata_value(row, ["used_for_kmer_screening", "kmer_screening"]),
                "used_for_probe_tiling": metadata_value(row, ["used_for_probe_tiling", "probe_tiling"]),
                "used_for_family_tiling": metadata_value(row, ["used_for_family_tiling", "family_tiling"]),
                "used_for_haplogroup_tiling": metadata_value(row, ["used_for_haplogroup_tiling", "haplogroup_tiling"]),
                "sequence_length": str(len(seq)),
                "n_fraction": f"{n_fraction:.6f}",
            }
            out_rows.append(clean)
            kept += 1

    fieldnames = [
        "sequence_id", "accession", "family", "genus", "species",
        "group", "haplogroup", "population",
        "is_representative", "used_for_kmer_screening",
        "used_for_probe_tiling", "used_for_family_tiling", "used_for_haplogroup_tiling",
        "sequence_length", "n_fraction",
    ]
    write_table(args.out_metadata, out_rows, fieldnames)

    print(f"[INFO] kept references   : {kept}", file=sys.stderr)
    print(f"[INFO] skipped references: {skipped}", file=sys.stderr)
    print(f"[INFO] out FASTA         : {args.out_fasta}", file=sys.stderr)
    print(f"[INFO] out metadata      : {args.out_metadata}", file=sys.stderr)


if __name__ == "__main__":
    main()
