#!/usr/bin/env python3
"""
02_generate_52bp_kmers.py

Generate k-mers from the cleaned mitochondrial reference panel.

This script follows the eMito design logic:
  - k-mer length: 52 bp by default
  - step size: 1 bp by default
  - only A/C/G/T k-mers are emitted by default

Example
-------
python scripts/design/02_generate_52bp_kmers.py \
  --fasta results/design/01_reference_panel/reference_panel.clean.fasta \
  --metadata results/design/01_reference_panel/reference_panel.clean.tsv \
  --k 52 \
  --step 1 \
  --out results/design/02_kmers/reference_panel.52bp_kmers.tsv.gz
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path
from typing import Dict, Iterator, List, Tuple


def open_text(path, mode="rt"):
    path = str(path)
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode, newline="" if "t" in mode else None)


def parse_fasta(path) -> Iterator[Tuple[str, str, str]]:
    name, desc, chunks = None, "", []
    with open_text(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, desc, "".join(chunks).upper().replace("U", "T")
                desc = line[1:].strip()
                name = desc.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, desc, "".join(chunks).upper().replace("U", "T")


def read_metadata(path) -> Dict[str, dict]:
    if not path:
        return {}
    with open_text(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {row["sequence_id"]: row for row in reader if row.get("sequence_id")}


def seq_id_from_header(name: str) -> str:
    return name.split("|")[0]


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def valid_kmer(seq: str, allow_ambiguous: bool) -> bool:
    if allow_ambiguous:
        return bool(re.fullmatch(r"[ACGTN]+", seq))
    return bool(re.fullmatch(r"[ACGT]+", seq))


def main() -> None:
    p = argparse.ArgumentParser(description="Generate 52 bp k-mers from mitochondrial reference sequences.")
    p.add_argument("--fasta", required=True, help="Clean reference panel FASTA.")
    p.add_argument("--metadata", default=None, help="Clean reference metadata TSV from 01_prepare_reference_panel.py.")
    p.add_argument("--k", type=int, default=52, help="k-mer length.")
    p.add_argument("--step", type=int, default=1, help="Sliding-window step size.")
    p.add_argument("--include-rc", action="store_true", help="Also emit reverse-complement k-mers.")
    p.add_argument("--allow-ambiguous", action="store_true", help="Allow k-mers containing N.")
    p.add_argument("--out", required=True, help="Output TSV or TSV.GZ.")
    args = p.parse_args()

    meta = read_metadata(args.metadata)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "kmer", "sequence_id", "start", "end", "strand",
        "family", "genus", "species", "group", "haplogroup", "population"
    ]

    n_records = 0
    n_kmers = 0

    with open_text(args.out, "wt") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        for header_id, desc, seq in parse_fasta(args.fasta):
            sequence_id = seq_id_from_header(header_id)
            row = meta.get(sequence_id, {})
            n_records += 1

            if len(seq) < args.k:
                continue

            for start in range(0, len(seq) - args.k + 1, args.step):
                kmer = seq[start:start + args.k].upper()
                if not valid_kmer(kmer, args.allow_ambiguous):
                    continue
                writer.writerow({
                    "kmer": kmer,
                    "sequence_id": sequence_id,
                    "start": start + 1,
                    "end": start + args.k,
                    "strand": "+",
                    "family": row.get("family", ""),
                    "genus": row.get("genus", ""),
                    "species": row.get("species", ""),
                    "group": row.get("group", ""),
                    "haplogroup": row.get("haplogroup", ""),
                    "population": row.get("population", ""),
                })
                n_kmers += 1

                if args.include_rc:
                    writer.writerow({
                        "kmer": revcomp(kmer),
                        "sequence_id": sequence_id,
                        "start": start + 1,
                        "end": start + args.k,
                        "strand": "-",
                        "family": row.get("family", ""),
                        "genus": row.get("genus", ""),
                        "species": row.get("species", ""),
                        "group": row.get("group", ""),
                        "haplogroup": row.get("haplogroup", ""),
                        "population": row.get("population", ""),
                    })
                    n_kmers += 1

    print(f"[INFO] records processed: {n_records}", file=sys.stderr)
    print(f"[INFO] kmers written    : {n_kmers}", file=sys.stderr)
    print(f"[INFO] output           : {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
