#!/usr/bin/env python3
"""
03_generate_52bp_probe_tiles.py

Generate 52 bp probe tiles from representative mitochondrial genomes.

This follows the eMito design logic:
  - probe length: 52 bp by default
  - step size: 5 bp by default
  - one high-quality representative mitogenome per target taxon is recommended

Example
-------
python scripts/design/03_generate_52bp_probe_tiles.py \
  --fasta results/design/01_reference_panel/reference_panel.clean.fasta \
  --metadata results/design/01_reference_panel/reference_panel.clean.tsv \
  --representative-column used_for_probe_tiling \
  --probe-len 52 \
  --step 5 \
  --out-fasta results/design/03_probe_tiles/representative_tiles.fasta \
  --out-metadata results/design/03_probe_tiles/representative_tiles.metadata.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path
from typing import Dict, Iterator, List, Tuple


TRUE_VALUES = {"1", "true", "t", "yes", "y", "representative", "use", "used"}


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


def write_fasta_record(handle, name: str, seq: str, width: int = 80) -> None:
    handle.write(f">{name}\n")
    for i in range(0, len(seq), width):
        handle.write(seq[i:i + width] + "\n")


def read_metadata(path) -> Dict[str, dict]:
    if not path:
        return {}
    with open_text(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {row["sequence_id"]: row for row in reader if row.get("sequence_id")}


def seq_id_from_header(name: str) -> str:
    return name.split("|")[0]


def is_true(value: str) -> bool:
    return str(value).strip().lower() in TRUE_VALUES


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def gc_percent(seq: str) -> float:
    acgt = sum(seq.count(x) for x in "ACGT")
    return 0.0 if acgt == 0 else 100.0 * (seq.count("G") + seq.count("C")) / acgt


def main() -> None:
    p = argparse.ArgumentParser(description="Generate 52 bp probe tiles from representative mitogenomes.")
    p.add_argument("--fasta", required=True, help="Clean reference panel FASTA.")
    p.add_argument("--metadata", default=None, help="Clean metadata TSV.")
    p.add_argument("--representative-column", default="used_for_probe_tiling",
                   help="Metadata column used to select representative sequences. If absent, all sequences are used.")
    p.add_argument("--representative-values", default="1,true,yes,y,representative,use,used",
                   help="Comma-separated values treated as true for representative selection.")
    p.add_argument("--probe-len", type=int, default=52, help="Probe length.")
    p.add_argument("--step", type=int, default=5, help="Tiling step size.")
    p.add_argument("--max-n", type=int, default=0, help="Maximum number of Ns allowed in a probe tile.")
    p.add_argument("--include-rc", action="store_true", help="Also output reverse-complement tiles.")
    p.add_argument("--component", default="representative_tile", help="Component label written to metadata.")
    p.add_argument("--out-fasta", required=True, help="Output probe FASTA.")
    p.add_argument("--out-metadata", required=True, help="Output probe metadata TSV.")
    args = p.parse_args()

    global TRUE_VALUES
    TRUE_VALUES = {x.strip().lower() for x in args.representative_values.split(",") if x.strip()}

    meta = read_metadata(args.metadata)
    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_metadata).parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "probe_id", "sequence", "length", "gc_content",
        "reference_id", "start", "end", "strand",
        "family", "genus", "species", "group", "haplogroup", "population",
        "component"
    ]

    n_refs = 0
    n_tiles = 0

    with open_text(args.out_fasta, "wt") as out_fa, open_text(args.out_metadata, "wt") as out_meta:
        writer = csv.DictWriter(out_meta, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        for header_id, desc, seq in parse_fasta(args.fasta):
            ref_id = seq_id_from_header(header_id)
            row = meta.get(ref_id, {})
            if meta and args.representative_column in row and not is_true(row.get(args.representative_column, "")):
                continue

            n_refs += 1
            for start0 in range(0, len(seq) - args.probe_len + 1, args.step):
                probe_seq = seq[start0:start0 + args.probe_len].upper()
                if probe_seq.count("N") > args.max_n or not re.fullmatch(r"[ACGTN]+", probe_seq):
                    continue

                probe_id = f"tile_{n_tiles + 1:09d}|ref={ref_id}|start={start0 + 1}|end={start0 + args.probe_len}|strand=+"
                write_fasta_record(out_fa, probe_id, probe_seq)
                writer.writerow({
                    "probe_id": probe_id,
                    "sequence": probe_seq,
                    "length": len(probe_seq),
                    "gc_content": f"{gc_percent(probe_seq):.4f}",
                    "reference_id": ref_id,
                    "start": start0 + 1,
                    "end": start0 + args.probe_len,
                    "strand": "+",
                    "family": row.get("family", ""),
                    "genus": row.get("genus", ""),
                    "species": row.get("species", ""),
                    "group": row.get("group", ""),
                    "haplogroup": row.get("haplogroup", ""),
                    "population": row.get("population", ""),
                    "component": args.component,
                })
                n_tiles += 1

                if args.include_rc:
                    rc_seq = revcomp(probe_seq)
                    probe_id = f"tile_{n_tiles + 1:09d}|ref={ref_id}|start={start0 + 1}|end={start0 + args.probe_len}|strand=-"
                    write_fasta_record(out_fa, probe_id, rc_seq)
                    writer.writerow({
                        "probe_id": probe_id,
                        "sequence": rc_seq,
                        "length": len(rc_seq),
                        "gc_content": f"{gc_percent(rc_seq):.4f}",
                        "reference_id": ref_id,
                        "start": start0 + 1,
                        "end": start0 + args.probe_len,
                        "strand": "-",
                        "family": row.get("family", ""),
                        "genus": row.get("genus", ""),
                        "species": row.get("species", ""),
                        "group": row.get("group", ""),
                        "haplogroup": row.get("haplogroup", ""),
                        "population": row.get("population", ""),
                        "component": args.component,
                    })
                    n_tiles += 1

    print(f"[INFO] representative references used: {n_refs}", file=sys.stderr)
    print(f"[INFO] probe tiles written         : {n_tiles}", file=sys.stderr)
    print(f"[INFO] output FASTA                : {args.out_fasta}", file=sys.stderr)
    print(f"[INFO] output metadata             : {args.out_metadata}", file=sys.stderr)


if __name__ == "__main__":
    main()
