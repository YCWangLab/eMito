#!/usr/bin/env python3
"""
07_generate_haplogroup_family_tiling.py

Generate supplementary unfiltered tiling probes at family and haplogroup levels.

This implements the eMito supplementary tiling logic:
  - one representative species/genome per family;
  - one mitochondrial genome per haplogroup/lineage for Homo and domestic taxa;
  - 52 bp probes with 5 bp step by default;
  - no GC/DUST/dimer filtering by default, because this step is meant to retain
    broad phylogenetic representation and capture sensitivity.

Example
-------
python scripts/design/07_generate_haplogroup_family_tiling.py \
  --fasta results/design/01_reference_panel/reference_panel.clean.fasta \
  --metadata results/design/01_reference_panel/reference_panel.clean.tsv \
  --levels family,haplogroup \
  --probe-len 52 \
  --step 5 \
  --out-fasta results/design/07_supplementary_tiling/family_haplogroup_tiling.fasta \
  --out-metadata results/design/07_supplementary_tiling/family_haplogroup_tiling.metadata.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
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


def read_metadata(path: str) -> Dict[str, dict]:
    with open_text(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {row["sequence_id"]: row for row in reader if row.get("sequence_id")}


def seq_id_from_header(name: str) -> str:
    return name.split("|")[0]


def is_true(value: str) -> bool:
    return str(value or "").strip().lower() in TRUE_VALUES


def gc_percent(seq: str) -> float:
    acgt = sum(seq.count(x) for x in "ACGT")
    return 0.0 if acgt == 0 else 100.0 * (seq.count("G") + seq.count("C")) / acgt


def select_representatives(meta: Dict[str, dict], levels: List[str], policy: str) -> Dict[str, List[str]]:
    """
    Return dict component -> selected sequence IDs.
    component names are family_tiling and haplogroup_tiling.
    """
    selected: Dict[str, List[str]] = defaultdict(list)

    if "family" in levels:
        # Prefer explicitly flagged rows, otherwise first sequence per family.
        flagged = [sid for sid, row in meta.items() if is_true(row.get("used_for_family_tiling", ""))]
        if flagged and policy in {"flagged", "flagged_then_first"}:
            selected["family_tiling"].extend(flagged)
        elif policy != "flagged":
            seen = set()
            for sid, row in meta.items():
                fam = row.get("family", "")
                if fam and fam not in seen:
                    selected["family_tiling"].append(sid)
                    seen.add(fam)

    if "haplogroup" in levels or "lineage" in levels:
        flagged = [sid for sid, row in meta.items() if is_true(row.get("used_for_haplogroup_tiling", ""))]
        if flagged and policy in {"flagged", "flagged_then_first"}:
            selected["haplogroup_tiling"].extend(flagged)
        elif policy != "flagged":
            seen = set()
            for sid, row in meta.items():
                hap = row.get("haplogroup", "") or row.get("lineage", "") or row.get("population", "")
                # This naturally restricts to Homo/domestic lineages if only those rows have haplogroup labels.
                if hap and hap not in seen:
                    selected["haplogroup_tiling"].append(sid)
                    seen.add(hap)

    return selected


def main() -> None:
    p = argparse.ArgumentParser(description="Generate supplementary family/haplogroup tiling probes.")
    p.add_argument("--fasta", required=True, help="Clean reference panel FASTA.")
    p.add_argument("--metadata", required=True, help="Clean reference metadata TSV.")
    p.add_argument("--levels", default="family,haplogroup", help="Comma-separated levels: family,haplogroup.")
    p.add_argument("--representative-policy", default="flagged_then_first",
                   choices=["flagged_then_first", "first", "flagged"],
                   help="How to select representatives for tiling.")
    p.add_argument("--probe-len", type=int, default=52, help="Probe length.")
    p.add_argument("--step", type=int, default=5, help="Tiling step size.")
    p.add_argument("--max-n", type=int, default=0, help="Maximum number of Ns allowed in a tile.")
    p.add_argument("--out-fasta", required=True, help="Output supplementary probe FASTA.")
    p.add_argument("--out-metadata", required=True, help="Output supplementary probe metadata TSV.")
    args = p.parse_args()

    levels = [x.strip().lower() for x in args.levels.split(",") if x.strip()]
    meta = read_metadata(args.metadata)
    seqs = {seq_id_from_header(h): s for h, _d, s in parse_fasta(args.fasta)}
    selected = select_representatives(meta, levels, args.representative_policy)

    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_metadata).parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "probe_id", "sequence", "length", "gc_content",
        "reference_id", "start", "end", "strand",
        "family", "genus", "species", "group", "haplogroup", "population",
        "component", "tiling_level"
    ]

    n_tiles = 0
    selected_count = 0

    with open_text(args.out_fasta, "wt") as out_fa, open_text(args.out_metadata, "wt") as out_meta:
        writer = csv.DictWriter(out_meta, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        for component, ids in selected.items():
            for ref_id in ids:
                if ref_id not in seqs:
                    print(f"[WARN] selected reference missing from FASTA: {ref_id}", file=sys.stderr)
                    continue
                selected_count += 1
                row = meta.get(ref_id, {})
                seq = seqs[ref_id]

                for start0 in range(0, len(seq) - args.probe_len + 1, args.step):
                    probe_seq = seq[start0:start0 + args.probe_len].upper()
                    if probe_seq.count("N") > args.max_n or not re.fullmatch(r"[ACGTN]+", probe_seq):
                        continue

                    n_tiles += 1
                    level = component.replace("_tiling", "")
                    probe_id = (
                        f"{component}_{n_tiles:09d}|ref={ref_id}|"
                        f"start={start0 + 1}|end={start0 + args.probe_len}|strand=+"
                    )
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
                        "component": component,
                        "tiling_level": level,
                    })

    print(f"[INFO] selected references: {selected_count}", file=sys.stderr)
    for comp, ids in selected.items():
        print(f"[INFO] {comp}: {len(ids)} references", file=sys.stderr)
    print(f"[INFO] supplementary tiles: {n_tiles}", file=sys.stderr)
    print(f"[INFO] output FASTA       : {args.out_fasta}", file=sys.stderr)
    print(f"[INFO] output metadata    : {args.out_metadata}", file=sys.stderr)


if __name__ == "__main__":
    main()
