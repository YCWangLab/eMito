#!/usr/bin/env python3
"""
08_merge_deduplicate_probes.py

Merge filtered informative probes and supplementary tiling probes, then remove
duplicate probe sequences.

This implements the final eMito design step:
  all filtered and supplementary probes are merged and deduplicated.

Example
-------
python scripts/design/08_merge_deduplicate_probes.py \
  --fasta results/design/06_overlap_reduced/informative_overlap_reduced.fasta \
          results/design/07_supplementary_tiling/family_haplogroup_tiling.fasta \
  --metadata results/design/06_overlap_reduced/informative_overlap_reduced.metadata.tsv \
             results/design/07_supplementary_tiling/family_haplogroup_tiling.metadata.tsv \
  --dedup-mode exact \
  --out-fasta probes/emito_probes.fasta.gz \
  --out-metadata probes/emito_probes.metadata.tsv.gz \
  --summary results/design/08_final/final_probe_summary.tsv
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


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def canonical_key(seq: str, mode: str) -> str:
    seq = seq.upper()
    if mode == "exact":
        return seq
    rc = revcomp(seq)
    return min(seq, rc)


def gc_percent(seq: str) -> float:
    acgt = sum(seq.count(x) for x in "ACGT")
    return 0.0 if acgt == 0 else 100.0 * (seq.count("G") + seq.count("C")) / acgt


def read_metadata(paths: List[str]) -> Dict[str, dict]:
    out: Dict[str, dict] = {}
    for path in paths or []:
        with open_text(path, "rt") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if not reader.fieldnames or "probe_id" not in reader.fieldnames:
                sys.exit(f"[ERROR] metadata file lacks probe_id column: {path}")
            for row in reader:
                row["_metadata_file"] = str(path)
                out[row["probe_id"]] = row
    return out


def uniq_join(values: List[str]) -> str:
    vals = []
    for v in values:
        if v is None:
            continue
        for part in str(v).split(";"):
            part = part.strip()
            if part and part not in vals:
                vals.append(part)
    return ";".join(vals)


def main() -> None:
    p = argparse.ArgumentParser(description="Merge and deduplicate eMito probe FASTA files.")
    p.add_argument("--fasta", nargs="+", required=True, help="Input probe FASTA files.")
    p.add_argument("--metadata", nargs="*", default=[], help="Input probe metadata TSV files.")
    p.add_argument("--dedup-mode", default="exact", choices=["exact", "canonical_rc"],
                   help="Deduplicate exact sequences only, or collapse reverse-complement pairs too.")
    p.add_argument("--id-prefix", default="eMito_probe", help="Prefix for final probe IDs.")
    p.add_argument("--out-fasta", required=True, help="Final deduplicated probe FASTA or FASTA.GZ.")
    p.add_argument("--out-metadata", required=True, help="Final deduplicated metadata TSV or TSV.GZ.")
    p.add_argument("--summary", default=None, help="Optional summary TSV.")
    args = p.parse_args()

    meta = read_metadata(args.metadata)

    merged: Dict[str, dict] = {}
    input_counts = defaultdict(int)

    for fasta_path in args.fasta:
        source_file = Path(fasta_path).name
        for probe_id, desc, seq in parse_fasta(fasta_path):
            seq = seq.upper()
            if not re.fullmatch(r"[ACGTN]+", seq):
                print(f"[WARN] skipped invalid probe sequence: {probe_id}", file=sys.stderr)
                continue
            input_counts[source_file] += 1
            key = canonical_key(seq, args.dedup_mode)
            row = meta.get(probe_id, {})
            if key not in merged:
                merged[key] = {
                    "sequence": seq,
                    "source_probe_ids": [],
                    "source_files": [],
                    "components": [],
                    "families": [],
                    "genera": [],
                    "species": [],
                    "groups": [],
                    "haplogroups": [],
                    "populations": [],
                    "reference_ids": [],
                    "selection_reasons": [],
                }
            m = merged[key]
            m["source_probe_ids"].append(probe_id)
            m["source_files"].append(source_file)
            m["components"].append(row.get("component", ""))
            m["families"].append(row.get("family", ""))
            m["genera"].append(row.get("genus", ""))
            m["species"].append(row.get("species", ""))
            m["groups"].append(row.get("group", ""))
            m["haplogroups"].append(row.get("haplogroup", ""))
            m["populations"].append(row.get("population", ""))
            m["reference_ids"].append(row.get("reference_id", ""))
            m["selection_reasons"].append(row.get("selection_reason", ""))

    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_metadata).parent.mkdir(parents=True, exist_ok=True)
    if args.summary:
        Path(args.summary).parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "probe_id", "sequence", "length", "gc_content",
        "n_merged_records", "source_probe_ids", "source_files",
        "components", "families", "genera", "species", "groups",
        "haplogroups", "populations", "reference_ids", "selection_reasons",
        "dedup_mode"
    ]

    component_counter = defaultdict(int)

    with open_text(args.out_fasta, "wt") as out_fa, open_text(args.out_metadata, "wt") as out_meta:
        writer = csv.DictWriter(out_meta, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        for idx, key in enumerate(sorted(merged), start=1):
            m = merged[key]
            final_id = f"{args.id_prefix}_{idx:09d}"
            seq = m["sequence"]
            write_fasta_record(out_fa, final_id, seq)
            components = uniq_join(m["components"])
            for comp in components.split(";"):
                if comp:
                    component_counter[comp] += 1

            writer.writerow({
                "probe_id": final_id,
                "sequence": seq,
                "length": len(seq),
                "gc_content": f"{gc_percent(seq):.4f}",
                "n_merged_records": len(m["source_probe_ids"]),
                "source_probe_ids": uniq_join(m["source_probe_ids"]),
                "source_files": uniq_join(m["source_files"]),
                "components": components,
                "families": uniq_join(m["families"]),
                "genera": uniq_join(m["genera"]),
                "species": uniq_join(m["species"]),
                "groups": uniq_join(m["groups"]),
                "haplogroups": uniq_join(m["haplogroups"]),
                "populations": uniq_join(m["populations"]),
                "reference_ids": uniq_join(m["reference_ids"]),
                "selection_reasons": uniq_join(m["selection_reasons"]),
                "dedup_mode": args.dedup_mode,
            })

    if args.summary:
        with open_text(args.summary, "wt") as out:
            out.write("section\tname\tvalue\n")
            for name, count in sorted(input_counts.items()):
                out.write(f"input_records\t{name}\t{count}\n")
            out.write(f"final_unique_probes\tall\t{len(merged)}\n")
            for comp, count in sorted(component_counter.items()):
                out.write(f"final_component_count\t{comp}\t{count}\n")

    print(f"[INFO] input files          : {len(args.fasta)}", file=sys.stderr)
    print(f"[INFO] final unique probes  : {len(merged)}", file=sys.stderr)
    print(f"[INFO] output FASTA         : {args.out_fasta}", file=sys.stderr)
    print(f"[INFO] output metadata      : {args.out_metadata}", file=sys.stderr)


if __name__ == "__main__":
    main()
