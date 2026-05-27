#!/usr/bin/env python3
"""
06_reduce_overlapping_probes.py

Reduce redundancy among filtered informative probes by grouping overlapping
probe coordinates and retaining the central probe from each overlap cluster.

This implements the eMito logic:
  filtered 52 bp probes are grouped by overlapping genomic coordinates on each
  reference genome, and only the central probe from each overlapping cluster is
  retained.

Example
-------
python scripts/design/06_reduce_overlapping_probes.py \
  --probes-fasta results/design/05_filtered_probes/filtered_informative_probes.fasta \
  --metadata results/design/05_filtered_probes/filtered_informative_probes.metadata.tsv \
  --out-fasta results/design/06_overlap_reduced/informative_overlap_reduced.fasta \
  --out-metadata results/design/06_overlap_reduced/informative_overlap_reduced.metadata.tsv \
  --cluster-summary results/design/06_overlap_reduced/overlap_clusters.tsv
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


def read_metadata(path: str) -> Tuple[Dict[str, dict], List[str]]:
    with open_text(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or "probe_id" not in reader.fieldnames:
            sys.exit("[ERROR] metadata must contain probe_id")
        rows = {row["probe_id"]: row for row in reader}
        return rows, list(reader.fieldnames)


def parse_coord_from_probe_id(probe_id: str) -> dict:
    out = {}
    for key in ["ref", "start", "end", "strand"]:
        m = re.search(rf"{key}=([^|]+)", probe_id)
        if m:
            out[key] = m.group(1)
    if "ref" in out:
        out["reference_id"] = out["ref"]
    return out


def to_int(value, field: str, probe_id: str) -> int:
    try:
        return int(float(value))
    except Exception:
        sys.exit(f"[ERROR] cannot parse integer {field}={value!r} for {probe_id}")


def main() -> None:
    p = argparse.ArgumentParser(description="Reduce overlapping probes by keeping central probe per overlap cluster.")
    p.add_argument("--probes-fasta", required=True, help="Filtered probe FASTA.")
    p.add_argument("--metadata", required=True, help="Filtered probe metadata TSV.")
    p.add_argument("--group-columns", default="reference_id,strand",
                   help="Comma-separated metadata columns used before overlap clustering.")
    p.add_argument("--out-fasta", required=True, help="Overlap-reduced probe FASTA.")
    p.add_argument("--out-metadata", required=True, help="Overlap-reduced metadata TSV.")
    p.add_argument("--cluster-summary", default=None, help="Optional cluster summary TSV.")
    args = p.parse_args()

    meta, meta_fields = read_metadata(args.metadata)
    seqs = {probe_id: seq for probe_id, _desc, seq in parse_fasta(args.probes_fasta)}

    records = []
    for probe_id, row in meta.items():
        if probe_id not in seqs:
            continue
        merged = dict(row)
        parsed = parse_coord_from_probe_id(probe_id)
        for k, v in parsed.items():
            merged.setdefault(k, v)
        if "reference_id" not in merged or not merged.get("reference_id"):
            sys.exit(f"[ERROR] missing reference_id for {probe_id}")
        if "start" not in merged or "end" not in merged:
            sys.exit(f"[ERROR] missing start/end coordinates for {probe_id}")

        start = to_int(merged["start"], "start", probe_id)
        end = to_int(merged["end"], "end", probe_id)
        if start > end:
            start, end = end, start
        merged["start"] = start
        merged["end"] = end
        merged.setdefault("strand", "+")
        merged["_seq"] = seqs[probe_id]
        merged["_mid"] = (start + end) / 2.0
        records.append(merged)

    group_columns = [x.strip() for x in args.group_columns.split(",") if x.strip()]
    grouped = defaultdict(list)
    for row in records:
        key = tuple(row.get(c, "") for c in group_columns)
        grouped[key].append(row)

    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_metadata).parent.mkdir(parents=True, exist_ok=True)
    if args.cluster_summary:
        Path(args.cluster_summary).parent.mkdir(parents=True, exist_ok=True)

    out_fields = list(dict.fromkeys(meta_fields + [
        "cluster_id", "cluster_size", "cluster_start", "cluster_end", "selected_by"
    ]))

    kept_rows = []
    cluster_rows = []
    cluster_index = 0

    for key, rows in sorted(grouped.items()):
        rows = sorted(rows, key=lambda r: (int(r["start"]), int(r["end"]), r["probe_id"]))
        cluster: List[dict] = []
        cluster_end = None

        def flush_cluster(cl: List[dict]):
            nonlocal cluster_index
            if not cl:
                return
            cluster_index += 1
            c_start = min(int(r["start"]) for r in cl)
            c_end = max(int(r["end"]) for r in cl)
            c_mid = (c_start + c_end) / 2.0
            selected = min(cl, key=lambda r: (abs(float(r["_mid"]) - c_mid), r["probe_id"]))
            selected = dict(selected)
            selected["cluster_id"] = f"overlap_cluster_{cluster_index:09d}"
            selected["cluster_size"] = len(cl)
            selected["cluster_start"] = c_start
            selected["cluster_end"] = c_end
            selected["selected_by"] = "central_probe_in_overlap_cluster"
            kept_rows.append(selected)
            cluster_rows.append({
                "cluster_id": selected["cluster_id"],
                "group_key": "|".join(map(str, key)),
                "cluster_size": len(cl),
                "cluster_start": c_start,
                "cluster_end": c_end,
                "selected_probe_id": selected["probe_id"],
                "all_probe_ids": ";".join(r["probe_id"] for r in cl),
            })

        for row in rows:
            if not cluster:
                cluster = [row]
                cluster_end = int(row["end"])
            elif int(row["start"]) <= int(cluster_end):
                cluster.append(row)
                cluster_end = max(int(cluster_end), int(row["end"]))
            else:
                flush_cluster(cluster)
                cluster = [row]
                cluster_end = int(row["end"])

        flush_cluster(cluster)

    with open_text(args.out_fasta, "wt") as out_fa, open_text(args.out_metadata, "wt") as out_meta:
        writer = csv.DictWriter(out_meta, fieldnames=out_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in kept_rows:
            write_fasta_record(out_fa, row["probe_id"], row["_seq"])
            writer.writerow({k: row.get(k, "") for k in out_fields})

    if args.cluster_summary:
        with open_text(args.cluster_summary, "wt") as out:
            fieldnames = ["cluster_id", "group_key", "cluster_size", "cluster_start", "cluster_end", "selected_probe_id", "all_probe_ids"]
            writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in cluster_rows:
                writer.writerow(row)

    print(f"[INFO] input probes       : {len(records)}", file=sys.stderr)
    print(f"[INFO] overlap clusters   : {len(cluster_rows)}", file=sys.stderr)
    print(f"[INFO] retained probes    : {len(kept_rows)}", file=sys.stderr)
    print(f"[INFO] output FASTA       : {args.out_fasta}", file=sys.stderr)
    print(f"[INFO] output metadata    : {args.out_metadata}", file=sys.stderr)


if __name__ == "__main__":
    main()
