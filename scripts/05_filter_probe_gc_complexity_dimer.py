from pathlib import Path
import textwrap, subprocess, sys, os, zipfile

script = r'''#!/usr/bin/env python3
"""
05_filter_probe_gc_complexity_dimer.py

Filter eMito candidate probes by GC content, DUST-style 3-mer complexity,
and probe-pool dimer score.

This script reproduces the logic of the original eMito probe-filtering
workflow used for files such as:

    *.filtered.probe.fasta

Original filtering logic
------------------------
For each input FASTA:

1. Calculate GC content:
       GC = (G + C) / (A + C + G + T) * 100

   This is equivalent to Bio.SeqUtils.gc_fraction(seq) * 100 for A/C/G/T
   probe sequences.

2. Calculate DUST-style complexity score:
       k = 3
       all 3-mers are counted
       Complexity = sum(Ct * (Ct - 1) / 2) / (len(seq) - 3)

   where Ct is the count of each 3-mer in the probe.

3. Calculate dimer score:
       k = 11 by default
       reverse-complement all probes
       count all 11-mers in the reverse-complement probe pool
       keep only reverse-complement 11-mers with frequency >= 2 by default

   For each original probe:
       Dimer_Score =
           sum(freq of each probe 11-mer in filtered RC 11-mer table)
           / number_of_probes * 100

4. Filter:
       GC between 40 and 60, inclusive
       Complexity between 0 and 2, inclusive

5. Sort GC/complexity-passing probes by Dimer_Score ascending.

6. Retain the lowest 85% by Dimer_Score.

Example
-------
python scripts/design/05_filter_probe_gc_complexity_dimer.py \
  --probes-fasta results/design/03_probe_tiles/representative_tiles.fasta \
  --probes-metadata results/design/03_probe_tiles/representative_tiles.metadata.tsv \
  --gc-min 40 \
  --gc-max 60 \
  --complexity-min 0 \
  --complexity-max 2 \
  --dimer-k 11 \
  --dimer-min-freq 2 \
  --keep-lowest-dimer-fraction 0.85 \
  --out-fasta results/design/05_filtered_probes/filtered_probes.fasta \
  --out-metadata results/design/05_filtered_probes/filtered_probes.metadata.tsv \
  --all-scores results/design/05_filtered_probes/all_probe_gc_complexity_dimer.tsv \
  --rejected results/design/05_filtered_probes/rejected_probes.tsv

Optional informative-kmer mode
------------------------------
If --informative-kmers is provided, the script can additionally require probes
to contain at least one species/genus-informative k-mer. This is optional
because the original input files (*.filtered.probe.fasta) were already
pre-screened upstream.

Use:

  --informative-kmers species_genus_specific_kmers.tsv.gz
  --require-informative-kmer

"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import sys
from collections import Counter, OrderedDict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple


# ---------------------------------------------------------------------
# Basic I/O
# ---------------------------------------------------------------------

def open_text(path: str | Path, mode: str = "rt"):
    """Open plain-text or gzipped text files."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, newline="" if "t" in mode else None)


def parse_fasta(path: str | Path) -> Iterator[Tuple[str, str, str]]:
    """
    Parse FASTA.

    Yields
    ------
    record_id : str
        First token after ">".
    description : str
        Full FASTA header without ">".
    sequence : str
        Uppercase sequence with U converted to T.
    """
    record_id: Optional[str] = None
    description = ""
    chunks: List[str] = []

    with open_text(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith(">"):
                if record_id is not None:
                    yield record_id, description, "".join(chunks).upper().replace("U", "T")

                description = line[1:].strip()
                record_id = description.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())

        if record_id is not None:
            yield record_id, description, "".join(chunks).upper().replace("U", "T")


def write_fasta_record(handle, record_id: str, seq: str, width: int = 80) -> None:
    """Write one FASTA record."""
    handle.write(f">{record_id}\n")
    for i in range(0, len(seq), width):
        handle.write(seq[i:i + width] + "\n")


def read_metadata(path: Optional[str | Path]) -> Tuple[Dict[str, dict], List[str]]:
    """
    Read probe metadata TSV.

    The metadata table must contain a 'probe_id' column if provided.
    """
    if not path:
        return {}, []

    with open_text(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            sys.exit(f"[ERROR] Empty metadata file: {path}")
        if "probe_id" not in reader.fieldnames:
            sys.exit(f"[ERROR] Metadata file must contain a 'probe_id' column: {path}")

        rows = OrderedDict()
        for row in reader:
            pid = row.get("probe_id", "")
            if pid:
                rows[pid] = row

        return rows, list(reader.fieldnames)


def read_informative_kmers(path: Optional[str | Path]) -> Set[str]:
    """
    Read informative k-mers from a TSV/TSV.GZ file.

    Expected column:
      kmer
    """
    kmers: Set[str] = set()
    if not path:
        return kmers

    with open_text(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or "kmer" not in reader.fieldnames:
            sys.exit(f"[ERROR] Informative k-mer table must contain a 'kmer' column: {path}")
        for row in reader:
            k = str(row.get("kmer", "")).strip().upper()
            if k:
                kmers.add(k)

    return kmers


# ---------------------------------------------------------------------
# Sequence utility functions
# ---------------------------------------------------------------------

def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


def is_acgtn(seq: str) -> bool:
    """Return True if sequence contains only A/C/G/T/N."""
    return bool(re.fullmatch(r"[ACGTN]+", seq.upper()))


def gc_content_percent(seq: str) -> float:
    """
    Calculate GC percentage.

    For standard A/C/G/T probe sequences, this is equivalent to:
      Bio.SeqUtils.gc_fraction(seq) * 100

    Ambiguous characters are ignored from the denominator.
    """
    seq = seq.upper()
    a = seq.count("A")
    c = seq.count("C")
    g = seq.count("G")
    t = seq.count("T")
    denom = a + c + g + t
    if denom == 0:
        return 0.0
    return round(((g + c) / denom) * 100.0, 4)


def complexity_score_original(seq: str, k: int = 3) -> float:
    """
    Calculate the original eMito DUST-style 3-mer complexity score.

    Original formula:
      kmers_list = [seq[i:i+k] for i in range(len(seq)-k+1)]
      Ct_dict = count of each kmer
      all_ct = sum(ct * (ct - 1) / 2 for ct in Ct_dict.values())
      score = all_ct / (len(seq) - 3)

    For the original workflow, k = 3. The denominator is intentionally
    len(seq) - 3, not the number of 3-mers len(seq) - 2.
    """
    seq = seq.upper()
    if len(seq) <= k:
        return 999.0

    kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    counts = Counter(kmers)
    all_ct = sum(ct * (ct - 1) / 2 for ct in counts.values())

    # Keep original denominator used in the user's previous script.
    denom = len(seq) - 3
    if denom <= 0:
        return 999.0

    return round(all_ct / denom, 4)


def iter_kmers(seq: str, k: int) -> Iterator[str]:
    """Yield all contiguous k-mers from seq."""
    seq = seq.upper()
    if len(seq) < k:
        return
    for i in range(len(seq) - k + 1):
        yield seq[i:i + k]


def count_rc_kmers_for_probe_pool(records: Sequence[Tuple[str, str]], k: int) -> Counter:
    """
    Count k-mers in the reverse-complement sequences of all probes.

    Parameters
    ----------
    records
        Sequence of (probe_id, sequence).
    k
        k-mer size, default is 11 in the original workflow.
    """
    rc_kmer_counts: Counter = Counter()

    for _probe_id, seq in records:
        rc = reverse_complement(seq)
        for kmer in iter_kmers(rc, k):
            rc_kmer_counts[kmer] += 1

    return rc_kmer_counts


def dimer_score_original(
    seq: str,
    filtered_rc_kmer_counts: Dict[str, int],
    probe_num: int,
    k: int = 11
) -> float:
    """
    Calculate original eMito dimer score for one probe.

    Original formula:
      seq_kmers = all original probe k-mers
      dimer_score = sum(filtered_rc_kmer_counts.get(kmer, 0) for kmer in seq_kmers)
      Dimer_Score = (dimer_score / probe_num) * 100
    """
    if probe_num <= 0:
        return 0.0

    score = 0
    for kmer in iter_kmers(seq, k):
        score += filtered_rc_kmer_counts.get(kmer, 0)

    return (score / probe_num) * 100.0


def contains_informative_kmer(seq: str, informative_kmers: Set[str], include_rc: bool = True) -> bool:
    """
    Return True if seq contains one informative k-mer.

    This is optional and not part of the original GC/complexity/dimer script,
    because the original input FASTA files were already upstream-filtered.
    """
    if not informative_kmers:
        return True

    k = len(next(iter(informative_kmers)))
    seqs = [seq.upper()]
    if include_rc:
        seqs.append(reverse_complement(seq))

    for s in seqs:
        for kmer in iter_kmers(s, k):
            if kmer in informative_kmers:
                return True

    return False


# ---------------------------------------------------------------------
# Main filtering logic
# ---------------------------------------------------------------------

def build_score_rows(
    records: Sequence[Tuple[str, str]],
    metadata: Dict[str, dict],
    informative_kmers: Set[str],
    require_informative_kmer: bool,
    dimer_k: int,
    dimer_min_freq: int,
    allow_n: bool
) -> List[dict]:
    """
    Calculate GC, complexity, and dimer score for all probes.

    Dimer score must be computed using the full probe pool, so all records
    are loaded before scoring.
    """
    probe_num = len(records)

    rc_kmer_counts = count_rc_kmers_for_probe_pool(records, dimer_k)
    filtered_rc_kmer_counts = {
        kmer: count
        for kmer, count in rc_kmer_counts.items()
        if count >= dimer_min_freq
    }

    rows: List[dict] = []

    for probe_id, seq in records:
        seq = seq.upper()
        valid_sequence = is_acgtn(seq) and (allow_n or "N" not in seq)
        gc = gc_content_percent(seq)
        complexity = complexity_score_original(seq, k=3)
        dimer = dimer_score_original(
            seq=seq,
            filtered_rc_kmer_counts=filtered_rc_kmer_counts,
            probe_num=probe_num,
            k=dimer_k,
        )

        has_info = contains_informative_kmer(seq, informative_kmers, include_rc=True)
        pass_info = has_info if require_informative_kmer else True

        row = dict(metadata.get(probe_id, {}))
        row["probe_id"] = probe_id
        row["sequence"] = seq
        row["length"] = len(seq)
        row["valid_sequence"] = int(valid_sequence)
        row["gc_content"] = f"{gc:.4f}"
        row["Complexity"] = f"{complexity:.4f}"
        row["Dimer_Score"] = f"{dimer:.8f}"
        row["contains_informative_kmer"] = int(has_info)
        row["pass_informative_kmer"] = int(pass_info)
        row["_gc_float"] = gc
        row["_complexity_float"] = complexity
        row["_dimer_float"] = dimer
        row["_valid_bool"] = valid_sequence
        row["_pass_info_bool"] = pass_info

        rows.append(row)

    return rows


def assign_filter_status(
    rows: List[dict],
    gc_min: float,
    gc_max: float,
    complexity_min: float,
    complexity_max: float,
    keep_lowest_dimer_fraction: float
) -> Tuple[List[dict], List[dict]]:
    """
    Apply original filtering logic:

      1. pass valid sequence;
      2. pass optional informative-kmer check;
      3. GC between gc_min and gc_max inclusive;
      4. Complexity between complexity_min and complexity_max inclusive;
      5. among these, keep the lowest dimer-score fraction.

    Returns
    -------
    kept_rows, rejected_rows
    """
    gc_complexity_pass: List[dict] = []
    rejected_pre_dimer: List[dict] = []

    for row in rows:
        reasons: List[str] = []

        if not row["_valid_bool"]:
            reasons.append("invalid_or_ambiguous_sequence")
        if not row["_pass_info_bool"]:
            reasons.append("no_informative_kmer")
        if not (gc_min <= row["_gc_float"] <= gc_max):
            reasons.append("gc_out_of_range")
        if not (complexity_min <= row["_complexity_float"] <= complexity_max):
            reasons.append("complexity_out_of_range")

        if reasons:
            row["filter_status"] = "rejected"
            row["reject_reason"] = ";".join(reasons)
            row["selection_reason"] = ""
            rejected_pre_dimer.append(row)
        else:
            gc_complexity_pass.append(row)

    gc_complexity_pass = sorted(
        gc_complexity_pass,
        key=lambda r: (r["_dimer_float"], r["probe_id"])
    )

    if not (0 < keep_lowest_dimer_fraction <= 1):
        sys.exit("[ERROR] --keep-lowest-dimer-fraction must be > 0 and <= 1")

    # Original behavior used int(0.85 * len(df_filtered)).
    n_keep = int(keep_lowest_dimer_fraction * len(gc_complexity_pass))

    kept_rows: List[dict] = []
    rejected_dimer: List[dict] = []

    for idx, row in enumerate(gc_complexity_pass):
        row["dimer_rank"] = idx + 1
        row["n_gc_complexity_pass"] = len(gc_complexity_pass)
        row["keep_lowest_dimer_fraction"] = keep_lowest_dimer_fraction

        if idx < n_keep:
            row["filter_status"] = "kept"
            row["reject_reason"] = ""
            row["selection_reason"] = (
                "gc_complexity_pass_lowest_"
                f"{keep_lowest_dimer_fraction:g}_dimer_fraction"
            )
            kept_rows.append(row)
        else:
            row["filter_status"] = "rejected"
            row["reject_reason"] = "high_dimer_score_top_fraction_removed"
            row["selection_reason"] = ""
            rejected_dimer.append(row)

    rejected_rows = rejected_pre_dimer + rejected_dimer
    return kept_rows, rejected_rows


def clean_output_row(row: dict, fieldnames: List[str]) -> dict:
    """Remove private helper fields and return row restricted to fieldnames."""
    return {k: row.get(k, "") for k in fieldnames}


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Filter candidate probes by original eMito GC, 3-mer complexity, "
            "and probe-pool dimer score logic."
        )
    )

    parser.add_argument(
        "--probes-fasta",
        required=True,
        help="Input candidate probe FASTA."
    )
    parser.add_argument(
        "--probes-metadata",
        default=None,
        help="Optional probe metadata TSV with a probe_id column."
    )

    parser.add_argument(
        "--gc-min",
        type=float,
        default=40.0,
        help="Minimum GC percentage, inclusive. Default: 40."
    )
    parser.add_argument(
        "--gc-max",
        type=float,
        default=60.0,
        help="Maximum GC percentage, inclusive. Default: 60."
    )
    parser.add_argument(
        "--complexity-min",
        type=float,
        default=0.0,
        help="Minimum complexity score, inclusive. Default: 0."
    )
    parser.add_argument(
        "--complexity-max",
        type=float,
        default=2.0,
        help="Maximum complexity score, inclusive. Default: 2."
    )
    parser.add_argument(
        "--dimer-k",
        type=int,
        default=11,
        help="k-mer size for dimer scoring. Default: 11."
    )
    parser.add_argument(
        "--dimer-min-freq",
        type=int,
        default=2,
        help="Minimum RC k-mer frequency retained in dimer scoring. Default: 2."
    )
    parser.add_argument(
        "--keep-lowest-dimer-fraction",
        type=float,
        default=0.85,
        help=(
            "Fraction of GC/complexity-passing probes retained after sorting "
            "by Dimer_Score ascending. Default: 0.85."
        )
    )
    parser.add_argument(
        "--allow-n",
        action="store_true",
        help="Allow probes containing N. By default, probes with N are rejected."
    )

    parser.add_argument(
        "--informative-kmers",
        default=None,
        help=(
            "Optional species/genus-specific k-mer table with a 'kmer' column. "
            "Not required for reproducing the original GC/complexity/dimer filter."
        )
    )
    parser.add_argument(
        "--require-informative-kmer",
        action="store_true",
        help=(
            "If set, require each retained probe to contain an informative k-mer. "
            "Only meaningful together with --informative-kmers."
        )
    )

    parser.add_argument(
        "--out-fasta",
        required=True,
        help="Output filtered probe FASTA."
    )
    parser.add_argument(
        "--out-metadata",
        required=True,
        help="Output filtered probe metadata TSV."
    )
    parser.add_argument(
        "--all-scores",
        default=None,
        help="Optional TSV containing scores and status for all probes."
    )
    parser.add_argument(
        "--rejected",
        default=None,
        help="Optional TSV containing rejected probes and rejection reasons."
    )

    args = parser.parse_args()

    if args.dimer_k <= 0:
        sys.exit("[ERROR] --dimer-k must be positive")
    if args.dimer_min_freq <= 0:
        sys.exit("[ERROR] --dimer-min-freq must be positive")

    metadata, metadata_fields = read_metadata(args.probes_metadata)

    records: List[Tuple[str, str]] = []
    duplicate_ids: Set[str] = set()
    seen_ids: Set[str] = set()

    for probe_id, _description, seq in parse_fasta(args.probes_fasta):
        if probe_id in seen_ids:
            duplicate_ids.add(probe_id)
        seen_ids.add(probe_id)
        records.append((probe_id, seq))

    if not records:
        sys.exit(f"[ERROR] No FASTA records found: {args.probes_fasta}")

    if duplicate_ids:
        print(
            f"[WARN] duplicated FASTA IDs detected: {len(duplicate_ids)}. "
            "Metadata matching may be ambiguous.",
            file=sys.stderr,
        )

    informative_kmers = read_informative_kmers(args.informative_kmers)
    if args.require_informative_kmer and not informative_kmers:
        sys.exit("[ERROR] --require-informative-kmer was set but --informative-kmers is empty or missing")

    print(f"[INFO] input probes                 : {len(records)}", file=sys.stderr)
    print(f"[INFO] metadata rows                : {len(metadata)}", file=sys.stderr)
    print(f"[INFO] informative kmers            : {len(informative_kmers)}", file=sys.stderr)
    print(f"[INFO] GC filter                    : {args.gc_min} <= GC <= {args.gc_max}", file=sys.stderr)
    print(
        f"[INFO] Complexity filter            : "
        f"{args.complexity_min} <= Complexity <= {args.complexity_max}",
        file=sys.stderr,
    )
    print(f"[INFO] Dimer k                      : {args.dimer_k}", file=sys.stderr)
    print(f"[INFO] Dimer RC kmer min frequency  : {args.dimer_min_freq}", file=sys.stderr)
    print(
        f"[INFO] Keep lowest dimer fraction   : {args.keep_lowest_dimer_fraction}",
        file=sys.stderr,
    )

    rows = build_score_rows(
        records=records,
        metadata=metadata,
        informative_kmers=informative_kmers,
        require_informative_kmer=args.require_informative_kmer,
        dimer_k=args.dimer_k,
        dimer_min_freq=args.dimer_min_freq,
        allow_n=args.allow_n,
    )

    kept_rows, rejected_rows = assign_filter_status(
        rows=rows,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        complexity_min=args.complexity_min,
        complexity_max=args.complexity_max,
        keep_lowest_dimer_fraction=args.keep_lowest_dimer_fraction,
    )

    # Output fields.
    score_fields = [
        "probe_id",
        "sequence",
        "length",
        "valid_sequence",
        "gc_content",
        "Complexity",
        "Dimer_Score",
        "contains_informative_kmer",
        "pass_informative_kmer",
        "dimer_rank",
        "n_gc_complexity_pass",
        "keep_lowest_dimer_fraction",
        "filter_status",
        "reject_reason",
        "selection_reason",
    ]

    # Preserve original metadata fields first, then append score/filter fields.
    out_fields = list(dict.fromkeys(metadata_fields + score_fields))
    if "probe_id" not in out_fields:
        out_fields.insert(0, "probe_id")

    private_fields = {
        "_gc_float",
        "_complexity_float",
        "_dimer_float",
        "_valid_bool",
        "_pass_info_bool",
    }

    final_fields = [f for f in out_fields if f not in private_fields]

    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_metadata).parent.mkdir(parents=True, exist_ok=True)
    if args.all_scores:
        Path(args.all_scores).parent.mkdir(parents=True, exist_ok=True)
    if args.rejected:
        Path(args.rejected).parent.mkdir(parents=True, exist_ok=True)

    # Write kept FASTA and metadata.
    with open_text(args.out_fasta, "wt") as out_fa, open_text(args.out_metadata, "wt") as out_meta:
        writer = csv.DictWriter(out_meta, fieldnames=final_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        for row in kept_rows:
            write_fasta_record(out_fa, row["probe_id"], row["sequence"])
            writer.writerow(clean_output_row(row, final_fields))

    # Write all scores.
    if args.all_scores:
        all_rows = sorted(
            kept_rows + rejected_rows,
            key=lambda r: (r.get("filter_status", ""), r.get("probe_id", ""))
        )
        with open_text(args.all_scores, "wt") as out:
            writer = csv.DictWriter(out, fieldnames=final_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in all_rows:
                writer.writerow(clean_output_row(row, final_fields))

    # Write rejected table.
    if args.rejected:
        rejected_fields = final_fields
        with open_text(args.rejected, "wt") as out:
            writer = csv.DictWriter(out, fieldnames=rejected_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rejected_rows:
                writer.writerow(clean_output_row(row, rejected_fields))

    print(f"[INFO] GC/complexity/dimer kept    : {len(kept_rows)}", file=sys.stderr)
    print(f"[INFO] rejected probes             : {len(rejected_rows)}", file=sys.stderr)
    print(f"[INFO] output FASTA                : {args.out_fasta}", file=sys.stderr)
    print(f"[INFO] output metadata             : {args.out_metadata}", file=sys.stderr)
    if args.all_scores:
        print(f"[INFO] all-score table             : {args.all_scores}", file=sys.stderr)
    if args.rejected:
        print(f"[INFO] rejected table              : {args.rejected}", file=sys.stderr)


if __name__ == "__main__":
    main()
'''

out_path = Path("/mnt/data/05_filter_probe_gc_complexity_dimer.py")
out_path.write_text(script, encoding="utf-8")
os.chmod(out_path, 0o755)

# compile check
res = subprocess.run([sys.executable, "-m", "py_compile", str(out_path)], capture_output=True, text=True, timeout=20)
print("Saved:", out_path)
print("py_compile:", "OK" if res.returncode == 0 else "FAILED")
if res.returncode != 0:
    print(res.stderr)

# also zip as scripts/design path
zip_path = Path("/mnt/data/05_filter_probe_gc_complexity_dimer_update.zip")
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
    z.write(out_path, Path("scripts/design/05_filter_probe_gc_complexity_dimer.py"))
print("Zip:", zip_path)
