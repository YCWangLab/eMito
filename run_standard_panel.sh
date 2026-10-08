#!/usr/bin/env bash
# Run explicitly selected modules; individual modules never auto-chain.
set -euo pipefail
CODE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PIPELINE_PYTHON:-python3}"
MAFFT="${MAFFT_BIN:-mafft}"
METADATA="${1:?Usage: bash run_standard_panel.sh METADATA.tsv NEW_OUTPUT_DIRECTORY}"
OUT="${2:?Supply a new output directory}"
FASTA_DIR="${FASTA_DIR:?Set FASTA_DIR to the accession FASTA directory}"
NAMES_DMP="${NAMES_DMP:?Set NAMES_DMP to NCBI names.dmp}"
NODES_DMP="${NODES_DMP:?Set NODES_DMP to NCBI nodes.dmp}"
THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-4}}"
MIN_GENOME_LENGTH="${MIN_GENOME_LENGTH:-10000}"
[[ "$MIN_GENOME_LENGTH" =~ ^[0-9]+$ && "$MIN_GENOME_LENGTH" =~ [1-9] ]] || { echo 'ERROR: MIN_GENOME_LENGTH must be a positive integer in bp' >&2; exit 1; }
export PYTHONUNBUFFERED=1
PYTHON="$(command -v "$PYTHON")"
MAFFT="$(command -v "$MAFFT")"
[[ "$METADATA" = /* && "$OUT" = /* ]] || { echo 'ERROR: metadata/output paths must be absolute' >&2; exit 1; }
for f in "$METADATA" "$NAMES_DMP" "$NODES_DMP"; do
    [[ -r "$f" ]] || { echo "ERROR: missing/unreadable: $f" >&2; exit 1; }
done
[[ -d "$FASTA_DIR" && -x "$PYTHON" && -x "$MAFFT" ]] || { echo 'ERROR: check FASTA directory, Python and MAFFT paths' >&2; exit 1; }
"$PYTHON" -c 'import sys; assert sys.version_info >= (3,9), "Python >=3.9 required"'
[[ ! -e "$OUT" ]] || { echo "ERROR: output already exists; use a NEW directory: $OUT" >&2; exit 1; }
mkdir -p "$OUT"
exec > >(tee "$OUT/workflow.log") 2>&1
echo "Code: $CODE_DIR"
echo "Metadata: $METADATA"
echo "Output: $OUT"
echo "Minimum genome length (bp): $MIN_GENOME_LENGTH"
echo '[1/5] prepare: alignment ON, circular windows ON'
"$PYTHON" "$CODE_DIR/scripts/eMito-prepare.py" \
  --metadata "$METADATA" --fasta-dir "$FASTA_DIR" \
  --names-dmp "$NAMES_DMP" --nodes-dmp "$NODES_DMP" \
  --mafft "$MAFFT" --threads "$THREADS" --align \
  --min-genome-length "$MIN_GENOME_LENGTH" --window-length 52 --probe-step 5 --kmer-step 1 \
  --output "$OUT/01_prepare"
echo '[2/5] taxa-generate: species/genus AND eligible subgroups'
"$PYTHON" "$CODE_DIR/scripts/eMito-taxa-generate.py" \
  --prepared "$OUT/01_prepare" --representative-fraction 0.75 --small-group-all-max 3 \
  --output "$OUT/02_taxa_generate"
echo '[3/5] access: independently assess EACH species/subgroup FASTA; NO collapse'
"$PYTHON" "$CODE_DIR/scripts/eMito-access.py" \
  --inputs "$OUT/02_taxa_generate" --gc-min 35 --gc-max 65 \
  --complexity-min 0 --complexity-max 2 --dimer-k 11 --dimer 0.15 \
  --output "$OUT/03_taxa_access"
echo '[4/5] node-generate: family tiling; NO access; NO collapse'
"$PYTHON" "$CODE_DIR/scripts/eMito-node-generate.py" \
  --prepared "$OUT/01_prepare" --node-rank family --window-length 52 --node-step 5 \
  --output "$OUT/04_node_generate"
echo '[5/5] merge: assessed taxa/subgroups + unassessed family tiling; exact dedup ON'
"$PYTHON" "$CODE_DIR/scripts/eMito-merge.py" \
  --inputs "$OUT/03_taxa_access" "$OUT/04_node_generate" --dedup \
  --output "$OUT/05_final_merge"
echo "SUCCESS: $OUT/05_final_merge/probe_sets/merged.fasta"
echo "Counts: $OUT/05_final_merge/merge_summary.json"
