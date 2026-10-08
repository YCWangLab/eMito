#!/usr/bin/env bash
#SBATCH --job-name=eMito_compare
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --output=eMito_compare_%j.out
#SBATCH --error=eMito_compare_%j.err
set -euo pipefail
CODE_DIR="${EMITO_CODE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
exec "${PIPELINE_PYTHON:-python3}" "$CODE_DIR/scripts/run_matched_comparison.py" \
  --threads "${SLURM_CPUS_PER_TASK:-4}" "$@"
