#!/usr/bin/env bash
#SBATCH --job-name=eMito
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --output=eMito_%j.out
#SBATCH --error=eMito_%j.err
set -euo pipefail
# Submit from the repository root, or export EMITO_CODE_DIR before sbatch.
# SLURM executes a spool copy, so BASH_SOURCE is not the repository path.
CODE_DIR="${EMITO_CODE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
exec bash "$CODE_DIR/run_standard_panel.sh" "$@"
