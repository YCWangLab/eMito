#!/bin/bash
#SBATCH --job-name=eMito
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G

set -Eeuo pipefail

# Activate the environment in which eMito and MAFFT are installed.
source /path/to/miniconda/etc/profile.d/conda.sh
conda activate /path/to/emito/environment

emito run \
  --metadata /path/to/metadata.tsv \
  --fasta-dir /path/to/fasta \
  --names-dmp /path/to/taxdump/names.dmp \
  --nodes-dmp /path/to/taxdump/nodes.dmp \
  --output-root /path/to/emito_output \
  --threads "${SLURM_CPUS_PER_TASK:-4}"
