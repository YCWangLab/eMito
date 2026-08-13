#!/bin/bash
#SBATCH --job-name=eMito_taxa_AC
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G

set -Eeuo pipefail

source /path/to/miniconda/etc/profile.d/conda.sh
conda activate /path/to/emito/environment

# taxa-generate -> per-species access -> per-species collapse -> final merge
emito run \
  --metadata /path/to/metadata.tsv \
  --fasta-dir /path/to/fasta \
  --names-dmp /path/to/taxdump/names.dmp \
  --nodes-dmp /path/to/taxdump/nodes.dmp \
  --output-root /path/to/emito_taxa_access_collapse \
  --modes taxa \
  --taxa-access \
  --taxa-collapse \
  --threads "${SLURM_CPUS_PER_TASK:-4}"
