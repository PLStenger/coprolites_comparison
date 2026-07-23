#!/bin/bash

#SBATCH --job-name=08_bracken_create_db
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -p smp
#SBATCH --mem=128G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/08_bracken_create_db.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/08_bracken_create_db.out"

module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics


# k29
bracken-build \
  -d /storage/groups/gdec/shared/Kraken_database/core_nt_k29 \
  -k 29 \
  -l 50 \
  -t 16

# k25
bracken-build \
  -d /storage/groups/gdec/shared/Kraken_database/core_nt_k25 \
  -k 25 \
  -l 50 \
  -t 16
