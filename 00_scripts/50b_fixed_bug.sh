#!/bin/bash

#SBATCH --job-name=rerun_bracken_k25_k29
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=36
#SBATCH --mem=200G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/50b_fixed_bug.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/50b_fixed_bug.out"

set -euo pipefail

# ==============================================================================
# Script : rerun_bracken_k25_k29.sh
# Purpose :
# Relancer UNIQUEMENT Bracken sur les reports Kraken2 deja existants
# pour k25 et k29, avec les MEMES noms de sorties que dans le pipeline
# principal afin d'ecraser/remplacer proprement les anciens fichiers.
# ==============================================================================

# ==============================================================================
# ENVIRONMENT
# ==============================================================================

module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics

# ==============================================================================
# RUNS / SAMPLES (identiques au pipeline principal)
# ==============================================================================

RUN_LABELS=("run1" "run2" "run3" "run4" "run5")
SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")

# ==============================================================================
# PATHS (identiques au pipeline principal)
# ==============================================================================

WORKDIR="/home/plstenge/coprolites_comparison"
OUT_ROOT="${WORKDIR}/11_per_run_analysis"

KRAKEN_K25_DIR="${OUT_ROOT}/05_kraken2_k25"
KRAKEN_K29_DIR="${OUT_ROOT}/05_kraken2_k29"

BRACKEN_DIR="${OUT_ROOT}/06_bracken"

# ==============================================================================
# DATABASES / PARAMS (identiques au pipeline principal)
# ==============================================================================

KRAKEN2_DB_K25="/storage/groups/gdec/shared/Kraken_database/core_nt_k25"
KRAKEN2_DB_K29="/storage/groups/gdec/shared/Kraken_database/core_nt_k29"

BRACKEN_READ_LEN=50
BRACKEN_LEVEL="S"
BRACKEN_THRESHOLD=10

# ==============================================================================
# OUTPUT DIRS
# ==============================================================================

mkdir -p "${BRACKEN_DIR}/k25" "${BRACKEN_DIR}/k29"

# ==============================================================================
# FUNCTION
# ==============================================================================

run_bracken_for_kmer() {
  local KRAKEN_ROOT="$1"
  local DB="$2"
  local OUT_ROOT_K="$3"
  local KLABEL="$4"

  echo "======================================================================="
  echo " Relance Bracken pour ${KLABEL}"
  echo " Kraken root : ${KRAKEN_ROOT}"
  echo " Output root : ${OUT_ROOT_K}"
  echo "======================================================================="

  for RUN_LABEL in "${RUN_LABELS[@]}"; do
    for SAMPLE in "${SAMPLES[@]}"; do
      UNIT="${SAMPLE}_${RUN_LABEL}"
      KRAKEN_DIR="${KRAKEN_ROOT}/${RUN_LABEL}/${SAMPLE}"
      OUT_DIR="${OUT_ROOT_K}/${RUN_LABEL}/${SAMPLE}"
      mkdir -p "${OUT_DIR}"

      MERGED_REPORT="${KRAKEN_DIR}/${UNIT}_merged.report"
      UNMERGED_REPORT="${KRAKEN_DIR}/${UNIT}_unmerged.report"

      if [[ -f "${MERGED_REPORT}" ]]; then
        echo " [Bracken ${KLABEL}] ${UNIT} (merged)"
        bracken \
          -d "${DB}" \
          -i "${MERGED_REPORT}" \
          -o "${OUT_DIR}/${UNIT}_merged.bracken" \
          -w "${OUT_DIR}/${UNIT}_merged.bracken.report" \
          -r "${BRACKEN_READ_LEN}" \
          -l "${BRACKEN_LEVEL}" \
          -t "${BRACKEN_THRESHOLD}"
      else
        echo " [SKIP ${KLABEL}] ${UNIT} (merged) : report absent -> ${MERGED_REPORT}"
      fi

      if [[ -f "${UNMERGED_REPORT}" ]]; then
        echo " [Bracken ${KLABEL}] ${UNIT} (unmerged)"
        bracken \
          -d "${DB}" \
          -i "${UNMERGED_REPORT}" \
          -o "${OUT_DIR}/${UNIT}_unmerged.bracken" \
          -w "${OUT_DIR}/${UNIT}_unmerged.bracken.report" \
          -r "${BRACKEN_READ_LEN}" \
          -l "${BRACKEN_LEVEL}" \
          -t "${BRACKEN_THRESHOLD}"
      else
        echo " [SKIP ${KLABEL}] ${UNIT} (unmerged) : report absent -> ${UNMERGED_REPORT}"
      fi
    done
  done
}

# ==============================================================================
# RUN
# ==============================================================================

run_bracken_for_kmer "${KRAKEN_K25_DIR}" "${KRAKEN2_DB_K25}" "${BRACKEN_DIR}/k25" "k25"
run_bracken_for_kmer "${KRAKEN_K29_DIR}" "${KRAKEN2_DB_K29}" "${BRACKEN_DIR}/k29" "k29"

echo ""
echo "======================================================================="
echo " Relance Bracken terminee pour k25 et k29"
echo " Les fichiers existants ont ete remplaces avec les memes noms."
echo "======================================================================="
