#!/bin/bash

#SBATCH --job-name=08_bracken_krona_mpa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -p smp
#SBATCH --mem=128G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/08_bracken_krona_mpa.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/08_bracken_krona_mpa.out"

# ==============================================================================
# Script : 09_bracken_krona_mpa_all_k.sh
# Author : Pierre-Louis Stenger
# Purpose :
# 1. Run Bracken on all Kraken2 reports for k35, k29 and k25
# 2. Keep separate output directories for each database size
# 3. Process merged and unmerged reports independently
# 4. Estimate abundances at species level (S) with read length r=50
# 5. Build Krona visualizations from Bracken reports
# 6. Export MPA files and combined TSV tables from Bracken reports
# ==============================================================================

# set -euo pipefail

module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics

# ==============================================================================
# PATHS
# ==============================================================================

BASE_DIR="/home/plstenge/coprolites_comparison"

#K35_KRAKEN_DIR="${BASE_DIR}/07_kraken2_k35"
K29_KRAKEN_DIR="${BASE_DIR}/07_kraken2_k29"
#K25_KRAKEN_DIR="${BASE_DIR}/07_kraken2_k25"

#K35_DB="/storage/groups/gdec/shared/Kraken_database/k2_core_nt_20250609"
K29_DB="/storage/groups/gdec/shared/Kraken_database/core_nt_k29"
#K25_DB="/storage/groups/gdec/shared/Kraken_database/core_nt_k25"

OUT_BASE="${BASE_DIR}/10_bracken"
#K35_OUT="${OUT_BASE}/k35"
K29_OUT="${OUT_BASE}/k29"
#K25_OUT="${OUT_BASE}/k25"

KRAKENTOOLS_DIR="${BASE_DIR}/08_krakentools/KrakenTools"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
READ_LEN=50
LEVEL="S"
THRESHOLD=10

#mkdir -p "${K35_OUT}" "${K29_OUT}" "${K25_OUT}"

# ==============================================================================
# STEP 0 : KrakenTools
# ==============================================================================

echo "======================================================================="
echo " STEP 0: Vérification/installation de KrakenTools"
echo "======================================================================="

if [[ ! -d "${KRAKENTOOLS_DIR}" ]]; then
    echo " KrakenTools non trouvé, installation..."
    mkdir -p "${BASE_DIR}/08_krakentools"
    cd "${BASE_DIR}/08_krakentools" || exit 1
    git clone https://github.com/jenniferlu717/KrakenTools.git
else
    echo " KrakenTools déjà présent dans : ${KRAKENTOOLS_DIR}"
fi

# ==============================================================================
# FUNCTION 1 : Bracken
# ==============================================================================

run_bracken_set() {
    local KRAKEN_DIR="$1"
    local DB_DIR="$2"
    local OUT_DIR="$3"
    local LABEL="$4"

    echo ""
    echo "======================================================================="
    echo " Bracken on ${LABEL}"
    echo "   Kraken dir : ${KRAKEN_DIR}"
    echo "   DB         : ${DB_DIR}"
    echo "   Output     : ${OUT_DIR}"
    echo "   Level      : ${LEVEL}"
    echo "   Read len   : ${READ_LEN}"
    echo "======================================================================="

    mkdir -p "${OUT_DIR}/merged"
    mkdir -p "${OUT_DIR}/unmerged"

    for report in "${KRAKEN_DIR}"/*_merged.report; do
        if [[ -f "${report}" ]]; then
            sample=$(basename "${report}" _merged.report)
            echo "  [merged] ${sample}"

            bracken \
                -d "${DB_DIR}" \
                -i "${report}" \
                -o "${OUT_DIR}/merged/${sample}_merged.bracken" \
                -w "${OUT_DIR}/merged/${sample}_merged.bracken.report" \
                -r "${READ_LEN}" \
                -l "${LEVEL}" \
                -t "${THRESHOLD}"
        fi
    done

    for report in "${KRAKEN_DIR}"/*_unmerged.report; do
        if [[ -f "${report}" ]]; then
            sample=$(basename "${report}" _unmerged.report)
            echo "  [unmerged] ${sample}"

            bracken \
                -d "${DB_DIR}" \
                -i "${report}" \
                -o "${OUT_DIR}/unmerged/${sample}_unmerged.bracken" \
                -w "${OUT_DIR}/unmerged/${sample}_unmerged.bracken.report" \
                -r "${READ_LEN}" \
                -l "${LEVEL}" \
                -t "${THRESHOLD}"
        fi
    done

    echo "  >>> Bracken termine pour ${LABEL} <<<"
}

# ==============================================================================
# FUNCTION 2 : Krona à partir des rapports Bracken
# ==============================================================================

run_krona_set() {
    local OUT_DIR="$1"
    local LABEL="$2"
    local KRONA_DIR="${OUT_DIR}/krona"

    mkdir -p "${KRONA_DIR}"

    echo ""
    echo "======================================================================="
    echo " Krona on ${LABEL}"
    echo "======================================================================="

    if ls "${OUT_DIR}/merged"/*.bracken.report >/dev/null 2>&1; then
        echo "  Generation du Krona pour les reads merged"
        ktImportTaxonomy \
            -t 5 \
            -m 3 \
            -o "${KRONA_DIR}/krona_merged.html" \
            "${OUT_DIR}/merged"/*.bracken.report
    else
        echo "  Aucun fichier merged .bracken.report trouvé pour ${LABEL}"
    fi

    if ls "${OUT_DIR}/unmerged"/*.bracken.report >/dev/null 2>&1; then
        echo "  Generation du Krona pour les reads unmerged"
        ktImportTaxonomy \
            -t 5 \
            -m 3 \
            -o "${KRONA_DIR}/krona_unmerged.html" \
            "${OUT_DIR}/unmerged"/*.bracken.report
    else
        echo "  Aucun fichier unmerged .bracken.report trouvé pour ${LABEL}"
    fi

    echo "  >>> Krona termine pour ${LABEL} <<<"
}

# ==============================================================================
# FUNCTION 3 : Export MPA + TSV combinés
# ==============================================================================

export_mpa_set() {
    local OUT_DIR="$1"
    local LABEL="$2"
    local MPA_DIR="${OUT_DIR}/mpa_tables"
    local MPA_MERGED_DIR="${MPA_DIR}/merged"
    local MPA_UNMERGED_DIR="${MPA_DIR}/unmerged"

    mkdir -p "${MPA_MERGED_DIR}" "${MPA_UNMERGED_DIR}"

    echo ""
    echo "======================================================================="
    echo " Export MPA/TSV for ${LABEL}"
    echo "======================================================================="

    cd "${OUT_DIR}" || exit 1

    declare -a merged_mpa_files=()
    declare -a unmerged_mpa_files=()

    for report in merged/*.bracken.report; do
        if [[ -f "${report}" ]]; then
            base=$(basename "${report}" .bracken.report)
            mpa_file="${MPA_MERGED_DIR}/${base}.mpa"
            echo "  → Conversion MPA merged : ${base}"
            python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" -r "${report}" -o "${mpa_file}"
            merged_mpa_files+=("${mpa_file}")
        fi
    done

    for report in unmerged/*.bracken.report; do
        if [[ -f "${report}" ]]; then
            base=$(basename "${report}" .bracken.report)
            mpa_file="${MPA_UNMERGED_DIR}/${base}.mpa"
            echo "  → Conversion MPA unmerged : ${base}"
            python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" -r "${report}" -o "${mpa_file}"
            unmerged_mpa_files+=("${mpa_file}")
        fi
    done

    if [[ ${#merged_mpa_files[@]} -gt 0 ]]; then
        echo "  → Combinaison TSV merged"
        python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" -i "${merged_mpa_files[@]}" -o "${MPA_MERGED_DIR}/combined_merged.tsv"
    else
        echo "  Aucun fichier MPA merged à combiner pour ${LABEL}"
    fi

    if [[ ${#unmerged_mpa_files[@]} -gt 0 ]]; then
        echo "  → Combinaison TSV unmerged"
        python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" -i "${unmerged_mpa_files[@]}" -o "${MPA_UNMERGED_DIR}/combined_unmerged.tsv"
    else
        echo "  Aucun fichier MPA unmerged à combiner pour ${LABEL}"
    fi

    echo "  >>> Export MPA termine pour ${LABEL} <<<"
}

# ==============================================================================
# RUN
# ==============================================================================

echo "======================================================================="
echo " STEP 1: Running Bracken on all Kraken result sets"
echo "======================================================================="
echo " Parameters: read length=${READ_LEN}, level=${LEVEL}, threshold=${THRESHOLD}"

run_bracken_set "${K35_KRAKEN_DIR}" "${K35_DB}" "${K35_OUT}" "k35"
run_bracken_set "${K29_KRAKEN_DIR}" "${K29_DB}" "${K29_OUT}" "k29"
run_bracken_set "${K25_KRAKEN_DIR}" "${K25_DB}" "${K25_OUT}" "k25"

echo ""
echo "======================================================================="
echo " STEP 2: Krona from Bracken reports"
echo "======================================================================="

run_krona_set "${K35_OUT}" "k35"
run_krona_set "${K29_OUT}" "k29"
run_krona_set "${K25_OUT}" "k25"

echo ""
echo "======================================================================="
echo " STEP 3: Export MPA and combined TSV tables"
echo "======================================================================="

export_mpa_set "${K35_OUT}" "k35"
export_mpa_set "${K29_OUT}" "k29"
export_mpa_set "${K25_OUT}" "k25"

# ==============================================================================
# SUMMARY
# ==============================================================================

echo ""
echo "======================================================================="
echo " ALL STEPS COMPLETED — Bracken + Krona + MPA"
echo "======================================================================="
echo "  k35 outputs : ${K35_OUT}"
echo "  k29 outputs : ${K29_OUT}"
echo "  k25 outputs : ${K25_OUT}"
echo "  Krona       : ${OUT_BASE}/k*/krona/"
echo "  MPA tables  : ${OUT_BASE}/k*/mpa_tables/"
echo "======================================================================="
