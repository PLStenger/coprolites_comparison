#!/bin/bash

#SBATCH --job-name=07_mpa_tables_cop
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/07_mpa_tables.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/07_mpa_tables.out"

# ==============================================================================
# Script : 08_mpa_tables.sh
# Author : Pierre-Louis Stenger
# Purpose :
# 1. Convertir les rapports Kraken2 (.report) en tables MPA (.mpa) avec KrakenTools
# 2. Combiner les MPA par k-mer (k35, k29, k25) en tables TSV
#    avec en-têtes d'échantillons explicites (cop408, cop410, ...)
# ==============================================================================

# set -euo pipefail

module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics

# ==============================================================================
# PATHS
# ==============================================================================

BASE_DIR="/home/plstenge/coprolites_comparison"

# Répertoires des résultats Kraken2 par k-mer
K35_DIR="${BASE_DIR}/07_kraken2_k35"
K29_DIR="${BASE_DIR}/07_kraken2_k29"
K25_DIR="${BASE_DIR}/07_kraken2_k25"

# Répertoire d'installation / usage de KrakenTools
KRAKENTOOLS_DIR="${BASE_DIR}/08_krakentools/KrakenTools"

# Répertoires de sortie pour les tables MPA/TSV
OUT_BASE="${BASE_DIR}/09_mpa_tables"
K35_OUT="${OUT_BASE}/k35"
K29_OUT="${OUT_BASE}/k29"
K25_OUT="${OUT_BASE}/k25"

mkdir -p "${K35_OUT}" "${K29_OUT}" "${K25_OUT}"

# ==============================================================================
# STEP 0 : Installation / présence de KrakenTools
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
# STEP 1 : Fonction utilitaire pour convertir .report -> .mpa puis combiner
# ==============================================================================

convert_and_combine() {
    local IN_DIR="$1"       # répertoire des .report
    local OUT_DIR="$2"      # répertoire des .mpa/.tsv
    local LABEL="$3"        # label du k-mer (k35/k29/k25)

    echo ""
    echo "---------------------------------------------------------------"
    echo " Traitement des rapports Kraken2 pour ${LABEL}"
    echo "   IN  : ${IN_DIR}"
    echo "   OUT : ${OUT_DIR}"
    echo "---------------------------------------------------------------"

    mkdir -p "${OUT_DIR}"

    cd "${IN_DIR}" || exit 1

    declare -a MPA_FILES=()

    for report in *.report; do
        if [[ -f "${report}" ]]; then
            base=$(basename "${report}" .report)
            mpa_file="${OUT_DIR}/${base}.mpa"

            echo "  → Conversion MPA pour ${base}"
            python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" -r "${report}" -o "${mpa_file}"
            MPA_FILES+=("${mpa_file}")
        fi
    done

    if [[ ${#MPA_FILES[@]} -gt 0 ]]; then
        echo "  → Combinaison de tous les fichiers MPA pour ${LABEL}"
        python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" -i "${MPA_FILES[@]}" -o "${OUT_DIR}/combined_all.tsv"

        echo "  → Table combinée générée : ${OUT_DIR}/combined_all.tsv"
        echo "     Les en-têtes colonnes correspondent aux noms de fichiers .mpa,"
        echo "     donc aux échantillons (cop408, cop410, etc.)."
    else
        echo "  AUCUN fichier MPA généré pour ${LABEL} (pas de .report trouvés)."
    fi
}

# ==============================================================================
# STEP 2 : Conversion / combinaison pour chaque k-mer
# ==============================================================================

echo "======================================================================="
echo " STEP 2: Création des tables MPA/TSV pour chaque k-mer"
echo "======================================================================="

convert_and_combine "${K35_DIR}" "${K35_OUT}" "k35"
convert_and_combine "${K29_DIR}" "${K29_OUT}" "k29"
convert_and_combine "${K25_DIR}" "${K25_OUT}" "k25"

# ==============================================================================
# RÉSUMÉ FINAL
# ==============================================================================

echo ""
echo "======================================================================="
echo " ALL STEPS COMPLETED — MPA tables"
echo "======================================================================="
echo "  k35 MPA : ${K35_OUT}/*.mpa"
echo "  k35 TSV : ${K35_OUT}/combined_all.tsv"
echo "  k29 MPA : ${K29_OUT}/*.mpa"
echo "  k29 TSV : ${K29_OUT}/combined_all.tsv"
echo "  k25 MPA : ${K25_OUT}/*.mpa"
echo "  k25 TSV : ${K25_OUT}/combined_all.tsv"
echo "======================================================================="
