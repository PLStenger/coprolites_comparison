#!/bin/bash
#SBATCH --job-name=mpa_k35
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p smp
#SBATCH --mem=10G
#SBATCH --time=01:00:00
#SBATCH --output=/home/plstenge/coprolites_comparison/00_scripts/mpa_k35.out
#SBATCH --error=/home/plstenge/coprolites_comparison/00_scripts/mpa_k35.err

set -euo pipefail

# ============================================================
# Conversion des rapports Kraken2 k35 uniquement en tables MPA
# ============================================================

WORKDIR="/home/plstenge/coprolites_comparison"

KRAKENTOOLS_DIR="${WORKDIR}/08_krakentools/KrakenTools"
K35_DIR="${WORKDIR}/09_kraken_following/k35"

# Nouveau dossier, pour ne pas écraser les résultats k35+k29
OUTDIR="${WORKDIR}/09_kraken_following/mpa_k35_only"

SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")

mkdir -p "${OUTDIR}"

# Vérifications
if [[ ! -f "${KRAKENTOOLS_DIR}/kreport2mpa.py" ]]; then
    echo "ERREUR : kreport2mpa.py introuvable dans ${KRAKENTOOLS_DIR}"
    exit 1
fi

if [[ ! -f "${KRAKENTOOLS_DIR}/combine_mpa.py" ]]; then
    echo "ERREUR : combine_mpa.py introuvable dans ${KRAKENTOOLS_DIR}"
    exit 1
fi

echo "============================================================"
echo "Conversion des rapports Kraken2 k35 en tables MPA"
echo "============================================================"

MPA_FILES=()

for SAMPLE in "${SAMPLES[@]}"; do

    REPORT="${K35_DIR}/${SAMPLE}_merged.report"
    MPA="${OUTDIR}/${SAMPLE}_k35.mpa"

    if [[ ! -s "${REPORT}" ]]; then
        echo "[WARN] Rapport absent ou vide : ${REPORT}"
        continue
    fi

    echo "[k35] Conversion : ${REPORT}"

    python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" \
        -r "${REPORT}" \
        -o "${MPA}"

    if [[ -s "${MPA}" ]]; then
        MPA_FILES+=("${MPA}")
        echo "[OK] Table créée : ${MPA}"
    else
        echo "[WARN] Table MPA vide : ${MPA}"
    fi
done

echo "============================================================"
echo "Fusion des tables MPA k35 des cinq échantillons"
echo "============================================================"

if [[ "${#MPA_FILES[@]}" -eq 0 ]]; then
    echo "ERREUR : aucune table MPA k35 créée."
    exit 1
fi

python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" \
    -i "${MPA_FILES[@]}" \
    -o "${OUTDIR}/combined_all_samples_k35.tsv"

echo "============================================================"
echo "Terminé."
echo "Tables individuelles : ${OUTDIR}/copXXX_k35.mpa"
echo "Table combinée finale :"
echo "${OUTDIR}/combined_all_samples_k35.tsv"
echo "============================================================"
