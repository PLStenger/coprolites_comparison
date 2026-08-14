#!/bin/bash
#SBATCH --job-name=mpa_k29_k25
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p smp
#SBATCH --mem=10G
#SBATCH --time=01:00:00
#SBATCH --output=/home/plstenge/coprolites_comparison/00_scripts/mpa_k29_k25.out
#SBATCH --error=/home/plstenge/coprolites_comparison/00_scripts/mpa_k29_k25.err

# ============================================================
# Création de tables MPA séparées :
#   1. k29 uniquement
#   2. k25 uniquement
#
# Aucune fusion ni somme entre k29, k25 et k35.
# ============================================================

WORKDIR="/home/plstenge/coprolites_comparison"

KRAKENTOOLS_DIR="${WORKDIR}/08_krakentools/KrakenTools"

K29_DIR="${WORKDIR}/09_kraken_following/k29"
K25_DIR="${WORKDIR}/09_kraken_following/k25"

OUTDIR="${WORKDIR}/09_kraken_following/mpa_k29_k25_only"
OUT_K29="${OUTDIR}/k29"
OUT_K25="${OUTDIR}/k25"

SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")

mkdir -p "${OUT_K29}" "${OUT_K25}"

# Vérification de KrakenTools
if [[ ! -f "${KRAKENTOOLS_DIR}/kreport2mpa.py" ]]; then
    echo "ERREUR : kreport2mpa.py introuvable dans ${KRAKENTOOLS_DIR}"
    exit 1
fi

if [[ ! -f "${KRAKENTOOLS_DIR}/combine_mpa.py" ]]; then
    echo "ERREUR : combine_mpa.py introuvable dans ${KRAKENTOOLS_DIR}"
    exit 1
fi

echo "============================================================"
echo "STEP 1 — Conversion des rapports k29 en MPA"
echo "============================================================"

MPA_K29_FILES=()

for SAMPLE in "${SAMPLES[@]}"; do

    REPORT="${K29_DIR}/${SAMPLE}_merged_unassigned_k35.report"
    MPA="${OUT_K29}/${SAMPLE}_k29_only.mpa"

    if [[ ! -s "${REPORT}" ]]; then
        echo "[WARN] Rapport k29 absent ou vide : ${REPORT}"
        continue
    fi

    echo "[k29] Conversion : ${REPORT}"

    python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" \
        -r "${REPORT}" \
        -o "${MPA}"

    if [[ -s "${MPA}" ]]; then
        MPA_K29_FILES+=("${MPA}")
    else
        echo "[WARN] MPA k29 vide : ${MPA}"
    fi
done

echo "============================================================"
echo "STEP 2 — Table combinée : k29 uniquement"
echo "============================================================"

if [[ "${#MPA_K29_FILES[@]}" -gt 0 ]]; then
    python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" \
        -i "${MPA_K29_FILES[@]}" \
        -o "${OUTDIR}/combined_all_samples_k29_only.tsv"

    echo "[OK] Table k29 : ${OUTDIR}/combined_all_samples_k29_only.tsv"
else
    echo "[WARN] Aucune table k29 disponible."
fi

echo "============================================================"
echo "STEP 3 — Conversion des rapports k25 en MPA"
echo "============================================================"

MPA_K25_FILES=()

for SAMPLE in "${SAMPLES[@]}"; do

    REPORT="${K25_DIR}/${SAMPLE}_merged_unassigned_k29.report"
    MPA="${OUT_K25}/${SAMPLE}_k25_only.mpa"

    if [[ ! -s "${REPORT}" ]]; then
        echo "[WARN] Rapport k25 absent ou vide : ${REPORT}"
        continue
    fi

    echo "[k25] Conversion : ${REPORT}"

    python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" \
        -r "${REPORT}" \
        -o "${MPA}"

    if [[ -s "${MPA}" ]]; then
        MPA_K25_FILES+=("${MPA}")
    else
        echo "[WARN] MPA k25 vide : ${MPA}"
    fi
done

echo "============================================================"
echo "STEP 4 — Table combinée : k25 uniquement"
echo "============================================================"

if [[ "${#MPA_K25_FILES[@]}" -gt 0 ]]; then
    python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" \
        -i "${MPA_K25_FILES[@]}" \
        -o "${OUTDIR}/combined_all_samples_k25_only.tsv"

    echo "[OK] Table k25 : ${OUTDIR}/combined_all_samples_k25_only.tsv"
else
    echo "[WARN] Aucune table k25 disponible."
fi

echo "============================================================"
echo "Terminé."
echo "MPA k29 par échantillon : ${OUT_K29}/"
echo "MPA k25 par échantillon : ${OUT_K25}/"
echo "Table
