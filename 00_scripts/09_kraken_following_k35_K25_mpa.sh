#!/bin/bash
#SBATCH --job-name=mpa_k35_k25
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p smp
#SBATCH --mem=20G
#SBATCH --time=02:00:00
#SBATCH --output=/home/plstenge/coprolites_comparison/00_scripts/mpa_k35_k25.out
#SBATCH --error=/home/plstenge/coprolites_comparison/00_scripts/mpa_k35_k25.err

# ============================================================
# Création de tables MPA à partir de :
#   - k35 : reads merged complets
#   - k25 : reads non assignés après k35 ET k29
#
# ATTENTION :
# Les résultats k29 sont volontairement exclus de la table finale.
# ============================================================

WORKDIR="/home/plstenge/coprolites_comparison"

KRAKENTOOLS_DIR="${WORKDIR}/08_krakentools/KrakenTools"

K35_DIR="${WORKDIR}/09_kraken_following/k35"
K25_DIR="${WORKDIR}/09_kraken_following/k25"

OUTDIR="${WORKDIR}/09_kraken_following/mpa_k35_k25"
MPA_K35="${OUTDIR}/k35"
MPA_K25="${OUTDIR}/k25"
MPA_PER_SAMPLE="${OUTDIR}/per_sample_combined"

SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")

mkdir -p "${MPA_K35}" "${MPA_K25}" "${MPA_PER_SAMPLE}"

# Vérifie KrakenTools
if [[ ! -f "${KRAKENTOOLS_DIR}/kreport2mpa.py" ]]; then
    echo "ERREUR : kreport2mpa.py introuvable : ${KRAKENTOOLS_DIR}"
    exit 1
fi

if [[ ! -f "${KRAKENTOOLS_DIR}/combine_mpa.py" ]]; then
    echo "ERREUR : combine_mpa.py introuvable : ${KRAKENTOOLS_DIR}"
    exit 1
fi

echo "============================================================"
echo "STEP 1 — Conversion des rapports k35 et k25 en MPA"
echo "============================================================"

for SAMPLE in "${SAMPLES[@]}"; do

    # k35 : classification initiale sur tous les reads merged
    K35_REPORT="${K35_DIR}/${SAMPLE}_merged.report"
    K35_MPA="${MPA_K35}/${SAMPLE}_merged_k35.mpa"

    if [[ -s "${K35_REPORT}" ]]; then
        echo "[k35] ${SAMPLE}"
        python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" \
            -r "${K35_REPORT}" \
            -o "${K35_MPA}"
    else
        echo "[WARN] Rapport k35 absent/vide : ${K35_REPORT}"
    fi

    # k25 : classification sur reads non assignés par k35 ET k29
    K25_REPORT="${K25_DIR}/${SAMPLE}_merged_unassigned_k29.report"
    K25_MPA="${MPA_K25}/${SAMPLE}_merged_unassigned_k29_k25.mpa"

    if [[ -s "${K25_REPORT}" ]]; then
        echo "[k25] ${SAMPLE}"
        python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" \
            -r "${K25_REPORT}" \
            -o "${K25_MPA}"
    else
        echo "[WARN] Rapport k25 absent/vide : ${K25_REPORT}"
    fi
done

echo "============================================================"
echo "STEP 2 — Somme des comptages k35 + k25 par échantillon"
echo "============================================================"

MPA_K35_ENV="${MPA_K35}" \
MPA_K25_ENV="${MPA_K25}" \
MPA_PER_SAMPLE_ENV="${MPA_PER_SAMPLE}" \
python3 << 'PYEOF'
import os

samples = ["cop408", "cop410", "cop412", "cop414", "cop417"]

mpa_k35_dir = os.environ["MPA_K35_ENV"]
mpa_k25_dir = os.environ["MPA_K25_ENV"]
out_dir = os.environ["MPA_PER_SAMPLE_ENV"]


def load_mpa(path):
    counts = {}

    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return counts

    with open(path, "r") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 2:
                continue

            taxon = fields[0]

            try:
                count = float(fields[1])
            except ValueError:
                continue

            counts[taxon] = counts.get(taxon, 0) + count

    return counts


for sample in samples:
    k35_mpa = os.path.join(
        mpa_k35_dir,
        f"{sample}_merged_k35.mpa"
    )

    k25_mpa = os.path.join(
        mpa_k25_dir,
        f"{sample}_merged_unassigned_k29_k25.mpa"
    )

    combined = {}

    for infile in [k35_mpa, k25_mpa]:
        for taxon, count in load_mpa(infile).items():
            combined[taxon] = combined.get(taxon, 0) + count

    if not combined:
        print(f"[WARN] Pas de données MPA pour {sample}.")
        continue

    outfile = os.path.join(
        out_dir,
        f"{sample}_k35_k25_only.mpa"
    )

    with open(outfile, "w") as out:
        for taxon in sorted(combined):
            count = combined[taxon]
            out.write(f"{taxon}\t{int(count) if count.is_integer() else count}\n")

    print(f"[OK] {sample}: {outfile}")

PYEOF

echo "============================================================"
echo "STEP 3 — Table finale tous échantillons"
echo "============================================================"

FINAL_MPA_FILES=()

for SAMPLE in "${SAMPLES[@]}"; do
    FILE="${MPA_PER_SAMPLE}/${SAMPLE}_k35_k25_only.mpa"

    if [[ -s "${FILE}" ]]; then
        FINAL_MPA_FILES+=("${FILE}")
    fi
done

if [[ "${#FINAL_MPA_FILES[@]}" -eq 0 ]]; then
    echo "ERREUR : aucune table MPA combinée disponible."
    exit 1
fi

python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" \
    -i "${FINAL_MPA_FILES[@]}" \
    -o "${OUTDIR}/combined_all_samples_k35_k25_only.tsv"

echo "============================================================"
echo "Terminé."
echo "Table finale :"
echo "${OUTDIR}/combined_all_samples_k35_k25_only.tsv"
echo "============================================================"
