#!/bin/bash
#SBATCH --job-name=mpa_k35_k29
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p smp
#SBATCH --mem=20G
#SBATCH --time=02:00:00
#SBATCH --output=/home/plstenge/coprolites_comparison/00_scripts/09_kraken_following_k35_29_mpa.out
#SBATCH --error=/home/plstenge/coprolites_comparison/00_scripts/09_kraken_following_k35_29_mpa.err

# ============================================================
# Génération des tables MPA à partir de :
# - k35 : analyses initiales sur tous les reads merged
# - k29 : analyses des reads non assignés par k35
#
# AUCUNE donnée k25 n'est utilisée.
# ============================================================

WORKDIR="/home/plstenge/coprolites_comparison"
KRAKENTOOLS_DIR="${WORKDIR}/08_krakentools/KrakenTools"

K35_DIR="${WORKDIR}/09_kraken_following/k35"
K29_DIR="${WORKDIR}/09_kraken_following/k29"

OUTDIR="${WORKDIR}/09_kraken_following/mpa_k35_k29"
MPA_K35="${OUTDIR}/k35"
MPA_K29="${OUTDIR}/k29"
MPA_PER_SAMPLE="${OUTDIR}/per_sample_combined"

SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")

mkdir -p "${MPA_K35}" "${MPA_K29}" "${MPA_PER_SAMPLE}"

# Vérification de KrakenTools
if [[ ! -f "${KRAKENTOOLS_DIR}/kreport2mpa.py" ]]; then
    echo "ERREUR : kreport2mpa.py introuvable : ${KRAKENTOOLS_DIR}"
    exit 1
fi

if [[ ! -f "${KRAKENTOOLS_DIR}/combine_mpa.py" ]]; then
    echo "ERREUR : combine_mpa.py introuvable : ${KRAKENTOOLS_DIR}"
    exit 1
fi

echo "============================================================"
echo "STEP 1 — Conversion des rapports Kraken vers MPA"
echo "============================================================"

for SAMPLE in "${SAMPLES[@]}"; do

    # k35 : rapport sur les reads merged complets
    K35_REPORT="${K35_DIR}/${SAMPLE}_merged.report"
    K35_MPA="${MPA_K35}/${SAMPLE}_merged.mpa"

    if [[ -s "${K35_REPORT}" ]]; then
        echo "[k35] Conversion : ${K35_REPORT}"
        python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" \
            -r "${K35_REPORT}" \
            -o "${K35_MPA}"
    else
        echo "[WARN] Rapport k35 absent/vide : ${K35_REPORT}"
    fi

    # k29 : rapport portant exclusivement sur les non-assignés de k35
    K29_REPORT="${K29_DIR}/${SAMPLE}_merged_unassigned_k35.report"
    K29_MPA="${MPA_K29}/${SAMPLE}_merged_unassigned_k35.mpa"

    if [[ -s "${K29_REPORT}" ]]; then
        echo "[k29] Conversion : ${K29_REPORT}"
        python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" \
            -r "${K29_REPORT}" \
            -o "${K29_MPA}"
    else
        echo "[WARN] Rapport k29 absent/vide : ${K29_REPORT}"
    fi
done

echo "============================================================"
echo "STEP 2 — Fusion k35 + k29 pour chaque échantillon"
echo "============================================================"

MPA_K35_ENV="${MPA_K35}" \
MPA_K29_ENV="${MPA_K29}" \
MPA_PER_SAMPLE_ENV="${MPA_PER_SAMPLE}" \
python3 << 'PYEOF'
import os

samples = ["cop408", "cop410", "cop412", "cop414", "cop417"]

mpa_k35 = os.environ["MPA_K35_ENV"]
mpa_k29 = os.environ["MPA_K29_ENV"]
out_dir = os.environ["MPA_PER_SAMPLE_ENV"]


def load_mpa(path):
    """Lit une table taxon<TAB>compte et retourne un dictionnaire."""
    counts = {}

    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return counts

    with open(path, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) < 2:
                continue

            taxon = fields[0]

            try:
                abundance = float(fields[1])
            except ValueError:
                continue

            counts[taxon] = counts.get(taxon, 0) + abundance

    return counts


for sample in samples:

    k35_file = os.path.join(mpa_k35, f"{sample}_merged.mpa")
    k29_file = os.path.join(mpa_k29, f"{sample}_merged_unassigned_k35.mpa")

    combined = {}

    for mpa_file in [k35_file, k29_file]:
        counts = load_mpa(mpa_file)

        for taxon, abundance in counts.items():
            combined[taxon] = combined.get(taxon, 0) + abundance

    if not combined:
        print(f"[WARN] Aucune donnée MPA trouvée pour {sample}.")
        continue

    output_file = os.path.join(out_dir, f"{sample}_k35_k29_following.mpa")

    with open(output_file, "w") as out:
        for taxon in sorted(combined):
            abundance = combined[taxon]

            # Écrit des entiers sans '.0' si les comptages sont entiers
            if abundance.is_integer():
                abundance = int(abundance)

            out.write(f"{taxon}\t{abundance}\n")

    print(f"[OK] {sample}: {output_file}")

PYEOF

echo "============================================================"
echo "STEP 3 — Table finale combinant les cinq échantillons"
echo "============================================================"

FINAL_MPA_FILES=()

for SAMPLE in "${SAMPLES[@]}"; do
    FILE="${MPA_PER_SAMPLE}/${SAMPLE}_k35_k29_following.mpa"

    if [[ -s "${FILE}" ]]; then
        FINAL_MPA_FILES+=("${FILE}")
    fi
done

if [[ "${#FINAL_MPA_FILES[@]}" -gt 0 ]]; then
    python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" \
        -i "${FINAL_MPA_FILES[@]}" \
        -o "${OUTDIR}/combined_all_samples_k35_k29.tsv"

    echo "[OK] Table finale créée :"
    echo "     ${OUTDIR}/combined_all_samples_k35_k29.tsv"
else
    echo "[ERREUR] Aucune table MPA par échantillon à combiner."
    exit 1
fi

echo "============================================================"
echo "Terminé : MPA k35 + k29 uniquement, sans k25."
echo "Dossier de sortie : ${OUTDIR}"
echo "============================================================"
