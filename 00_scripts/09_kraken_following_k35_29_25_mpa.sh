#!/bin/bash

#SBATCH --job-name=09_kraken_following
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/09_kraken_following.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/09_kraken_following.out"

# ==============================================================================
# Script : 09_kraken_following_k35_29_25_mpa.sh
# Author : Pierre-Louis Stenger
# Purpose :
#   Pipeline "kraken following" (cascade k35 -> k29 -> k25) applique
#   UNIQUEMENT aux echantillons "grouped" (cop408, cop410, cop412, cop414,
#   cop417) et UNIQUEMENT sur les reads mergés issus de fastp
#   (clean_<sample>_grouped_dedup_fastp_merged.fastq.gz).
#
#   Logique de cascade :
#     1) Classification Kraken2 avec la base k35 (la plus stricte / plus
#        sensible taxonomiquement) sur les reads merged de chaque echantillon.
#     2) Les reads NON assignes (unclassified, taxid=0) par k35 sont extraits
#        et RE-analyses (et seulement eux) avec la base k29.
#     3) Les reads encore NON assignes par k29 sont extraits et RE-analyses
#        (et seulement eux) avec la base k25.
#     4) A la fin, chaque read merged de chaque echantillon a ete tente
#        d'assignation successivement par k35, puis k29 (sur les restes),
#        puis k25 (sur les restes des restes).
#
#   Sorties produites :
#     - 09_kraken_following/k35/<sample>_merged.kraken(.report)
#     - 09_kraken_following/k29/<sample>_merged_unassigned_k35.kraken(.report)
#     - 09_kraken_following/k25/<sample>_merged_unassigned_k29.kraken(.report)
#     - 09_kraken_following/mpa/k35|k29|k25/<sample>.mpa
#     - 09_kraken_following/mpa/per_sample_combined/<sample>_following.mpa
#         (somme des comptages k35+k29+k25 pour chaque echantillon = la
#          table mpa "finale" de la cascade pour cet echantillon)
#     - 09_kraken_following/mpa/combined_all_samples_following.tsv
#         (table finale groupant TOUS les echantillons "grouped", cascade
#          complete k35+k29+k25)
#     - 09_kraken_following/combined_kraken_for_mapdamage/<sample>_following.kraken
#         (concatenation des 3 fichiers .kraken bruts k35+k29+k25 : un seul
#          fichier par echantillon contenant l'assignation taxonomique finale
#          de CHAQUE read merged, quel que soit le niveau de k-mer qui l'a
#          assigne. Ce fichier + le fastq merged original peuvent etre
#          utilises directement par extract_kraken_reads.py (KrakenTools)
#          pour recuperer TOUS les reads d'une espece donnee -- assignes
#          par k35, k29 OU k25 -- en une seule extraction, avant de lancer
#          mapDamage, comme dans 52_mapdamage.sh)
# ==============================================================================

module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics

# ==============================================================================
# PATHS & PARAMETRES
# ==============================================================================

WORKDIR="/home/plstenge/coprolites_comparison"

ENTREE="${WORKDIR}/06_fastp"

KRAKEN2_DB_K35="/storage/groups/gdec/shared/Kraken_database/k2_core_nt_20250609"
KRAKEN2_DB_K29="/storage/groups/gdec/shared/Kraken_database/core_nt_k29"
KRAKEN2_DB_K25="/storage/groups/gdec/shared/Kraken_database/core_nt_k25"

SORTIE="${WORKDIR}/09_kraken_following"
SORTIE_K35="${SORTIE}/k35"
SORTIE_K29="${SORTIE}/k29"
SORTIE_K25="${SORTIE}/k25"
SORTIE_TMP="${SORTIE}/tmp_unassigned_fastq"
SORTIE_MPA="${SORTIE}/mpa"
SORTIE_MPA_K35="${SORTIE_MPA}/k35"
SORTIE_MPA_K29="${SORTIE_MPA}/k29"
SORTIE_MPA_K25="${SORTIE_MPA}/k25"
SORTIE_MPA_PER_SAMPLE="${SORTIE_MPA}/per_sample_combined"
SORTIE_COMBINED_KRAKEN="${SORTIE}/combined_kraken_for_mapdamage"

KRAKENTOOLS_DIR="${WORKDIR}/08_krakentools/KrakenTools"

THREADS="${SLURM_CPUS_PER_TASK:-36}"

# Echantillons "grouped" uniquement
SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")

mkdir -p "${SORTIE_K35}" "${SORTIE_K29}" "${SORTIE_K25}" "${SORTIE_TMP}" \
         "${SORTIE_MPA_K35}" "${SORTIE_MPA_K29}" "${SORTIE_MPA_K25}" \
         "${SORTIE_MPA_PER_SAMPLE}" "${SORTIE_COMBINED_KRAKEN}"

# ==============================================================================
# STEP 0 : Verification / installation de KrakenTools
# ==============================================================================

echo "======================================================================="
echo " STEP 0: Verification/installation de KrakenTools"
echo "======================================================================="

if [[ ! -d "${KRAKENTOOLS_DIR}" ]]; then
    echo " KrakenTools non trouve, installation..."
    mkdir -p "${WORKDIR}/08_krakentools"
    cd "${WORKDIR}/08_krakentools" || exit 1
    git clone https://github.com/jenniferlu717/KrakenTools.git
else
    echo " KrakenTools deja present dans : ${KRAKENTOOLS_DIR}"
fi

cd "${WORKDIR}" || exit 1

# ==============================================================================
# FONCTION : extraire les reads NON classes (taxid 0 / colonne 1 = "U")
# d'un fichier .kraken, a partir d'un fichier fastq source, vers un nouveau
# fastq (utilise ensuite comme entree pour le niveau de k-mer suivant).
# ==============================================================================

extract_unassigned() {
    local KRAKEN_FILE="$1"   # fichier .kraken du niveau precedent
    local SEQ_FILE="$2"      # fastq source du niveau precedent
    local OUT_FASTQ="$3"     # fastq de sortie (reads non classes uniquement)

    python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
        -k "${KRAKEN_FILE}" \
        -s "${SEQ_FILE}" \
        -t 0 \
        -o "${OUT_FASTQ}" \
        --include-children \
        --fastq-output
}

# ==============================================================================
# STEP 1 : Cascade Kraken2 k35 -> k29 -> k25 sur les reads merged, par
#          echantillon "grouped"
# ==============================================================================

echo "======================================================================="
echo " STEP 1: Cascade Kraken2 (k35 -> k29 -> k25) sur les reads merged"
echo "         Echantillons grouped : ${SAMPLES[*]}"
echo "======================================================================="

for SAMPLE in "${SAMPLES[@]}"; do

    echo ""
    echo "-----------------------------------------------------------------"
    echo " Echantillon : ${SAMPLE} (grouped, merged)"
    echo "-----------------------------------------------------------------"

    MERGED_FASTQ="${ENTREE}/clean_${SAMPLE}_grouped_dedup_fastp_merged.fastq.gz"

    if [[ ! -s "${MERGED_FASTQ}" ]]; then
        echo " ATTENTION: fichier merged introuvable pour ${SAMPLE} (${MERGED_FASTQ}), on saute."
        continue
    fi

    # --------------------------------------------------------------------
    # NIVEAU 1 : k35 sur le fastq merged complet
    # --------------------------------------------------------------------

    K35_KRAKEN="${SORTIE_K35}/${SAMPLE}_merged.kraken"
    K35_REPORT="${SORTIE_K35}/${SAMPLE}_merged.report"

    echo " [k35] Classification Kraken2 sur ${SAMPLE} (merged, complet)"
    kraken2 --conf 0.2 --db "${KRAKEN2_DB_K35}" --threads "${THREADS}" \
        --output "${K35_KRAKEN}" --report "${K35_REPORT}" "${MERGED_FASTQ}"

    N_UNCLASS_K35=$(awk -F'\t' '$1=="U"' "${K35_KRAKEN}" | wc -l)
    echo " [k35] Reads non assignes apres k35 : ${N_UNCLASS_K35}"

    # --------------------------------------------------------------------
    # NIVEAU 2 : k29, uniquement sur les reads non assignes par k35
    # --------------------------------------------------------------------

    UNASSIGNED_K35_FASTQ="${SORTIE_TMP}/${SAMPLE}_unassigned_after_k35.fastq"
    K29_KRAKEN="${SORTIE_K29}/${SAMPLE}_merged_unassigned_k35.kraken"
    K29_REPORT="${SORTIE_K29}/${SAMPLE}_merged_unassigned_k35.report"

    if [[ "${N_UNCLASS_K35}" -gt 0 ]]; then
        echo " [k29] Extraction des reads non assignes par k35"
        extract_unassigned "${K35_KRAKEN}" "${MERGED_FASTQ}" "${UNASSIGNED_K35_FASTQ}"

        if [[ -s "${UNASSIGNED_K35_FASTQ}" ]]; then
            echo " [k29] Classification Kraken2 sur ${SAMPLE} (reads non assignes par k35 uniquement)"
            kraken2 --conf 0.2 --db "${KRAKEN2_DB_K29}" --threads "${THREADS}" \
                --output "${K29_KRAKEN}" --report "${K29_REPORT}" "${UNASSIGNED_K35_FASTQ}"
        else
            echo " [k29] Aucun read extrait (fichier vide), etape ignoree."
            : > "${K29_KRAKEN}"
            : > "${K29_REPORT}"
        fi
    else
        echo " [k29] Aucun read non assigne par k35, etape ignoree."
        : > "${K29_KRAKEN}"
        : > "${K29_REPORT}"
    fi

    N_UNCLASS_K29=0
    if [[ -s "${K29_KRAKEN}" ]]; then
        N_UNCLASS_K29=$(awk -F'\t' '$1=="U"' "${K29_KRAKEN}" | wc -l)
    fi
    echo " [k29] Reads non assignes apres k29 : ${N_UNCLASS_K29}"

    # --------------------------------------------------------------------
    # NIVEAU 3 : k25, uniquement sur les reads non assignes par k29
    # --------------------------------------------------------------------

    UNASSIGNED_K29_FASTQ="${SORTIE_TMP}/${SAMPLE}_unassigned_after_k29.fastq"
    K25_KRAKEN="${SORTIE_K25}/${SAMPLE}_merged_unassigned_k29.kraken"
    K25_REPORT="${SORTIE_K25}/${SAMPLE}_merged_unassigned_k29.report"

    if [[ "${N_UNCLASS_K29}" -gt 0 ]]; then
        echo " [k25] Extraction des reads non assignes par k29"
        extract_unassigned "${K29_KRAKEN}" "${UNASSIGNED_K35_FASTQ}" "${UNASSIGNED_K29_FASTQ}"

        if [[ -s "${UNASSIGNED_K29_FASTQ}" ]]; then
            echo " [k25] Classification Kraken2 sur ${SAMPLE} (reads non assignes par k35 ET k29 uniquement)"
            kraken2 --conf 0.2 --db "${KRAKEN2_DB_K25}" --threads "${THREADS}" \
                --output "${K25_KRAKEN}" --report "${K25_REPORT}" "${UNASSIGNED_K29_FASTQ}"
        else
            echo " [k25] Aucun read extrait (fichier vide), etape ignoree."
            : > "${K25_KRAKEN}"
            : > "${K25_REPORT}"
        fi
    else
        echo " [k25] Aucun read non assigne par k29, etape ignoree."
        : > "${K25_KRAKEN}"
        : > "${K25_REPORT}"
    fi

    N_UNCLASS_K25=0
    if [[ -s "${K25_KRAKEN}" ]]; then
        N_UNCLASS_K25=$(awk -F'\t' '$1=="U"' "${K25_KRAKEN}" | wc -l)
    fi
    echo " [k25] Reads non assignes apres k25 (jamais assignes) : ${N_UNCLASS_K25}"

    # --------------------------------------------------------------------
    # Fichier .kraken combine (k35 + k29 + k25) pour mapDamage ulterieur
    # Chaque read merged de l'echantillon apparait UNE SEULE FOIS, avec
    # l'assignation du niveau qui l'a effectivement classe (ou "U" si
    # jamais classe, meme apres k25).
    # --------------------------------------------------------------------

    COMBINED_KRAKEN="${SORTIE_COMBINED_KRAKEN}/${SAMPLE}_following.kraken"
    cat "${K35_KRAKEN}" "${K29_KRAKEN}" "${K25_KRAKEN}" > "${COMBINED_KRAKEN}" 2>/dev/null || true

    echo " [combine] Fichier kraken combine (k35+k29+k25) : ${COMBINED_KRAKEN}"
    echo "           -> a utiliser avec le fastq merged original (${MERGED_FASTQ})"
    echo "              dans extract_kraken_reads.py pour recuperer, pour un"
    echo "              taxid donne, TOUS les reads assignes quel que soit le"
    echo "              niveau k35/k29/k25 (cf. logique de 52_mapdamage.sh)."

    # Nettoyage des fastq temporaires intermediaires
    rm -f "${UNASSIGNED_K35_FASTQ}" "${UNASSIGNED_K29_FASTQ}"

    echo " >>> Termine : ${SAMPLE} (cascade k35 -> k29 -> k25) <<<"

done

echo ""
echo "Cascade Kraken2 (k35 -> k29 -> k25) terminee pour tous les echantillons grouped."

# ==============================================================================
# STEP 2 : Tables MPA par niveau de k-mer (KrakenTools kreport2mpa.py)
# ==============================================================================

echo "======================================================================="
echo " STEP 2: Generation des tables MPA par niveau (k35, k29, k25)"
echo "======================================================================="

convert_reports_to_mpa() {
    local IN_DIR="$1"    # repertoire des .report
    local OUT_DIR="$2"   # repertoire de sortie des .mpa
    local LABEL="$3"     # k35 / k29 / k25

    echo ""
    echo " -> Niveau ${LABEL} : ${IN_DIR} -> ${OUT_DIR}"

    for report in "${IN_DIR}"/*.report; do
        [[ -s "${report}" ]] || continue
        base=$(basename "${report}" .report)
        mpa_file="${OUT_DIR}/${base}.mpa"
        echo "    Conversion MPA : ${base}"
        python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" -r "${report}" -o "${mpa_file}"
    done
}

convert_reports_to_mpa "${SORTIE_K35}" "${SORTIE_MPA_K35}" "k35"
convert_reports_to_mpa "${SORTIE_K29}" "${SORTIE_MPA_K29}" "k29"
convert_reports_to_mpa "${SORTIE_K25}" "${SORTIE_MPA_K25}" "k25"

# ==============================================================================
# STEP 3 : Table MPA "finale" par echantillon = somme des comptages
#          k35 + k29 (sur restes) + k25 (sur restes des restes)
# ==============================================================================

echo "======================================================================="
echo " STEP 3: Fusion (somme) des comptages k35+k29+k25 par echantillon"
echo "======================================================================="

SORTIE_MPA_K35_ENV="${SORTIE_MPA_K35}" SORTIE_MPA_K29_ENV="${SORTIE_MPA_K29}" \
SORTIE_MPA_K25_ENV="${SORTIE_MPA_K25}" SORTIE_MPA_PER_SAMPLE_ENV="${SORTIE_MPA_PER_SAMPLE}" \
python3 << 'PYEOF'
import glob
import os

samples = ["cop408", "cop410", "cop412", "cop414", "cop417"]

mpa_dirs = {
    "k35": os.environ["SORTIE_MPA_K35_ENV"],
    "k29": os.environ["SORTIE_MPA_K29_ENV"],
    "k25": os.environ["SORTIE_MPA_K25_ENV"],
}
out_dir = os.environ["SORTIE_MPA_PER_SAMPLE_ENV"]

def load_mpa(path):
    counts = {}
    if not os.path.isfile(path):
        return counts
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            taxon, value = parts[0], parts[1]
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    continue
            counts[taxon] = counts.get(taxon, 0) + value
    return counts

for sample in samples:
    combined = {}
    found_any = False
    for level, d in mpa_dirs.items():
        # noms possibles selon le niveau
        candidates = [
            os.path.join(d, f"{sample}_merged.mpa"),
            os.path.join(d, f"{sample}_merged_unassigned_k35.mpa"),
            os.path.join(d, f"{sample}_merged_unassigned_k29.mpa"),
        ]
        for c in candidates:
            if os.path.isfile(c):
                found_any = True
                cts = load_mpa(c)
                for taxon, val in cts.items():
                    combined[taxon] = combined.get(taxon, 0) + val

    if not found_any:
        print(f"  [WARN] Aucune table mpa trouvee pour {sample}, on saute.")
        continue

    out_path = os.path.join(out_dir, f"{sample}_following.mpa")
    with open(out_path, "w") as out:
        for taxon in sorted(combined.keys()):
            out.write(f"{taxon}\t{combined[taxon]}\n")
    print(f"  -> Table mpa finale (k35+k29+k25) ecrite : {out_path}")
PYEOF

echo ""
echo "Fusion par echantillon terminee."

# ==============================================================================
# STEP 4 : Table finale groupant TOUS les echantillons grouped
#          (combine_mpa.py sur les tables per_sample_combined)
# ==============================================================================

echo "======================================================================="
echo " STEP 4: Table finale combinee — tous echantillons grouped (k35+k29+k25)"
echo "======================================================================="

FINAL_MPA_FILES=()
for SAMPLE in "${SAMPLES[@]}"; do
    F="${SORTIE_MPA_PER_SAMPLE}/${SAMPLE}_following.mpa"
    [[ -s "${F}" ]] && FINAL_MPA_FILES+=("${F}")
done

if [[ ${#FINAL_MPA_FILES[@]} -gt 0 ]]; then
    python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" \
        -i "${FINAL_MPA_FILES[@]}" \
        -o "${SORTIE_MPA}/combined_all_samples_following.tsv"
    echo " -> Table finale : ${SORTIE_MPA}/combined_all_samples_following.tsv"
else
    echo " AUCUNE table mpa finale disponible, combine_mpa non execute."
fi

# On garde aussi les tables combinees par niveau (comme dans 07_mpa_tables-2.sh)
for LEVEL_DIR in "${SORTIE_MPA_K35}:k35" "${SORTIE_MPA_K29}:k29" "${SORTIE_MPA_K25}:k25"; do
    IFS=":" read -r DIR LABEL <<< "${LEVEL_DIR}"
    MPA_FILES=("${DIR}"/*.mpa)
    if [[ -e "${MPA_FILES[0]}" ]]; then
        python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" -i "${MPA_FILES[@]}" -o "${DIR}/combined_all.tsv"
        echo " -> Table combinee ${LABEL} : ${DIR}/combined_all.tsv"
    fi
done

# ==============================================================================
# RESUME FINAL
# ==============================================================================

echo ""
echo "======================================================================="
echo " ALL STEPS COMPLETED — Kraken following (k35 -> k29 -> k25)"
echo "======================================================================="
echo " Kraken2 k35 (initial)                 : ${SORTIE_K35}"
echo " Kraken2 k29 (restes de k35)           : ${SORTIE_K29}"
echo " Kraken2 k25 (restes de k29)           : ${SORTIE_K25}"
echo " MPA par niveau                        : ${SORTIE_MPA}/{k35,k29,k25}/*.mpa"
echo " MPA final par echantillon (k35+29+25) : ${SORTIE_MPA_PER_SAMPLE}/*_following.mpa"
echo " Table finale tous echantillons        : ${SORTIE_MPA}/combined_all_samples_following.tsv"
echo " Fichiers .kraken combines (mapDamage) : ${SORTIE_COMBINED_KRAKEN}/*_following.kraken"
echo "======================================================================="
echo ""
echo "Pour mapDamage : utiliser <sample>_following.kraken (combined_kraken_for_mapdamage)"
echo "avec le fastq merged original clean_<sample>_grouped_dedup_fastp_merged.fastq.gz"
echo "dans extract_kraken_reads.py -t <TaxID> pour recuperer TOUS les reads du taxon,"
echo "quel que soit le niveau (k35/k29/k25) qui les a assignes, en une seule extraction."
