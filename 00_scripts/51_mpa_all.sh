#!/bin/bash

#SBATCH --job-name=51_mpa_all
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=36
#SBATCH --mem=100G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/51_mpa_all.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/51_mpa_all.out"

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics


################################################################################
# 51_mpa_all.sh
#
# Objectifs :
#   1. Convertir les rapports Kraken2 et Bracken en fichiers .mpa.
#   2. Inclure k25, k29, k35 ET pracken.
#   3. Combiner les .mpa dans une table finale unique.
#   4. Construire des noms de colonnes standardises :
#
#      k25_merged_cop408_run1_bracken
#      k25_merged_cop408_run1_kraken
#      pracken_merged_cop408_run1_kraken
#
# Important :
#   - pracken est traite comme une base Kraken supplementaire.
#   - Des colonnes pracken_*_bracken ne seront ajoutees que si des fichiers
#     *.bracken.report existent effectivement dans 06_bracken/pracken/.
################################################################################

set -Eeuo pipefail

trap 'echo "[ERREUR] Ligne ${LINENO} : ${BASH_COMMAND}" >&2' ERR

################################################################################
# 0. ENVIRONNEMENT
################################################################################

module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics

################################################################################
# 1. PARAMETRES ET CHEMINS
################################################################################

WORKDIR="/home/plstenge/coprolites_comparison"
OUT_ROOT="${WORKDIR}/11_per_run_analysis"

KRAKENTOOLS_DIR="${WORKDIR}/08_krakentools/KrakenTools"

KRAKEN_K25_DIR="${OUT_ROOT}/05_kraken2_k25"
KRAKEN_K29_DIR="${OUT_ROOT}/05_kraken2_k29"
KRAKEN_K35_DIR="${OUT_ROOT}/05_kraken2_k35"
KRAKEN_PRACKEN_DIR="${OUT_ROOT}/05_kraken2_pracken"

BRACKEN_DIR="${OUT_ROOT}/06_bracken"

MPA_DIR="${OUT_ROOT}/07_mpa_tables"

FINAL_MPA="${MPA_DIR}/combined_all_k25_k29_k35_pracken_merged_unmerged_bracken_kraken.tsv"

RUN_LABELS=(
    "run1"
    "run2"
    "run3"
    "run4"
    "run5"
)

SAMPLES=(
    "cop408"
    "cop410"
    "cop412"
    "cop414"
    "cop417"
)

KLABELS=(
    "k25"
    "k29"
    "k35"
    "pracken"
)

MODES=(
    "merged"
    "unmerged"
)

################################################################################
# 2. FONCTIONS
################################################################################

log() {
    echo "[$(date '+%F %T')] $*"
}

die() {
    echo "[ERREUR] $*" >&2
    exit 1
}

check_file() {
    [[ -f "$1" ]] || die "Fichier introuvable : $1"
}

check_dir() {
    [[ -d "$1" ]] || die "Repertoire introuvable : $1"
}

check_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Commande introuvable : $1"
}

################################################################################
# 3. VERIFICATIONS INITIALES
################################################################################

log "======================================================================="
log " VERIFICATIONS INITIALES"
log "======================================================================="

check_cmd python3

check_file "${KRAKENTOOLS_DIR}/kreport2mpa.py"
check_file "${KRAKENTOOLS_DIR}/combine_mpa.py"

check_dir "${KRAKEN_K25_DIR}"
check_dir "${KRAKEN_K29_DIR}"
check_dir "${KRAKEN_K35_DIR}"
check_dir "${KRAKEN_PRACKEN_DIR}"

log "KrakenTools     : ${KRAKENTOOLS_DIR}"
log "MPA_DIR         : ${MPA_DIR}"
log "Table finale    : ${FINAL_MPA}"
log "Kraken pracken  : ${KRAKEN_PRACKEN_DIR}"

################################################################################
# 4. CREATION DES REPERTOIRES DE SORTIE MPA
################################################################################

log "======================================================================="
log " CREATION DES REPERTOIRES MPA"
log "======================================================================="

for KLABEL in "${KLABELS[@]}"; do
    for MODE in "${MODES[@]}"; do
        mkdir -p "${MPA_DIR}/${KLABEL}/${MODE}"
    done
done

################################################################################
# 5. CONVERSION RAPPORTS -> MPA
################################################################################

# Usage :
# make_mpa_files "kraken"  "k25"     "/.../05_kraken2_k25"
# make_mpa_files "bracken" "k25"     "/.../06_bracken"
#
# Convention Kraken :
# <REPORT_ROOT>/<run>/<sample>/<sample_run>_<merged|unmerged>.report
#
# Convention Bracken :
# <REPORT_ROOT>/<klabel>/<run>/<sample>/<sample_run>_<merged|unmerged>.bracken.report

make_mpa_files() {

    local SOURCE="$1"
    local KLABEL="$2"
    local REPORT_ROOT="$3"

    local N_CREATED=0
    local N_MISSING=0

    log "======================================================================="
    log " CONVERSION MPA : ${SOURCE} / ${KLABEL}"
    log "======================================================================="

    for RUN in "${RUN_LABELS[@]}"; do
        for SAMPLE in "${SAMPLES[@]}"; do

            UNIT="${SAMPLE}_${RUN}"

            for MODE in "${MODES[@]}"; do

                if [[ "${SOURCE}" == "kraken" ]]; then
                    REPORT_FILE="${REPORT_ROOT}/${RUN}/${SAMPLE}/${UNIT}_${MODE}.report"
                elif [[ "${SOURCE}" == "bracken" ]]; then
                    REPORT_FILE="${REPORT_ROOT}/${KLABEL}/${RUN}/${SAMPLE}/${UNIT}_${MODE}.bracken.report"
                else
                    die "Source inconnue : ${SOURCE}"
                fi

                MPA_FILE="${MPA_DIR}/${KLABEL}/${MODE}/${UNIT}_${MODE}.${SOURCE}.mpa"

                if [[ ! -s "${REPORT_FILE}" ]]; then
                    ((N_MISSING+=1))
                    continue
                fi

                log "[MPA ${SOURCE}/${KLABEL}] ${UNIT} (${MODE})"

                python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" \
                    -r "${REPORT_FILE}" \
                    -o "${MPA_FILE}"

                [[ -s "${MPA_FILE}" ]] || die \
                    "Le fichier MPA est vide ou absent apres conversion : ${MPA_FILE}"

                ((N_CREATED+=1))

            done
        done
    done

    log "[BILAN ${SOURCE}/${KLABEL}]"
    log "  MPA generes       : ${N_CREATED}"
    log "  Rapports absents  : ${N_MISSING}"
}

################################################################################
# 6. GENERATION DES MPA KRAKEN
################################################################################

log "======================================================================="
log " GENERATION DES MPA KRAKEN"
log "======================================================================="

make_mpa_files "kraken" "k25" "${KRAKEN_K25_DIR}"
make_mpa_files "kraken" "k29" "${KRAKEN_K29_DIR}"
make_mpa_files "kraken" "k35" "${KRAKEN_K35_DIR}"
make_mpa_files "kraken" "pracken" "${KRAKEN_PRACKEN_DIR}"

################################################################################
# 7. GENERATION DES MPA BRACKEN
################################################################################

log "======================================================================="
log " GENERATION DES MPA BRACKEN"
log "======================================================================="

if [[ -d "${BRACKEN_DIR}" ]]; then

    make_mpa_files "bracken" "k25" "${BRACKEN_DIR}"
    make_mpa_files "bracken" "k29" "${BRACKEN_DIR}"
    make_mpa_files "bracken" "k35" "${BRACKEN_DIR}"

    # Bracken pracken est optionnel : il n'est lance que si les rapports existent.
    N_PRACKEN_BRACKEN_REPORTS=$(
        find "${BRACKEN_DIR}/pracken" \
            -type f \
            -name "*.bracken.report" \
            2>/dev/null \
        | wc -l
    )

    if [[ "${N_PRACKEN_BRACKEN_REPORTS}" -gt 0 ]]; then
        log "Rapports Bracken pracken detectes : ${N_PRACKEN_BRACKEN_REPORTS}"
        make_mpa_files "bracken" "pracken" "${BRACKEN_DIR}"
    else
        log "[INFO] Aucun rapport Bracken pracken detecte."
        log "[INFO] Seules les colonnes pracken_*_kraken seront ajoutees."
    fi

else
    log "[ATTENTION] Repertoire Bracken absent : ${BRACKEN_DIR}"
    log "[ATTENTION] Aucune colonne Bracken ne sera ajoutee."
fi

################################################################################
# 8. CONTROLE SPECIFIQUE PRACKEN
################################################################################

log "======================================================================="
log " CONTROLE DES FICHIERS MPA PRACKEN"
log "======================================================================="

N_PRACKEN_MERGED=$(
    find "${MPA_DIR}/pracken/merged" \
        -maxdepth 1 \
        -type f \
        -name "*.kraken.mpa" \
    | wc -l
)

N_PRACKEN_UNMERGED=$(
    find "${MPA_DIR}/pracken/unmerged" \
        -maxdepth 1 \
        -type f \
        -name "*.kraken.mpa" \
    | wc -l
)

log "MPA Kraken pracken merged   : ${N_PRACKEN_MERGED}"
log "MPA Kraken pracken unmerged : ${N_PRACKEN_UNMERGED}"

if [[ "${N_PRACKEN_MERGED}" -eq 0 && "${N_PRACKEN_UNMERGED}" -eq 0 ]]; then
    die "Aucun .kraken.mpa pracken produit. Verifier ${KRAKEN_PRACKEN_DIR}"
fi

log "Liste des MPA pracken :"
find "${MPA_DIR}/pracken" \
    -type f \
    -name "*.mpa" \
    | sort

################################################################################
# 9. CONSTRUCTION DE LA LISTE ORDONNEE DES MPA
################################################################################

log "======================================================================="
log " CONSTRUCTION DE LA LISTE DES FICHIERS MPA"
log "======================================================================="

declare -a MPA_FILES=()
declare -a HEADER_COLS=()

# Ordre final :
# k25 -> k29 -> k35 -> pracken
# merged -> unmerged
# bracken -> kraken
# cop408 -> cop417
# run1 -> run5
for KLABEL in "${KLABELS[@]}"; do
    for MODE in "${MODES[@]}"; do
        for SOURCE in "bracken" "kraken"; do

            MPA_SUBDIR="${MPA_DIR}/${KLABEL}/${MODE}"

            while IFS= read -r MPA_FILE; do

                MPA_FILES+=("${MPA_FILE}")

                BASENAME_FILE=$(basename "${MPA_FILE}")

                # cop408_run1_merged.kraken.mpa -> cop408
                SAMPLE=$(echo "${BASENAME_FILE}" | sed -E 's/^([^_]+)_.*$/\1/')

                # cop408_run1_merged.kraken.mpa -> run1
                RUN=$(echo "${BASENAME_FILE}" | sed -E 's/^[^_]+_([^_]+)_.*$/\1/')

                # cop408_run1_merged.kraken.mpa -> kraken
                PARSED_SOURCE=$(
                    echo "${BASENAME_FILE}" \
                    | sed -E 's/^.*\.(bracken|kraken)\.mpa$/\1/'
                )

                # Exemple final :
                # pracken_merged_cop408_run1_kraken
                COLNAME="${KLABEL}_${MODE}_${SAMPLE}_${RUN}_${PARSED_SOURCE}"

                HEADER_COLS+=("${COLNAME}")

            done < <(
                find "${MPA_SUBDIR}" \
                    -maxdepth 1 \
                    -type f \
                    -name "*.${SOURCE}.mpa" \
                    | sort
            )

        done
    done
done

[[ "${#MPA_FILES[@]}" -gt 0 ]] || die \
    "Aucun fichier MPA n'a ete detecte dans ${MPA_DIR}"

[[ "${#MPA_FILES[@]}" -eq "${#HEADER_COLS[@]}" ]] || die \
    "Incoherence entre fichiers MPA et noms de colonnes"

log "Nombre total de fichiers MPA : ${#MPA_FILES[@]}"

log "Colonnes pracken prevues :"
printf '%s\n' "${HEADER_COLS[@]}" | grep '^pracken_' || true

################################################################################
# 10. COMBINAISON DES MPA
################################################################################

log "======================================================================="
log " COMBINAISON DES FICHIERS MPA"
log "======================================================================="

rm -f "${FINAL_MPA}" "${FINAL_MPA}.tmp"

python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" \
    -i "${MPA_FILES[@]}" \
    -o "${FINAL_MPA}"

[[ -s "${FINAL_MPA}" ]] || die \
    "La table finale est absente ou vide : ${FINAL_MPA}"

################################################################################
# 11. REECRITURE DU HEADER
################################################################################

log "======================================================================="
log " REECRITURE DU HEADER"
log "======================================================================="

HEADER_LINE=$(head -n 1 "${FINAL_MPA}")

IFS=$'\t' read -r -a OLD_HEADER <<< "${HEADER_LINE}"

EXPECTED_NCOL=$(( ${#MPA_FILES[@]} + 1 ))

if [[ "${#OLD_HEADER[@]}" -ne "${EXPECTED_NCOL}" ]]; then
    die "Nombre de colonnes invalide apres combine_mpa.py : \
${#OLD_HEADER[@]} detectees ; ${EXPECTED_NCOL} attendues"
fi

{
    # Premiere colonne : Classification
    printf '%s' "${OLD_HEADER[0]}"

    # Colonnes correspondant exactement a l'ordre de MPA_FILES
    for COLNAME in "${HEADER_COLS[@]}"; do
        printf '\t%s' "${COLNAME}"
    done

    printf '\n'

    # Toutes les lignes sauf le header d'origine
    tail -n +2 "${FINAL_MPA}"

} > "${FINAL_MPA}.tmp"

mv "${FINAL_MPA}.tmp" "${FINAL_MPA}"

################################################################################
# 12. CONTROLES FINAUX
################################################################################

log "======================================================================="
log " CONTROLES FINAUX"
log "======================================================================="

FINAL_NCOL=$(head -n 1 "${FINAL_MPA}" | awk -F'\t' '{print NF}')

N_PRACKEN_COLS=$(
    head -n 1 "${FINAL_MPA}" \
    | tr '\t' '\n' \
    | grep -c '^pracken_' || true
)

log "Table finale : ${FINAL_MPA}"
log "Nombre total de colonnes : ${FINAL_NCOL}"
log "Nombre de colonnes pracken : ${N_PRACKEN_COLS}"

if [[ "${N_PRACKEN_COLS}" -eq 0 ]]; then
    die "Aucune colonne pracken dans le fichier final"
fi

log "Colonnes pracken effectivement presentes :"

head -n 1 "${FINAL_MPA}" \
    | tr '\t' '\n' \
    | grep '^pracken_'

log "======================================================================="
log " TERMINE AVEC SUCCES"
log "======================================================================="
