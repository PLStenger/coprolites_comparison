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

#!/bin/bash

#!/bin/bash

###############################################################################
# Objectif :
#  1) Générer les fichiers MPA à partir des reports Bracken et Kraken
#  2) Combiner tous les MPA dans une table unique avec header informatif
#     Ordre colonnes :
#       k25/k29/k35 -> merged/unmerged -> cop4* -> run* -> bracken/kraken
###############################################################################

set -euo pipefail

########################
# PARAMÈTRES GLOBAUX  #
########################

WORKDIR="/home/plstenge/coprolites_comparison"
OUT_ROOT="${WORKDIR}/11_per_run_analysis"

KRAKEN_K25_DIR="${OUT_ROOT}/05_kraken2_k25"
KRAKEN_K29_DIR="${OUT_ROOT}/05_kraken2_k29"
KRAKEN_K35_DIR="${OUT_ROOT}/05_kraken2_k35"

BRACKEN_DIR="${OUT_ROOT}/06_bracken"

MPA_DIR="${OUT_ROOT}/07_mpa_tables"
KRAKENTOOLS_DIR="${WORKDIR}/08_krakentools/KrakenTools"

FINAL_MPA="${MPA_DIR}/combined_all_k25_k29_k35_merged_unmerged_bracken_kraken.tsv"

# Même jeux que dans ton pipeline
RUN_LABELS=("run1" "run2" "run3" "run4" "run5")
SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")

######################################
# 1. Générer les MPA Bracken/Kraken #
######################################

echo "== Génération des fichiers MPA Bracken + Kraken =="

mkdir -p "${MPA_DIR}/k25/merged" "${MPA_DIR}/k25/unmerged"
mkdir -p "${MPA_DIR}/k29/merged" "${MPA_DIR}/k29/unmerged"
mkdir -p "${MPA_DIR}/k35/merged" "${MPA_DIR}/k35/unmerged"

gen_mpa_for_source() {
  local SOURCE="$1"        # "bracken" ou "kraken"
  local KLABEL="$2"        # k25/k29/k35
  local REPORT_ROOT="$3"   # racine des reports (06_bracken/k*/ ou 05_kraken2_k*/)
  local MPA_OUT_ROOT="$4"  # racine des .mpa (07_mpa_tables/k*/)

  for run in "${RUN_LABELS[@]}"; do
    for sample in "${SAMPLES[@]}"; do
      UNIT="${sample}_${run}"

      # merged
      if [[ "${SOURCE}" == "bracken" ]]; then
        REPORT_MERGED="${REPORT_ROOT}/${KLABEL}/${run}/${sample}/${UNIT}_merged.bracken.report"
      else
        REPORT_MERGED="${REPORT_ROOT}/${run}/${sample}/${UNIT}_merged.report"
      fi

      if [[ -f "${REPORT_MERGED}" ]]; then
        MPA_MERGED="${MPA_OUT_ROOT}/${KLABEL}/merged/${UNIT}_merged.${SOURCE}.mpa"
        echo "  [${SOURCE} ${KLABEL}] ${UNIT} (merged) -> ${MPA_MERGED}"
        python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" -r "${REPORT_MERGED}" -o "${MPA_MERGED}"
      fi

      # unmerged
      if [[ "${SOURCE}" == "bracken" ]]; then
        REPORT_UNMERGED="${REPORT_ROOT}/${KLABEL}/${run}/${sample}/${UNIT}_unmerged.bracken.report"
      else
        REPORT_UNMERGED="${REPORT_ROOT}/${run}/${sample}/${UNIT}_unmerged.report"
      fi

      if [[ -f "${REPORT_UNMERGED}" ]]; then
        MPA_UNMERGED="${MPA_OUT_ROOT}/${KLABEL}/unmerged/${UNIT}_unmerged.${SOURCE}.mpa"
        echo "  [${SOURCE} ${KLABEL}] ${UNIT} (unmerged) -> ${MPA_UNMERGED}"
        python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" -r "${REPORT_UNMERGED}" -o "${MPA_UNMERGED}"
      fi

    done
  done
}

# Bracken : reports dans 06_bracken/k*/run*/cop*/
gen_mpa_for_source "bracken" "k25" "${BRACKEN_DIR}" "${MPA_DIR}"
gen_mpa_for_source "bracken" "k29" "${BRACKEN_DIR}" "${MPA_DIR}"
gen_mpa_for_source "bracken" "k35" "${BRACKEN_DIR}" "${MPA_DIR}"

# Kraken : reports dans 05_kraken2_k*/run*/cop*/
gen_mpa_for_source "kraken" "k25" "${KRAKEN_K25_DIR}" "${MPA_DIR}"
gen_mpa_for_source "kraken" "k29" "${KRAKEN_K29_DIR}" "${MPA_DIR}"
gen_mpa_for_source "kraken" "k35" "${KRAKEN_K35_DIR}" "${MPA_DIR}"

######################################
# 2. Construire la liste des MPA    #
######################################

echo "== Construction de la liste ordonnée des .mpa =="

declare -a MPA_FILES=()

for KLABEL in k25 k29 k35; do
  for MODE in merged unmerged; do
    MPA_SUBDIR="${MPA_DIR}/${KLABEL}/${MODE}"

    # Bracken d'abord
    if [[ -d "${MPA_SUBDIR}" ]]; then
      while IFS= read -r f; do
        MPA_FILES+=("$f")
      done < <(find "${MPA_SUBDIR}" -maxdepth 1 -type f -name "*.bracken.mpa" | sort)
    fi

    # Kraken ensuite
    if [[ -d "${MPA_SUBDIR}" ]]; then
      while IFS= read -r f; do
        MPA_FILES+=("$f")
      done < <(find "${MPA_SUBDIR}" -maxdepth 1 -type f -name "*.kraken.mpa" | sort)
    fi

  done
done

if [[ ${#MPA_FILES[@]} -eq 0 ]]; then
  echo "[ERREUR] Aucun fichier .mpa généré dans ${MPA_DIR}/k*/merged|unmerged (.bracken.mpa / .kraken.mpa)."
  exit 1
fi

echo "[INFO] Nombre de fichiers MPA détectés : ${#MPA_FILES[@]}"
printf '%s\n' "${MPA_FILES[@]}"

########################################
# 3. Construire le header informatif   #
########################################

echo "== Construction du header informatif =="

declare -a HEADER_COLS=()

for f in "${MPA_FILES[@]}"; do
  basename_f=$(basename "$f")  # ex: cop408_run1_merged.bracken.mpa

  # klabel : k25, k29, k35
  klabel=$(echo "$f" | awk -F'/' '{for(i=1;i<=NF;i++) if($i ~ /^k[0-9]+$/){print $i; break}}')

  # mode : merged ou unmerged
  mode=$(echo "$f" | awk -F'/' '{for(i=1;i<=NF;i++) if($i == "merged" || $i == "unmerged"){print $i; break}}')

  # sample : cop408, etc.
  sample=$(echo "$basename_f" | sed -E 's/^([^_]+)_.*$/\1/')

  # run : run1, run2, etc.
  run=$(echo "$basename_f" | sed -E 's/^[^_]+_([^_]+)_.*$/\1/')

  # source : bracken ou kraken
  source=$(echo "$basename_f" | sed -E 's/^.*_([^_]+)\.mpa$/\1/')

  colname="${klabel}_${mode}_${sample}_${run}_${source}"
  HEADER_COLS+=("$colname")
done

echo "[INFO] Nombre de colonnes dans le header : ${#HEADER_COLS[@]}"

#####################################################
# 4. Combiner toutes les MPA avec combine_mpa.py   #
#####################################################

echo "== Combinaison des fichiers MPA avec combine_mpa.py =="

python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" \
  -i "${MPA_FILES[@]}" \
  -o "${FINAL_MPA}"

if [[ ! -f "${FINAL_MPA}" ]]; then
  echo "[ERREUR] Le fichier final ${FINAL_MPA} n'a pas été créé."
  exit 1
fi

#####################################################
# 5. Réécrire le header de la table finale         #
#####################################################

echo "== Réécriture du header de la table finale =="

tmp_out="${FINAL_MPA}.tmp"

{
  IFS=$'\n' read -r header_line
  IFS=$'\t' read -r -a header_fields <<< "$header_line"

  new_header=("${header_fields[0]}")  # colonne de classification

  for col in "${HEADER_COLS[@]}"; do
    new_header+=("$col")
  done

  printf '%s' "${new_header[0]}"
  for ((i=1; i<${#new_header[@]}; i++)); do
    printf '\t%s' "${new_header[i]}"
  done
  printf '\n'

  cat

} < "${FINAL_MPA}" > "${tmp_out}"

mv "${tmp_out}" "${FINAL_MPA}"

echo "[INFO] Table finale écrite dans : ${FINAL_MPA}"
echo "[INFO] Ordre des colonnes : k25/k29/k35 -> merged/unmerged -> cop4* -> run* -> bracken/kraken"
echo "[INFO] Exemple de nom de colonne : k25_merged_cop408_run1_bracken"
