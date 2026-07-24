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

###############################################################################
# Script : 99_combine_all_mpa_bracken_kraken.sh
# Objectif :
#  - Construire une table unique MPA combinant Bracken + Kraken
#  - Ordre des colonnes :
#      k25/k29/k35 -> merged/unmerged -> cop4* -> run* -> bracken/kraken
#  - Header informatif :
#      k25_merged_cop408_run1_bracken, k25_merged_cop408_run1_kraken, etc.
###############################################################################

########################
# PARAMÈTRES GLOBAUX  #
########################

WORKDIR="/home/plstenge/coprolites_comparison"
OUT_ROOT="${WORKDIR}/11_per_run_analysis"
MPA_DIR="${OUT_ROOT}/07_mpa_tables"

# Répertoire où est installé KrakenTools (comme dans ton pipeline)
KRAKENTOOLS_DIR="${WORKDIR}/08_krakentools/KrakenTools"

# Fichier de sortie final
FINAL_MPA="${MPA_DIR}/combined_all_k25_k29_k35_merged_unmerged_bracken_kraken.tsv"

######################################
# 1. Construire la liste des MPA    #
######################################

# Ordre voulu :
#  - k25 puis k29 puis k35
#  - merged puis unmerged
#  - cop4* puis run* (on laisse l’ordre donné par le nom de fichier)
#  - pour chaque condition, bracken puis kraken

declare -a MPA_FILES=()

for KLABEL in k25 k29 k35; do
  for MODE in merged unmerged; do
    # On suppose que les fichiers suivent la convention :
    #  - .../07_mpa_tables/k25/merged/cop408_run1_merged.bracken.mpa
    #  - .../07_mpa_tables/k25/merged/cop408_run1_merged.kraken.mpa
    #
    # On récupère d'abord les MPA Bracken, puis les MPA Kraken.
    BRACKEN_MPA_DIR="${MPA_DIR}/${KLABEL}/${MODE}"
    KRAKEN_MPA_DIR="${MPA_DIR}/${KLABEL}/${MODE}"

    # Bracken d'abord
    if [[ -d "${BRACKEN_MPA_DIR}" ]]; then
      while IFS= read -r f; do
        MPA_FILES+=("$f")
      done < <(find "${BRACKEN_MPA_DIR}" -maxdepth 1 -type f -name "*.bracken.mpa" | sort)
    fi

    # Kraken ensuite
    if [[ -d "${KRAKEN_MPA_DIR}" ]]; then
      while IFS= read -r f; do
        MPA_FILES+=("$f")
      done < <(find "${KRAKEN_MPA_DIR}" -maxdepth 1 -type f -name "*.kraken.mpa" | sort)
    fi

  done
done

if [[ ${#MPA_FILES[@]} -eq 0 ]]; then
  echo "[ERREUR] Aucun fichier .mpa trouvé dans ${MPA_DIR}/k*/merged|unmerged avec suffixes .bracken.mpa ou .kraken.mpa"
  exit 1
fi

echo "[INFO] Nombre de fichiers MPA détectés : ${#MPA_FILES[@]}"
printf '%s\n' "${MPA_FILES[@]}"

########################################
# 2. Construire le header informatif   #
########################################

declare -a HEADER_COLS=()

for f in "${MPA_FILES[@]}"; do
  # Exemple de chemin :
  #   /home/.../07_mpa_tables/k25/merged/cop408_run1_merged.bracken.mpa
  # ou
  #   /home/.../07_mpa_tables/k25/unmerged/cop410_run2_unmerged.kraken.mpa

  basename_f=$(basename "$f")  # cop408_run1_merged.bracken.mpa

  # klabel : k25, k29, k35
  klabel=$(echo "$f" | awk -F'/' '{for(i=1;i<=NF;i++) if($i ~ /^k[0-9]+$/){print $i; break}}')

  # mode : merged ou unmerged
  mode=$(echo "$f" | awk -F'/' '{for(i=1;i<=NF;i++) if($i == "merged" || $i == "unmerged"){print $i; break}}')

  # sample : cop408, etc. (avant le premier "_")
  sample=$(echo "$basename_f" | sed -E 's/^([^_]+)_.*$/\1/')

  # run : run1, run2, etc. (entre premier et deuxième "_")
  run=$(echo "$basename_f" | sed -E 's/^[^_]+_([^_]+)_.*$/\1/')

  # source : bracken ou kraken (dans le suffixe .bracken.mpa / .kraken.mpa)
  # On dépouille le suffixe après le dernier "_" et avant l'extension .mpa
  source=$(echo "$basename_f" | sed -E 's/^.*_([^_]+)\.mpa$/\1/')

  # Construire le nom de colonne dans l’ordre souhaité :
  #  k25/k29/k35 -> merged/unmerged -> cop4* -> run* -> bracken/kraken
  colname="${klabel}_${mode}_${sample}_${run}_${source}"

  HEADER_COLS+=("$colname")
done

echo "[INFO] Nombre de colonnes dans le header : ${#HEADER_COLS[@]}"

#####################################################
# 3. Combiner toutes les MPA avec combine_mpa.py   #
#####################################################

echo "[INFO] Combinaison des fichiers MPA avec KrakenTools/combine_mpa.py ..."
python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" \
  -i "${MPA_FILES[@]}" \
  -o "${FINAL_MPA}"

if [[ ! -f "${FINAL_MPA}" ]]; then
  echo "[ERREUR] Le fichier final ${FINAL_MPA} n'a pas été créé."
  exit 1
fi

#####################################################
# 4. Réécrire la première ligne de header          #
#####################################################

tmp_out="${FINAL_MPA}.tmp"

{
  # Lire la première ligne (header original)
  IFS=$'\n' read -r header_line

  # Séparer sur tabulation
  IFS=$'\t' read -r -a header_fields <<< "$header_line"

  # header_fields[0] = colonne de classification, qu'on laisse telle quelle.
  # Les colonnes suivantes sont les colonnes Sample (noms de fichiers .mpa).
  # On les remplace par HEADER_COLS dans le même ordre que MPA_FILES.
  new_header=("${header_fields[0]}")

  for col in "${HEADER_COLS[@]}"; do
    new_header+=("$col")
  done

  # Écrire le nouveau header
  printf '%s' "${new_header[0]}"
  for ((i=1; i<${#new_header[@]}; i++)); do
    printf '\t%s' "${new_header[i]}"
  done
  printf '\n'

  # Puis recopier le reste du fichier tel quel
  cat

} < "${FINAL_MPA}" > "${tmp_out}"

mv "${tmp_out}" "${FINAL_MPA}"

echo "[INFO] Table finale écrite dans : ${FINAL_MPA}"
echo "[INFO] Header reconstruit avec des colonnes :"
echo "       k25/k29/k35 -> merged/unmerged -> cop4* -> run* -> bracken/kraken"
echo "       Exemple : k25_merged_cop408_run1_bracken, k25_merged_cop408_run1_kraken"

###############################################################################
# Fin du script
###############################################################################
