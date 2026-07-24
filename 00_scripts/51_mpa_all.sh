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

# Racine des résultats 
WORKDIR="/home/plstenge/coprolites_comparison"
OUT_ROOT="${WORKDIR}/11_per_run_analysis"
MPA_DIR="${OUT_ROOT}/07_mpa_tables"

# Fichier de sortie final
FINAL_MPA="${MPA_DIR}/combined_all_kmerged_unmerged.tsv"

# 1. Construire la liste ordonnée de tous les fichiers .mpa
# Ordre : k25/k29/k35, puis merged/unmerged, puis cop4*, puis run*
declare -a MPA_FILES=()

for KLABEL in k25 k29 k35; do
  for MODE in merged unmerged; do
    # On force un ordre cop4* puis run* via sort
    while IFS= read -r f; do
      MPA_FILES+=("$f")
    done < <(find "${MPA_DIR}/${KLABEL}/${MODE}" -maxdepth 1 -type f -name "*.mpa" \
             | sort)
  done
done

if [[ ${#MPA_FILES[@]} -eq 0 ]]; then
  echo "Aucun fichier .mpa trouvé, vérifie les chemins (${MPA_DIR}/k*/merged|unmerged)."
  exit 1
fi

echo "Nombre de fichiers MPA détectés : ${#MPA_FILES[@]}"
printf '%s\n' "${MPA_FILES[@]}"

# 2. Construire le header informatif pour chaque fichier MPA
# Format souhaité : k*_*merged_*cop4*_*run*_*bracken
# ajustable en fonction de la provenance (Bracken ou Kraken)
declare -a HEADER_COLS=()

for f in "${MPA_FILES[@]}"; do
  # Exemple de chemin:
  # /.../07_mpa_tables/k25/merged/cop408_run1_merged.mpa
  #                               ^^^^^^^^^^^^^^^^^^^^^^ basename

  basename_f=$(basename "$f")  # cop408_run1_merged.mpa
  klabel=$(echo "$f" | awk -F'/' '{for(i=1;i<=NF;i++) if($i ~ /^k[0-9]+$/){print $i; break}}')
  mode=$(echo "$f" | awk -F'/' '{for(i=1;i<=NF;i++) if($i == "merged" || $i == "unmerged"){print $i; break}}')

  sample=$(echo "$basename_f" | sed -E 's/^([^_]+)_.*$/\1/')
  run=$(echo "$basename_f"    | sed -E 's/^[^_]+_([^_]+)_.*$/\1/')

  # suffixe, à adapter si tu veux "bracken" ou "kraken"
  suffix="bracken"

  # Construire le nom de colonne dans l’ordre que tu souhaites
  colname="${klabel}_${mode}_${sample}_${run}_${suffix}"
  HEADER_COLS+=("$colname")
done

# 3. Utiliser KrakenTools/combine_mpa.py pour combiner toutes les colonnes
# Adapte KRAKENTOOLS_DIR si nécessaire
KRAKENTOOLS_DIR="${WORKDIR}/08_krakentools/KrakenTools"

echo "Combinaison des fichiers MPA avec KrakenTools/combine_mpa.py ..."
python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" \
  -i "${MPA_FILES[@]}" \
  -o "${FINAL_MPA}"

if [[ ! -f "${FINAL_MPA}" ]]; then
  echo "Erreur : le fichier final ${FINAL_MPA} n'a pas été créé."
  exit 1
fi

# 4. Réécrire la ligne de header avec les noms informatifs
# combine_mpa.py met en principe les noms de fichiers .mpa en header.
# On les remplace par HEADER_COLS dans le même ordre.

tmp_out="${FINAL_MPA}.tmp"

{
  # lire la première ligne (header original)
  IFS=$'\n' read -r header_line
  # séparer sur tabulation
  IFS=$'\t' read -r -a header_fields <<< "$header_line"

  # header_fields[0] est en général la colonne de classification (par ex. '#Classification')
  # On laisse tel quel, et on remplace uniquement les colonnes Sample.
  new_header=("${header_fields[0]}")

  # Construire les nouvelles colonnes à partir de HEADER_COLS
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

echo "Table finale écrite dans : ${FINAL_MPA}"
echo "Header reconstruit avec des colonnes de la forme : k*_merged/unmerged_cop4*_run*_bracken"
