#!/bin/bash
#SBATCH --job-name=coprolites_merged
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=/home/plstenge/coprolites_comparison/00_scripts/coprolites_merged.err
#SBATCH --output=/home/plstenge/coprolites_comparison/00_scripts/coprolites_merged.out

set -euo pipefail

############################################
# CONFIGURATION GLOBALE
############################################

BASEDIR="/home/plstenge/coprolites_comparison"

BBDUK="/home/plstenge/bbmap/bbduk.sh"
CLUMPIFY="/home/plstenge/bbmap/clumpify.sh"
PHIX="/home/plstenge/bbmap/resources/phix174_ill.ref.fa.gz"

KRAKEN2DB="/home/plstenge/k2_corent_2025-06-09"
KRAKENTOOLSDIR="${BASEDIR}/08_kraken2/KrakenTools"

THREADS=36

echo
echo "Configuration :"
echo "  BASEDIR         = ${BASEDIR}"
echo "  KRAKEN2DB       = ${KRAKEN2DB}"
echo "  KRAKENTOOLSDIR  = ${KRAKENTOOLSDIR}"
echo "  THREADS         = ${THREADS}"
echo

############################################
# DÉFINITION DU CONTRÔLE NÉGATIF (NTC)
############################################

# Fichiers NTC tels que fournis
NTC_RAW_DIR="${BASEDIR}/01_raw_data_merged/Lot_NTC"
NTC_R1_SOURCE="${BASEDIR}/01_raw_data_merged/Lot_NTC/NTC_cop_R1.fastq.gz"
NTC_R2_SOURCE="${BASEDIR}/01_raw_data_merged/Lot_NTC/NTC_cop_R2.fastq.gz"
NTC_SAMPLE="NTCcop"

echo
echo "Dossier et noms du contrôle négatif :"
echo "  NTC_RAW_DIR = ${NTC_RAW_DIR}"
echo "  NTC_R1      = ${NTC_R1_SOURCE}"
echo "  NTC_R2      = ${NTC_R2_SOURCE}"
echo

############################################
# TAPE 0 : ORGANISATION DES DONNÉES
# - 5 lots de runs
# - + 1 lot NTC qui sera traité comme les autres
############################################

echo
echo "===== TAPE 0 : Organisation des données - 5 LOTS + NTC ====="

# Liste des lots "classiques"
declare -a SOURCE_LOTS=(
  "Lot1_illu-4_R1_Ps4_150_Default"
  "Lot2_Run1_R2_Ps6_150_no_filter"
  "Lot4_Run2_R2_Ps6_150_no_filter"
  "Lot6_Run3_R3_Ps8_150_no_filter"
  "Lot9_Run4_R3_Ps8_75_no_filter"
)

# Dictionnaires associant lot → chemin source & mode
declare -A LOT_SOURCE
declare -A LOT_MODE

LOT_SOURCE["Lot1_illu-4_R1_Ps4_150_Default"]="${BASEDIR}/01_raw_data/Lot1_Illumina_R1"
LOT_MODE["Lot1_illu-4_R1_Ps4_150_Default"]="flat"

LOT_SOURCE["Lot2_Run1_R2_Ps6_150_no_filter"]="/storage/groups/gdec/shared/paleo/E1531/final_run1_2025-03-20/AV24-1601_E1531_Ps5_Lane1_Ps6_Lane2"
LOT_MODE["Lot2_Run1_R2_Ps6_150_no_filter"]="subdirs"

LOT_SOURCE["Lot4_Run2_R2_Ps6_150_no_filter"]="/storage/groups/gdec/shared/paleo/E1531/final_run2_2025-04-14/AV24-1601_E1531_Ps5_Ps6_14-04-2025"
LOT_MODE["Lot4_Run2_R2_Ps6_150_no_filter"]="subdirs"

LOT_SOURCE["Lot6_Run3_R3_Ps8_150_no_filter"]="/storage/groups/gdec/shared/paleo/E1531/final_run3_2025-10-08/AV24-1601_E1531_Ps7_Ps8"
LOT_MODE["Lot6_Run3_R3_Ps8_150_no_filter"]="subdirs"

LOT_SOURCE["Lot9_Run4_R3_Ps8_75_no_filter"]="/storage/groups/gdec/shared/paleo/E1531/final_run4_2025-11-04/AV24-1601_E1531_Ps7_Ps8_04-11-2025"
LOT_MODE["Lot9_Run4_R3_Ps8_75_no_filter"]="subdirs"

# Ajout d’un lot dédié au NTC dans la même logique
NTC_LOT="Lot_NTC"
SOURCE_LOTS+=("${NTC_LOT}")
LOT_SOURCE["${NTC_LOT}"]="${NTC_RAW_DIR}"
LOT_MODE["${NTC_LOT}"]="flat"

# Création des dossiers principaux
mkdir -p "${BASEDIR}/01_raw_data_merged"
mkdir -p "${BASEDIR}/00_scripts"
mkdir -p "${BASEDIR}/00_intermediate_merged"

shopt -s nullglob

for lot in "${SOURCE_LOTS[@]}"; do
  SRC_DIR="${LOT_SOURCE[$lot]}"
  MODE="${LOT_MODE[$lot]}"
  DEST_DIR="${BASEDIR}/01_raw_data_merged/${lot}"

  echo
  echo "------------------------------------------"
  echo "Organisation du ${lot}"
  echo " Source      : ${SRC_DIR}"
  echo " Destination : ${DEST_DIR}"
  echo " Mode        : ${MODE}"
  echo "------------------------------------------"

  if [[ ! -d "${SRC_DIR}" ]]; then
    echo "ATTENTION : répertoire source introuvable pour ${lot} : ${SRC_DIR}"
    continue
  fi

  mkdir -p "${DEST_DIR}"

  if [[ "${MODE}" == "flat" ]]; then
    echo "Création de liens symboliques pour ${lot} (mode flat)..."
    for fq in "${SRC_DIR}"/*cop*R{1,2}.fastq.gz; do
      base="$(basename "${fq}")"
      if [[ -e "${DEST_DIR}/${base}" ]]; then
        echo "  Lien déjà présent pour ${base}, on saute."
      else
        ln -s "${fq}" "${DEST_DIR}/${base}"
        echo "  Lien créé : ${base}"
      fi
    done

  elif [[ "${MODE}" == "subdirs" ]]; then
    echo "Recherche des sous-dossiers copXXX dans ${SRC_DIR}..."
    for d in "${SRC_DIR}"/{0..9}cop{0..9}{0..9}{0..9}; do
      [[ -d "${d}" ]] || continue
      foldername="$(basename "${d}")"
      sampleid="${foldername}"

      R1_SRC="${d}/${foldername}_R1.fastq.gz"
      R2_SRC="${d}/${foldername}_R2.fastq.gz"

      if [[ -f "${R1_SRC}" && -f "${R2_SRC}" ]]; then
        R1_DEST="${DEST_DIR}/${sampleid}_R1.fastq.gz"
        R2_DEST="${DEST_DIR}/${sampleid}_R2.fastq.gz"

        if [[ -e "${R1_DEST}" && -e "${R2_DEST}" ]]; then
          echo "  Liens déjà présents pour ${sampleid} dans ${lot}, on saute."
        else
          ln -s "${R1_SRC}" "${R1_DEST}"
          ln -s "${R2_SRC}" "${R2_DEST}"
          echo "  ${sampleid} lié : ${foldername}_R1/R2.fastq.gz → ${sampleid}_R1/R2.fastq.gz"
        fi
      else
        echo "  Fichiers R1/R2 manquants pour ${foldername} dans ${SRC_DIR}"
      fi
    done
  fi
done

shopt -u nullglob

echo
echo "Organisation des données terminée."

############################################
# INTÉGRATION DU CONTRÔLE NTC DANS LE PIPELINE
############################################

echo
echo "===== INTÉGRATION DU CONTRÔLE NÉGATIF NTC ====="

# On s'assure que les fichiers NTC existent dans 01_raw_data_merged/Lot_NTC
mkdir -p "${NTC_RAW_DIR}"

if [[ -f "${NTC_R1_SOURCE}" && -f "${NTC_R2_SOURCE}" ]]; then
  echo "Fichiers NTC présents dans ${NTC_RAW_DIR}."
else
  echo "ATTENTION : fichiers NTC R1/R2 introuvables dans ${NTC_RAW_DIR}."
fi

############################################
# ACTIVATION ENVIRONNEMENT CONDA
############################################

echo
echo "===== Activation environnement conda (metagenomics) ====="

module load conda/4.12.0 || true
source ~/.bashrc
conda activate metagenomics

echo "Environnement activé : metagenomics"
echo

# Vérification rapide de Krona (optionnel)
if command -v ktImportTaxonomy &>/dev/null; then
  echo "Krona disponible."
else
  echo "Krona non trouvé mais les analyses continueront."
fi

############################################
# CRÉATION ARBORESCENCE GLOBALE
############################################

echo
echo "===== Création de l’arborescence ====="

mkdir -p "${BASEDIR}/02_quality_check_raw_merged"
mkdir -p "${BASEDIR}/03_bbduk_merged"
mkdir -p "${BASEDIR}/04_fastuniq_merged"
mkdir -p "${BASEDIR}/05_clumpify_merged"
mkdir -p "${BASEDIR}/06_fastp_merged_reads"
mkdir -p "${BASEDIR}/06_fastp_unmerged_reads"
mkdir -p "${BASEDIR}/07_quality_check_clean_merged_reads"
mkdir -p "${BASEDIR}/07_quality_check_clean_unmerged_reads"
mkdir -p "${BASEDIR}/08_kraken2_merged_reads"
mkdir -p "${BASEDIR}/08_kraken2_unmerged_reads"
mkdir -p "${BASEDIR}/09_krona_merged_reads"
mkdir -p "${BASEDIR}/09_krona_unmerged_reads"
mkdir -p "${BASEDIR}/10_mpa_tables_merged_reads"
mkdir -p "${BASEDIR}/10_mpa_tables_unmerged_reads"
mkdir -p "${BASEDIR}/11_summary_tables"
mkdir -p "${BASEDIR}/12_mapdamage_merged_reads"
mkdir -p "${BASEDIR}/12_mapdamage_unmerged_reads"

echo "Arborescence créée."

############################################
# TAPE 1 : NETTOYAGE PAR LOT
# (BBDuk, FastUniq, Clumpify)
# → à activer si tu veux refaire cette partie
############################################

echo
echo "===== TAPE 1 : Nettoyage des reads par lot (placeholder) ====="
echo "Cette section est supposée contenir BBDuk / FastUniq / Clumpify."
echo "Le script évite de refaire si les outputs existent déjà."
echo

for lot in "${SOURCE_LOTS[@]}"; do
  RAW_DIR="${BASEDIR}/01_raw_data_merged/${lot}"
  CLUMP_DIR="${BASEDIR}/05_clumpify_merged/${lot}"
  mkdir -p "${CLUMP_DIR}"

  # Exemple de garde-fou : si un fichier clumpify existe déjà, on considère le lot comme traité
  if ls "${CLUMP_DIR}"/*dedupclumpifyR1.fastq.gz &>/dev/null; then
    echo "  Lot ${lot} déjà clumpifié, on saute la TAPE 1 pour ce lot."
    continue
  fi

  echo "  [TODO] Nettoyage BBDuk/FastUniq/Clumpify pour ${lot} (à adapter/compléter si nécessaire)."
done

############################################
# TAPE 2 : MERGE PAR RUN PUIS CONCATÉNATION
############################################

echo
echo "===== TAPE 2 : Merge R1/R2 par run (fastp) puis concaténation ====="

INTERMEDIATE_DIR="${BASEDIR}/00_intermediate_merged"
mkdir -p "${INTERMEDIATE_DIR}"

# Échantillons finaux (copXXX + NTC)
declare -a ALL_SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417" "${NTC_SAMPLE}")

echo
echo "---- TAPE 2a : Merge R1/R2 pour chaque lot ----"

for lot in "${SOURCE_LOTS[@]}"; do
  echo
  echo "Fastp merge pour ${lot}..."
  INPUT_DIR="${BASEDIR}/05_clumpify_merged/${lot}"
  OUTPUT_DIR="${INTERMEDIATE_DIR}/${lot}"
  mkdir -p "${OUTPUT_DIR}"

  if [[ ! -d "${INPUT_DIR}" ]]; then
    echo "  Pas de données clumpify pour ${lot}, on saute."
    continue
  fi

  shopt -s nullglob
  for R1 in "${INPUT_DIR}"/*cleancop*dedupclumpifyR1.fastq.gz; do
    R2="${R1/R1.fastq.gz/R2.fastq.gz}"
    [[ -f "${R2}" ]] || continue
    base="$(basename "${R1}" "R1.fastq.gz")"

    MERGED_OUT="${OUTPUT_DIR}/${base}merged.fastq.gz"
    UNMERGED_R1_OUT="${OUTPUT_DIR}/${base}unmergedR1.fastq.gz"
    UNMERGED_R2_OUT="${OUTPUT_DIR}/${base}unmergedR2.fastq.gz"
    JSON_OUT="${OUTPUT_DIR}/${base}fastp.json"
    HTML_OUT="${OUTPUT_DIR}/${base}fastp.html"

    # Garde-fou : si merged existe déjà, on ne refait pas fastp
    if [[ -f "${MERGED_OUT}" && -f "${UNMERGED_R1_OUT}" && -f "${UNMERGED_R2_OUT}" ]]; then
      echo "  fastp déjà exécuté pour ${base}, on saute."
      continue
    fi

    echo "  Merge de ${base}..."
    fastp \
      -i "${R1}" -I "${R2}" \
      --merge \
      --merged_out "${MERGED_OUT}" \
      --out1 "${UNMERGED_R1_OUT}" \
      --out2 "${UNMERGED_R2_OUT}" \
      --json "${JSON_OUT}" \
      --html "${HTML_OUT}" \
      --thread 4 \
      --length_required 30 \
      --qualified_quality_phred 20
  done
  shopt -u nullglob
done

echo
echo "---- TAPE 2b : Concaténation des fichiers merged et unmerged entre runs ----"

FINAL_MERGED_DIR="${BASEDIR}/06_fastp_merged_reads"
FINAL_UNMERGED_DIR="${BASEDIR}/06_fastp_unmerged_reads"
mkdir -p "${FINAL_MERGED_DIR}" "${FINAL_UNMERGED_DIR}"

for sample in "${ALL_SAMPLES[@]}"; do
  echo
  echo "Concaténation pour ${sample}"

  MERGED_FINAL="${FINAL_MERGED_DIR}/${sample}_final_merged.fastq.gz"
  UNM_R1_FINAL="${FINAL_UNMERGED_DIR}/${sample}_final_unmerged_R1.fastq.gz"
  UNM_R2_FINAL="${FINAL_UNMERGED_DIR}/${sample}_final_unmerged_R2.fastq.gz"

  # Si les fichiers finaux existent déjà, on ne refait rien
  if [[ -f "${MERGED_FINAL}" && -f "${UNM_R1_FINAL}" && -f "${UNM_R2_FINAL}" ]]; then
    echo "  Fichiers finaux déjà présents pour ${sample}, on saute."
    continue
  fi

  merged_files=()
  unmerged_r1_files=()
  unmerged_r2_files=()

  for lot in "${SOURCE_LOTS[@]}"; do
    lotdir="${INTERMEDIATE_DIR}/${lot}"
    [[ -d "${lotdir}" ]] || continue

    # merged
    for mf in "${lotdir}"/*"${sample}"*cleancop*dedupclumpifymerged.fastq.gz; do
      [[ -f "${mf}" ]] || continue
      merged_files+=("${mf}")
      echo "  Trouvé merged $(basename "${mf}") dans ${lot}"
    done

    # unmerged
    for uf1 in "${lotdir}"/*"${sample}"*cleancop*dedupclumpifyunmergedR1.fastq.gz; do
      [[ -f "${uf1}" ]] || continue
      uf2="${uf1/R1.fastq.gz/R2.fastq.gz}"
      [[ -f "${uf2}" ]] || continue
      unmerged_r1_files+=("${uf1}")
      unmerged_r2_files+=("${uf2}")
      echo "  Trouvé unmerged $(basename "${uf1}") dans ${lot}"
    done
  done

  if (( ${#merged_files[@]} > 0 )); then
    echo "  Concaténation de ${#merged_files[@]} fichiers merged..."
    cat "${merged_files[@]}" > "${MERGED_FINAL}"
    echo "  Créé ${MERGED_FINAL}"
  else
    echo "  Aucun fichier merged trouvé pour ${sample}."
  fi

  if (( ${#unmerged_r1_files[@]} > 0 )); then
    echo "  Concaténation de ${#unmerged_r1_files[@]} paires unmerged..."
    cat "${unmerged_r1_files[@]}" > "${UNM_R1_FINAL}"
    cat "${unmerged_r2_files[@]}" > "${UNM_R2_FINAL}"
    echo "  Créé ${UNM_R1_FINAL} et ${UNM_R2_FINAL}"
  else
    echo "  Aucun fichier unmerged trouvé pour ${sample}."
  fi
done

############################################
# TAPE 3 : QC MERGED & UNMERGED
############################################

echo
echo "===== TAPE 3 : Contrôle qualité MERGED & UNMERGED ====="

QC_MERGED_DIR="${BASEDIR}/07_quality_check_clean_merged_reads"
QC_UNMERGED_DIR="${BASEDIR}/07_quality_check_clean_unmerged_reads"
mkdir -p "${QC_MERGED_DIR}" "${QC_UNMERGED_DIR}"

echo
echo "FastQC sur les fichiers merged..."
for FILE in "${FINAL_MERGED_DIR}"/*.fastq.gz; do
  [[ -f "${FILE}" ]] || continue
  base="$(basename "${FILE}")"
  if [[ -f "${QC_MERGED_DIR}/${base/.fastq.gz/_fastqc.zip}" ]]; then
    echo "  FastQC déjà présent pour ${base}, on saute."
    continue
  fi
  echo "  FastQC sur ${base}..."
  fastqc "${FILE}" -o "${QC_MERGED_DIR}" -t 4
done

if [[ ! -f "${QC_MERGED_DIR}/multiqc_merged_reads.html" ]]; then
  echo "MultiQC merged..."
  cd "${QC_MERGED_DIR}"
  multiqc . -n "multiqc_merged_reads.html" --force
fi

echo
echo "FastQC sur les fichiers unmerged..."
for FILE in "${FINAL_UNMERGED_DIR}"/*.fastq.gz; do
  [[ -f "${FILE}" ]] || continue
  base="$(basename "${FILE}")"
  if [[ -f "${QC_UNMERGED_DIR}/${base/.fastq.gz/_fastqc.zip}" ]]; then
    echo "  FastQC déjà présent pour ${base}, on saute."
    continue
  fi
  echo "  FastQC sur ${base}..."
  fastqc "${FILE}" -o "${QC_UNMERGED_DIR}" -t 4
done

if [[ ! -f "${QC_UNMERGED_DIR}/multiqc_unmerged_reads.html" ]]; then
  echo "MultiQC unmerged..."
  cd "${QC_UNMERGED_DIR}"
  multiqc . -n "multiqc_unmerged_reads.html" --force
fi

############################################
# TAPE 4 : KRAKEN2 MERGED & UNMERGED
############################################

echo
echo "===== TAPE 4a : Kraken2 - MERGED READS ====="

KRAKEN_MERGED_DIR="${BASEDIR}/08_kraken2_merged_reads"
mkdir -p "${KRAKEN_MERGED_DIR}"

for sample in "${ALL_SAMPLES[@]}"; do
  MERGED="${FINAL_MERGED_DIR}/${sample}_final_merged.fastq.gz"
  [[ -f "${MERGED}" ]] || { echo "  Fichier merged absent pour ${sample}, on saute."; continue; }

  OUT_KRAKEN="${KRAKEN_MERGED_DIR}/${sample}_merged.kraken"
  OUT_REPORT="${KRAKEN_MERGED_DIR}/${sample}_merged.report"

  if [[ -f "${OUT_KRAKEN}" && -f "${OUT_REPORT}" ]]; then
    echo "  Kraken2 merged déjà fait pour ${sample}, on saute."
    continue
  fi

  echo "  Kraken2 pour ${sample} merged..."
  kraken2 \
    --confidence 0.2 \
    --db "${KRAKEN2DB}" \
    --threads "${THREADS}" \
    --output "${OUT_KRAKEN}" \
    --report "${OUT_REPORT}" \
    "${MERGED}"
done

echo
echo "===== TAPE 4b : Kraken2 - UNMERGED READS ====="

KRAKEN_UNMERGED_DIR="${BASEDIR}/08_kraken2_unmerged_reads"
mkdir -p "${KRAKEN_UNMERGED_DIR}"

for sample in "${ALL_SAMPLES[@]}"; do
  R1="${FINAL_UNMERGED_DIR}/${sample}_final_unmerged_R1.fastq.gz"
  R2="${FINAL_UNMERGED_DIR}/${sample}_final_unmerged_R2.fastq.gz"
  if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
    echo "  Fichiers unmerged absents pour ${sample}, on saute."
    continue
  fi

  OUT_KRAKEN="${KRAKEN_UNMERGED_DIR}/${sample}_unmerged.kraken"
  OUT_REPORT="${KRAKEN_UNMERGED_DIR}/${sample}_unmerged.report"

  if [[ -f "${OUT_KRAKEN}" && -f "${OUT_REPORT}" ]]; then
    echo "  Kraken2 unmerged déjà fait pour ${sample}, on saute."
    continue
  fi

  echo "  Kraken2 pour ${sample} unmerged..."
  kraken2 \
    --confidence 0.2 \
    --paired \
    --db "${KRAKEN2DB}" \
    --threads "${THREADS}" \
    --output "${OUT_KRAKEN}" \
    --report "${OUT_REPORT}" \
    "${R1}" "${R2}"
done

############################################
# TAPE 5 : KRONA MERGED & UNMERGED
############################################

echo
echo "===== TAPE 5a : Krona - MERGED READS ====="

KRONA_MERGED_DIR="${BASEDIR}/09_krona_merged_reads"
mkdir -p "${KRONA_MERGED_DIR}"
cd "${KRAKEN_MERGED_DIR}"

if [[ ! -f "${KRONA_MERGED_DIR}/all_samples_merged_krona.html" ]]; then
  if ls *_merged.report &>/dev/null; then
    echo "  Génération Krona combiné merged..."
    ktImportTaxonomy -t 5 -m 3 -o "${KRONA_MERGED_DIR}/all_samples_merged_krona.html" *_merged.report
  fi
fi

for report in *_merged.report; do
  [[ -f "${report}" ]] || continue
  base="$(basename "${report}" ".report")"
  out="${KRONA_MERGED_DIR}/${base}_krona.html"
  if [[ -f "${out}" ]]; then
    echo "  Krona merged déjà présent pour ${base}, on saute."
    continue
  fi
  echo "  Génération Krona pour ${base}..."
  ktImportTaxonomy -t 5 -m 3 -o "${out}" "${report}"
done

echo
echo "===== TAPE 5b : Krona - UNMERGED READS ====="

KRONA_UNMERGED_DIR="${BASEDIR}/09_krona_unmerged_reads"
mkdir -p "${KRONA_UNMERGED_DIR}"
cd "${KRAKEN_UNMERGED_DIR}"

if [[ ! -f "${KRONA_UNMERGED_DIR}/all_samples_unmerged_krona.html" ]]; then
  if ls *_unmerged.report &>/dev/null; then
    echo "  Génération Krona combiné unmerged..."
    ktImportTaxonomy -t 5 -m 3 -o "${KRONA_UNMERGED_DIR}/all_samples_unmerged_krona.html" *_unmerged.report
  fi
fi

for report in *_unmerged.report; do
  [[ -f "${report}" ]] || continue
  base="$(basename "${report}" ".report")"
  out="${KRONA_UNMERGED_DIR}/${base}_krona.html"
  if [[ -f "${out}" ]]; then
    echo "  Krona unmerged déjà présent pour ${base}, on saute."
    continue
  fi
  echo "  Génération Krona pour ${base}..."
  ktImportTaxonomy -t 5 -m 3 -o "${out}" "${report}"
done

############################################
# TAPE 6 : TABLES MPA MERGED & UNMERGED
############################################

echo
echo "===== TAPE 6a : Tables MPA - MERGED READS ====="

if [[ ! -d "${KRAKENTOOLSDIR}" ]]; then
  echo "Installation de KrakenTools..."
  mkdir -p "${BASEDIR}/08_kraken2"
  cd "${BASEDIR}/08_kraken2"
  git clone https://github.com/jenniferlu717/KrakenTools.git
fi

MPA_MERGED_DIR="${BASEDIR}/10_mpa_tables_merged_reads"
mkdir -p "${MPA_MERGED_DIR}"
cd "${KRAKEN_MERGED_DIR}"

if [[ ! -f "${MPA_MERGED_DIR}/combined_all_merged.tsv" ]]; then
  declare -a mpa_files_merged=()
  for report in *_merged.report; do
    [[ -f "${report}" ]] || continue
    base="$(basename "${report}" ".report")"
    mpafile="${MPA_MERGED_DIR}/${base}.mpa"
    if [[ -f "${mpafile}" ]]; then
      echo "  MPA déjà présent pour ${base}, on saute conversion."
    else
      echo "  Conversion de ${base}..."
      python3 "${KRAKENTOOLSDIR}/kreport2mpa.py" -r "${report}" -o "${mpafile}"
    fi
    mpa_files_merged+=("${mpafile}")
  done

  if (( ${#mpa_files_merged[@]} > 0 )); then
    echo "  Combinaison de tous les fichiers MPA merged..."
    python3 "${KRAKENTOOLSDIR}/combine_mpa.py" \
      -i "${mpa_files_merged[@]}" \
      -o "${MPA_MERGED_DIR}/combined_all_merged.tsv"
  fi
fi

echo
echo "===== TAPE 6b : Tables MPA - UNMERGED READS ====="

MPA_UNMERGED_DIR="${BASEDIR}/10_mpa_tables_unmerged_reads"
mkdir -p "${MPA_UNMERGED_DIR}"
cd "${KRAKEN_UNMERGED_DIR}"

if [[ ! -f "${MPA_UNMERGED_DIR}/combined_all_unmerged.tsv" ]]; then
  declare -a mpa_files_unmerged=()
  for report in *_unmerged.report; do
    [[ -f "${report}" ]] || continue
    base="$(basename "${report}" ".report")"
    mpafile="${MPA_UNMERGED_DIR}/${base}.mpa"
    if [[ -f "${mpafile}" ]]; then
      echo "  MPA déjà présent pour ${base}, on saute conversion."
    else
      echo "  Conversion de ${base}..."
      python3 "${KRAKENTOOLSDIR}/kreport2mpa.py" -r "${report}" -o "${mpafile}"
    fi
    mpa_files_unmerged+=("${mpafile}")
  done

  if (( ${#mpa_files_unmerged[@]} > 0 )); then
    echo "  Combinaison de tous les fichiers MPA unmerged..."
    python3 "${KRAKENTOOLSDIR}/combine_mpa.py" \
      -i "${mpa_files_unmerged[@]}" \
      -o "${MPA_UNMERGED_DIR}/combined_all_unmerged.tsv"
  fi
fi

############################################
# TAPE 7 : TABLEAU RÉCAPITULATIF
############################################

echo
echo "===== TAPE 7 : Tableau récapitulatif ====="

SUMMARY_TABLE="${BASEDIR}/11_summary_tables/sequences_summary_merged.tsv"

if [[ ! -f "${SUMMARY_TABLE}" ]]; then
  echo -e "Sample\tType\tNb_sequences\tLongueur_moyenne\tGC_percent" > "${SUMMARY_TABLE}"

  extract_stats () {
    local file="$1"
    local sample="$2"
    local type="$3"

    [[ -f "${file}" ]] || return 0

    local nbseq
    nbseq=$(zcat "${file}" 2>/dev/null | echo $(( $(wc -l) / 4 )))
    if (( nbseq == 0 )); then
      echo -e "${sample}\t${type}\t0\t0\t0" >> "${SUMMARY_TABLE}"
      return 0
    fi

    local stats
    stats=$(zcat "${file}" 2>/dev/null | \
      awk '
        NR%4==2 {
          len=length($0);
          totallen+=len;
          gccount+=gsub(/G|C/,"",$0);
          atcount+=gsub(/A|T/,"",$0);
          count++;
        }
        END {
          if (count>0) {
            avglen=totallen/count;
            gcperc=(gccount/(gccount+atcount))*100;
            printf("%.1f\t%.2f",avglen,gcperc);
          } else {
            printf("0\t0");
          }
        }')
    echo -e "${sample}\t${type}\t${nbseq}\t${stats}" >> "${SUMMARY_TABLE}"
  }

  echo "Calcul des statistiques..."
  for sample in "${ALL_SAMPLES[@]}"; do
    echo "  ${sample}..."

    merged="${FINAL_MERGED_DIR}/${sample}_final_merged.fastq.gz"
    [[ -f "${merged}" ]] && extract_stats "${merged}" "${sample}" "merged"

    unmerged_r1="${FINAL_UNMERGED_DIR}/${sample}_final_unmerged_R1.fastq.gz"
    [[ -f "${unmerged_r1}" ]] && extract_stats "${unmerged_r1}" "${sample}" "unmerged_R1"

    unmerged_r2="${FINAL_UNMERGED_DIR}/${sample}_final_unmerged_R2.fastq.gz"
    [[ -f "${unmerged_r2}" ]] && extract_stats "${unmerged_r2}" "${sample}" "unmerged_R2"
  done

  echo "Tableau récapitulatif créé : ${SUMMARY_TABLE}"
fi

############################################
# TAPE 8 : MAPDAMAGE (UNMERGED READS)
############################################

echo
echo "===== TAPE 8 : MapDamage - Analyse des dommages (UNMERGED) ====="

module load conda/4.12.0 || true
source ~/.bashrc
conda activate mapdamage_py39

echo "Environnement MapDamage activé."

KRAKEN_UNMERGED="${BASEDIR}/08_kraken2_unmerged_reads"
FINAL_UNMERGED="${BASEDIR}/06_fastp_unmerged_reads"

if [[ ! -d "${KRAKEN_UNMERGED}" ]]; then
  echo "ERREUR : dossier ${KRAKEN_UNMERGED} introuvable."
  exit 1
fi
if [[ ! -d "${FINAL_UNMERGED}" ]]; then
  echo "ERREUR : dossier ${FINAL_UNMERGED} introuvable."
  exit 1
fi

echo "${KRAKEN_UNMERGED}"
echo "${FINAL_UNMERGED}"

# Références génomiques (nouveaux chemins)
REF_BASE="/home/plstenge/genomes"

declare -A FIXED_GENOMES
declare -A TAXONS

FIXED_GENOMES["Bos_taurus"]="${REF_BASE}/Bos_taurus/Bos_taurus.ARS-UCD1.3.dna.toplevel.fixed.fa"
TAXONS["Bos_taurus"]="9903"

FIXED_GENOMES["Ovis_aries"]="${REF_BASE}/Ovis_aries/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fixed.fa"
TAXONS["Ovis_aries"]="9940"

FIXED_GENOMES["Capra_hircus"]="${REF_BASE}/Capra_hircus/Capra_hircus.ARS1.dna.toplevel.fixed.fa"
TAXONS["Capra_hircus"]="9925"

# Ajouter ici les autres (Phragmites, Alnus, Corylus) si tu as des .fixed.fa

DAMAGE_UNMERGED_BASE="${BASEDIR}/12_mapdamage_unmerged_reads"
mkdir -p "${DAMAGE_UNMERGED_BASE}"

LOGFILE="${BASEDIR}/00_scripts/mapdamage_$(date +%Y%m%d%H%M%S).txt"
touch "${LOGFILE}"

MAPPINGINFO="${BASEDIR}/11_summary_tables/mapping_bwa_info_unmerged.tsv"
mkdir -p "${BASEDIR}/11_summary_tables"
if [[ ! -f "${MAPPINGINFO}" ]]; then
  echo -e "Sample\tSpecies\tType\tTotalReads\tMappedReads\tMappingRate" > "${MAPPINGINFO}"
fi

echo "$(date) Script MapDamage UNMERGED démarré" | tee -a "${LOGFILE}"

# Fonctions utilitaires
calculate_mapping_rate () {
  local bamfile="$1"
  local samplename="$2"
  local species="$3"
  local type="$4"

  [[ -f "${bamfile}" ]] || return 0

  local total mapped rate
  total=$(samtools view -c "${bamfile}" 2>/dev/null || echo 0)
  mapped=$(samtools view -c -F 4 "${bamfile}" 2>/dev/null || echo 0)
  rate="0.0"
  if (( total > 0 )); then
    rate=$(echo "scale=1; ${mapped}*100/${total}" | bc 2>/dev/null || echo "0.0")
  fi
  echo -e "${samplename}\t${species}\t${type}\t${total}\t${mapped}\t${rate}" >> "${MAPPINGINFO}"
  echo "  Stats ${samplename} ${species} : total=${total}, mapped=${mapped}, rate=${rate}" | tee -a "${LOGFILE}"
}

run_mapdamage_safe () {
  local bamfile="$1"
  local reffasta="$2"
  local outputdir="$3"

  if [[ ! -s "${bamfile}" ]]; then
    echo "  Fichier BAM vide ou inexistant : ${bamfile} (MapDamage sauté)" | tee -a "${LOGFILE}"
    return 0
  fi
  mkdir -p "${outputdir}"
  echo -n "  MapDamage en cours..." | tee -a "${LOGFILE}"
  if mapDamage -i "${bamfile}" -r "${reffasta}" --folder "${outputdir}" --no-stats >> "${LOGFILE}" 2>&1; then
    echo " OK" | tee -a "${LOGFILE}"
  else
    echo " échec (non bloquant)" | tee -a "${LOGFILE}"
  fi
}

shopt -s nullglob

if ls "${KRAKEN_UNMERGED}"/*.kraken &>/dev/null; then
  for KRAKENFILE in "${KRAKEN_UNMERGED}"/*_unmerged.kraken; do
    KRAKENBASE="$(basename "${KRAKENFILE}" ".kraken")"
    SAMPLE="${KRAKENBASE%_unmerged}"

    R1="${FINAL_UNMERGED}/${SAMPLE}_final_unmerged_R1.fastq.gz"
    R2="${FINAL_UNMERGED}/${SAMPLE}_final_unmerged_R2.fastq.gz"

    if [[ ! -f "${R1}" || ! -f "${R2}" || ! -s "${R1}" || ! -s "${R2}" ]]; then
      echo
      echo "  Fichiers unmerged manquants ou vides pour ${SAMPLE}, on saute."
      continue
    fi

    echo
    echo "===== ${SAMPLE} unmerged - $(ls -lh "${R1}" | awk '{print $5}') / $(ls -lh "${R2}" | awk '{print $5}') ====="

    for GROUP in "${!TAXONS[@]}"; do
      REF="${FIXED_GENOMES[${GROUP}]}"
      TAXID="${TAXONS[${GROUP}]}"

      if [[ -z "${REF}" || ! -f "${REF}" ]]; then
        echo "  Génome non disponible pour ${GROUP} (${REF}), on saute." | tee -a "${LOGFILE}"
        continue
      fi

      OUTDIR="${DAMAGE_UNMERGED_BASE}/${GROUP}_${SAMPLE}"
      mkdir -p "${OUTDIR}"
      OUTR1="${OUTDIR}/${KRAKENBASE}_${GROUP}_R1.fastq"
      OUTR2="${OUTDIR}/${KRAKENBASE}_${GROUP}_R2.fastq"
      BAM_SORTED="${OUTDIR}/${KRAKENBASE}_${GROUP}.sorted.bam"

      echo
      echo "  ${GROUP} (TaxID ${TAXID})" | tee -a "${LOGFILE}"

      # Skip complet si BAM déjà existant et indexé + MapDamage déjà fait
      if [[ -f "${BAM_SORTED}" && -f "${BAM_SORTED}.bai" && -d "${OUTDIR}/mapDamage" ]]; then
        echo "    BAM + MapDamage déjà présents pour ${SAMPLE} / ${GROUP}, on saute." | tee -a "${LOGFILE}"
        continue
      fi

      # Extraction des reads (KrakenTools)
      if [[ ! -s "${OUTR1}" || ! -s "${OUTR2}" ]]; then
        rm -f "${OUTR1}" "${OUTR2}" 2>/dev/null || true
        echo "    Extraction reads paired-end..." | tee -a "${LOGFILE}"
        python3 "${KRAKENTOOLSDIR}/extract_kraken_reads.py" \
          -k "${KRAKENFILE}" \
          -s "${R1}" \
          -s2 "${R2}" \
          -t "${TAXID}" \
          -o "${OUTR1}" \
          -o2 "${OUTR2}" \
          --fastq-output >> "${LOGFILE}" 2>&1
      else
        echo "    Fichiers extraits déjà présents pour ${GROUP}/${SAMPLE}, on saute extraction." | tee -a "${LOGFILE}"
      fi

      if [[ ! -s "${OUTR1}" || ! -s "${OUTR2}" ]]; then
        echo "    Aucun read extrait pour ${GROUP}, on nettoie et on continue." | tee -a "${LOGFILE}"
        rm -f "${OUTR1}" "${OUTR2}" 2>/dev/null || true
        continue
      fi

      READCOUNT=$(grep -c "^@" "${OUTR1}" 2>/dev/null || echo 0)
      echo "    ${READCOUNT} paires extraites." | tee -a "${LOGFILE}"

      # Indexation BWA si besoin
      FIXEDREF="${REF}"
      if [[ ! -f "${FIXEDREF}.bwt" ]]; then
        echo "    Indexation BWA (peut prendre du temps)..." | tee -a "${LOGFILE}"
        if timeout 1800 bwa index "${FIXEDREF}" >> "${LOGFILE}" 2>&1; then
          echo "    Index BWA créé." | tee -a "${LOGFILE}"
        else
          echo "    Index BWA incomplet pour ${FIXEDREF}, on continue néanmoins." | tee -a "${LOGFILE}"
        fi
      else
        echo "    Index BWA existant." | tee -a "${LOGFILE}"
      fi

      # Mapping BWA
      if [[ ! -f "${BAM_SORTED}" || ! -f "${BAM_SORTED}.bai" ]]; then
        echo "    Mapping BWA paired-end..." | tee -a "${LOGFILE}"
        bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "${FIXEDREF}" "${OUTR1}" > "${OUTDIR}/${KRAKENBASE}_${GROUP}_R1.sai" 2>> "${LOGFILE}"
        bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "${FIXEDREF}" "${OUTR2}" > "${OUTDIR}/${KRAKENBASE}_${GROUP}_R2.sai" 2>> "${LOGFILE}"
        bwa sampe "${FIXEDREF}" \
          "${OUTDIR}/${KRAKENBASE}_${GROUP}_R1.sai" \
          "${OUTDIR}/${KRAKENBASE}_${GROUP}_R2.sai" \
          "${OUTR1}" "${OUTR2}" 2>> "${LOGFILE}" | \
          samtools view -bS - 2>> "${LOGFILE}" | \
          samtools sort -o "${BAM_SORTED}" - 2>> "${LOGFILE}"

        if [[ ! -f "${BAM_SORTED}" ]]; then
          echo "    Erreur lors du mapping, on saute ${GROUP} pour ${SAMPLE}." | tee -a "${LOGFILE}"
          continue
        fi

        samtools index "${BAM_SORTED}" 2>> "${LOGFILE}"
        echo "    Mapping terminé." | tee -a "${LOGFILE}"
      else
        echo "    BAM déjà existant pour ${GROUP}/${SAMPLE}, on saute mapping." | tee -a "${LOGFILE}"
      fi

      calculate_mapping_rate "${BAM_SORTED}" "${SAMPLE}" "${GROUP}" "unmerged"

      # MapDamage
      run_mapdamage_safe "${BAM_SORTED}" "${FIXEDREF}" "${OUTDIR}/mapDamage"

      # Nettoyage minimal (garder BAM+index)
      rm -f "${OUTDIR}/${KRAKENBASE}_${GROUP}_R1.sai" \
            "${OUTDIR}/${KRAKENBASE}_${GROUP}_R2.sai" \
            "${OUTR1}" "${OUTR2}" 2>/dev/null || true
    done
  done
else
  echo "ERREUR : Aucun fichier .kraken trouvé dans ${KRAKEN_UNMERGED}."
  exit 1
fi

shopt -u nullglob

echo
echo "===== PIPELINE MAPDAMAGE TERMINÉ ====="
echo "Date : $(date)"
echo
echo "Statistiques de mapping : ${MAPPINGINFO}"
echo "Log complet : ${LOGFILE}"
echo "Fichiers BAM MapDamage : ${DAMAGE_UNMERGED_BASE}"
echo

mail -s "Pipeline Coprolites MapDamage - Terminé" pierrelouis.stenger@gmail.com < "${LOGFILE}" 2>/dev/null || true

exit 0
