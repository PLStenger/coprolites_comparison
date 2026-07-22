#!/bin/bash

#SBATCH --job-name=50_per_run_pipeline
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=36
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/50_per_run_pipeline.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/50_per_run_pipeline.out"

# ==============================================================================
# Script : 11_per_run_pipeline.sh
# Author : Pierre-Louis Stenger
# Purpose :
#   Reproduire EXACTEMENT le pipeline des scripts 01 -> 08 (BBDuk, FastUniq,
#   Clumpify, fastp, Kraken2 k25/k29/k35, Bracken, export MPA), mais SANS
#   concatener les runs entre eux. Chaque run (RUN1..RUN5) et chaque
#   echantillon (cop408, cop410, cop412, cop414, cop417) est traite comme une
#   unite independante, de bout en bout, avec les memes parametres que les
#   scripts d'origine.
#
#   Tous les fichiers de sortie encodent systematiquement :
#     <sample>_<run_label>_<merged|unmerged>
#   pour ne jamais confondre avec les anciens resultats "grouped" (concatenes).
#
#   Nouvelle arborescence de sortie, isolee des anciens resultats :
#     11_per_run_analysis/
#       01_bbduk/<run>/<sample>/...
#       02_fastuniq/<run>/<sample>/...
#       03_clumpify/<run>/<sample>/...
#       04_fastp/<run>/<sample>/...
#       05_kraken2_k25/<run>/<sample>/...
#       05_kraken2_k29/<run>/<sample>/...
#       05_kraken2_k35/<run>/<sample>/...
#       06_bracken/k25|k29|k35/<run>/<sample>/...
#       07_mpa_tables/k25|k29|k35/merged|unmerged/...
# ==============================================================================

#set -uo pipefail

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics

# ==============================================================================
# PATHS — Input run directories (identiques au script 01_group_cop_fastqc.sh)
# ==============================================================================

RUN1="/storage/groups/gdec/shared_paleo/Illumina/Coprolites_Illumina"
RUN2="/storage/groups/gdec/shared_paleo/E1531_final/run1_20250320_AV241601_E1531_Ps5Lane1_Ps6Lane2"
RUN3="/storage/groups/gdec/shared_paleo/E1531_final/run2_20250414_AV241601_E1531_Ps5_Ps6_11022026_CORRECTED"
RUN4="/storage/groups/gdec/shared_paleo/E1531_final/run3_20251008_AV241601_E1531_Ps7_Ps8"
RUN5="/storage/groups/gdec/shared_paleo/E1531_final/run4_20251104_AV241601_E1531_Ps7_Ps8_04112025"

RUN_LABELS=("run1" "run2" "run3" "run4" "run5")
ALL_RUNS=("${RUN1}" "${RUN2}" "${RUN3}" "${RUN4}" "${RUN5}")

SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")

# ==============================================================================
# OUTPUT ROOT — dossier dedie, pour ne jamais melanger avec les anciens
# resultats "grouped" (concatenes)
# ==============================================================================

WORKDIR="/home/plstenge/coprolites_comparison"
OUT_ROOT="${WORKDIR}/11_per_run_analysis"

BBDUK_DIR="${OUT_ROOT}/01_bbduk"
FASTUNIQ_DIR="${OUT_ROOT}/02_fastuniq"
CLUMPIFY_DIR="${OUT_ROOT}/03_clumpify"
FASTP_DIR="${OUT_ROOT}/04_fastp"
KRAKEN_K25_DIR="${OUT_ROOT}/05_kraken2_k25"
KRAKEN_K29_DIR="${OUT_ROOT}/05_kraken2_k29"
KRAKEN_K35_DIR="${OUT_ROOT}/05_kraken2_k35"
BRACKEN_DIR="${OUT_ROOT}/06_bracken"
MPA_DIR="${OUT_ROOT}/07_mpa_tables"

KRAKENTOOLS_DIR="${WORKDIR}/08_krakentools/KrakenTools"

THREADS="${SLURM_CPUS_PER_TASK:-36}"

# ==============================================================================
# TOOLS / DB (identiques aux scripts d'origine)
# ==============================================================================

BBDUK="/home/plstenge/bbmap/bbduk.sh"
PHIX="/home/plstenge/bbmap/resources/phix174_ill.ref.fa.gz"
CLUMPIFY="clumpify.sh"

KRAKEN2_DB_K25="/storage/groups/gdec/shared/Kraken_database/core_nt_k25"
KRAKEN2_DB_K29="/storage/groups/gdec/shared/Kraken_database/core_nt_k29"
KRAKEN2_DB_K35="/storage/groups/gdec/shared/Kraken_database/k2_core_nt_20250609"

BRACKEN_READ_LEN=50
BRACKEN_LEVEL="S"
BRACKEN_THRESHOLD=10

# ==============================================================================
# HELPER : find_fastq (identique au script 01)
# ==============================================================================

find_fastq() {
  local run_dir="$1"
  local sample="$2"
  local read="$3"
  find "${run_dir}" -maxdepth 2 -type f \
    -name "*${sample}*${read}*.fastq.gz" \
    2>/dev/null | head -n 1
}

# ==============================================================================
# STEP 0 : Creer l'arborescence de sortie
# ==============================================================================

echo "======================================================================="
echo " STEP 0: Creating dedicated output directories under ${OUT_ROOT}"
echo "======================================================================="

mkdir -p "${BBDUK_DIR}" "${FASTUNIQ_DIR}" "${CLUMPIFY_DIR}" "${FASTP_DIR}"
mkdir -p "${KRAKEN_K25_DIR}" "${KRAKEN_K29_DIR}" "${KRAKEN_K35_DIR}"
mkdir -p "${BRACKEN_DIR}/k25" "${BRACKEN_DIR}/k29" "${BRACKEN_DIR}/k35"
mkdir -p "${MPA_DIR}/k25/merged" "${MPA_DIR}/k25/unmerged"
mkdir -p "${MPA_DIR}/k29/merged" "${MPA_DIR}/k29/unmerged"
mkdir -p "${MPA_DIR}/k35/merged" "${MPA_DIR}/k35/unmerged"

if [[ ! -d "${KRAKENTOOLS_DIR}" ]]; then
  echo " KrakenTools non trouve, installation..."
  mkdir -p "${WORKDIR}/08_krakentools"
  ( cd "${WORKDIR}/08_krakentools" && git clone https://github.com/jenniferlu717/KrakenTools.git )
else
  echo " KrakenTools deja present dans : ${KRAKENTOOLS_DIR}"
fi

echo " Done."

# ==============================================================================
# STEP 1 : Boucle principale — pour chaque RUN x SAMPLE, pipeline complet
#          BBDuk -> FastUniq -> Clumpify -> fastp -> Kraken2 (k25/k29/k35)
# ==============================================================================

echo "======================================================================="
echo " STEP 1: Per-run, per-sample pipeline (BBDuk -> FastUniq -> Clumpify -> fastp -> Kraken2)"
echo "======================================================================="

for i in "${!ALL_RUNS[@]}"; do
  RUN_DIR="${ALL_RUNS[$i]}"
  RUN_LABEL="${RUN_LABELS[$i]}"

  for SAMPLE in "${SAMPLES[@]}"; do

    R1_RAW=$(find_fastq "${RUN_DIR}" "${SAMPLE}" "R1")
    R2_RAW=$(find_fastq "${RUN_DIR}" "${SAMPLE}" "R2")

    if [[ -z "${R1_RAW}" || -z "${R2_RAW}" ]]; then
      echo " [SKIP] ${SAMPLE} / ${RUN_LABEL} : fichier(s) R1/R2 non trouve(s) dans ${RUN_DIR}"
      continue
    fi

    UNIT="${SAMPLE}_${RUN_LABEL}"
    echo ""
    echo "-----------------------------------------------------------------"
    echo " Unite de traitement : ${UNIT}"
    echo " R1 brut : ${R1_RAW}"
    echo " R2 brut : ${R2_RAW}"
    echo "-----------------------------------------------------------------"

    # ------------------------------------------------------------------
    # 1a. BBDuk — decontamination phiX + trimming (memes parametres
    #     que 02_bbduk-2.sh : ktrim=rl k=23 mink=11 hdist=1)
    # ------------------------------------------------------------------
    BBDUK_OUT="${BBDUK_DIR}/${RUN_LABEL}/${SAMPLE}"
    mkdir -p "${BBDUK_OUT}"

    CLEAN_R1="${BBDUK_OUT}/clean_${UNIT}_R1.fastq.gz"
    CLEAN_R2="${BBDUK_OUT}/clean_${UNIT}_R2.fastq.gz"

    echo " [BBDuk] ${UNIT}"
    "${BBDUK}" \
      in1="${R1_RAW}" \
      in2="${R2_RAW}" \
      out1="${CLEAN_R1}" \
      out2="${CLEAN_R2}" \
      ref="${PHIX}" \
      ktrim=rl \
      k=23 \
      mink=11 \
      hdist=1 \
      threads="${THREADS}" \
      stats="${BBDUK_OUT}/${UNIT}_bbduk_stats.txt"

    if [[ ! -s "${CLEAN_R1}" || ! -s "${CLEAN_R2}" ]]; then
      echo " ERREUR: BBDuk n'a pas produit de sortie valide pour ${UNIT}, on saute."
      continue
    fi

    # ------------------------------------------------------------------
    # 1b. FastUniq — deduplication (memes parametres que 03_fastuniq-3.sh)
    #     Necessite une decompression temporaire (.fastq non compresse)
    # ------------------------------------------------------------------
    FASTUNIQ_OUT="${FASTUNIQ_DIR}/${RUN_LABEL}/${SAMPLE}"
    FASTUNIQ_TMP="${FASTUNIQ_OUT}/tmp"
    mkdir -p "${FASTUNIQ_OUT}" "${FASTUNIQ_TMP}"

    R1_TMP="${FASTUNIQ_TMP}/${UNIT}_R1.fastq"
    R2_TMP="${FASTUNIQ_TMP}/${UNIT}_R2.fastq"
    LISTFILE="${FASTUNIQ_TMP}/${UNIT}.list"

    echo " [FastUniq] decompression temporaire pour ${UNIT}"
    gzip -dc "${CLEAN_R1}" > "${R1_TMP}"
    gzip -dc "${CLEAN_R2}" > "${R2_TMP}"
    echo -e "${R1_TMP}\n${R2_TMP}" > "${LISTFILE}"

    DEDUP_R1="${FASTUNIQ_OUT}/clean_${UNIT}_dedup_R1.fastq"
    DEDUP_R2="${FASTUNIQ_OUT}/clean_${UNIT}_dedup_R2.fastq"

    echo " [FastUniq] deduplication pour ${UNIT}"
    fastuniq \
      -i "${LISTFILE}" \
      -t q \
      -o "${DEDUP_R1}" \
      -p "${DEDUP_R2}"

    rm -f "${R1_TMP}" "${R2_TMP}" "${LISTFILE}"
    rmdir "${FASTUNIQ_TMP}" 2>/dev/null || true

    if [[ ! -s "${DEDUP_R1}" || ! -s "${DEDUP_R2}" ]]; then
      echo " ERREUR: FastUniq n'a pas produit de sortie valide pour ${UNIT}, on saute."
      continue
    fi

    # ------------------------------------------------------------------
    # 1c. Clumpify — dedupe=t (memes parametres que 04_clumpify-4.sh)
    # ------------------------------------------------------------------
    CLUMPIFY_OUT="${CLUMPIFY_DIR}/${RUN_LABEL}/${SAMPLE}"
    mkdir -p "${CLUMPIFY_OUT}"

    CLUMP_R1="${CLUMPIFY_OUT}/clean_${UNIT}_dedup_clumpify_R1.fastq.gz"
    CLUMP_R2="${CLUMPIFY_OUT}/clean_${UNIT}_dedup_clumpify_R2.fastq.gz"

    echo " [Clumpify] ${UNIT}"
    "${CLUMPIFY}" \
      in="${DEDUP_R1}" \
      in2="${DEDUP_R2}" \
      out="${CLUMP_R1}" \
      out2="${CLUMP_R2}" \
      dedupe=t \
      threads="${THREADS}"

    if [[ ! -s "${CLUMP_R1}" || ! -s "${CLUMP_R2}" ]]; then
      echo " ERREUR: Clumpify n'a pas produit de sortie valide pour ${UNIT}, on saute."
      continue
    fi

    # ------------------------------------------------------------------
    # 1d. fastp — trimming + merge (memes parametres que 05_fastp-5.sh)
    # ------------------------------------------------------------------
    FASTP_OUT="${FASTP_DIR}/${RUN_LABEL}/${SAMPLE}"
    mkdir -p "${FASTP_OUT}"

    FASTP_UNMERGED_R1="${FASTP_OUT}/${UNIT}_fastp_unmerged_R1.fastq.gz"
    FASTP_UNMERGED_R2="${FASTP_OUT}/${UNIT}_fastp_unmerged_R2.fastq.gz"
    FASTP_MERGED="${FASTP_OUT}/${UNIT}_fastp_merged.fastq.gz"
    FASTP_HTML="${FASTP_OUT}/${UNIT}_fastp_report.html"
    FASTP_JSON="${FASTP_OUT}/${UNIT}_fastp_report.json"

    echo " [fastp] ${UNIT}"
    fastp \
      --in1 "${CLUMP_R1}" --in2 "${CLUMP_R2}" \
      --out1 "${FASTP_UNMERGED_R1}" --out2 "${FASTP_UNMERGED_R2}" \
      --merged_out "${FASTP_MERGED}" \
      --length_required 20 \
      --cut_front --cut_tail \
      --cut_window_size 4 \
      --cut_mean_quality 10 \
      --n_base_limit 5 \
      --unqualified_percent_limit 40 \
      --complexity_threshold 30 \
      --qualified_quality_phred 20 \
      --low_complexity_filter \
      --trim_poly_x \
      --poly_x_min_len 10 \
      --merge --correction \
      --overlap_len_require 10 \
      --overlap_diff_limit 5 \
      --overlap_diff_percent_limit 20 \
      --html "${FASTP_HTML}" \
      --json "${FASTP_JSON}" \
      --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
      --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
      --detect_adapter_for_pe \
      --thread "${THREADS}"

    # ------------------------------------------------------------------
    # 1e. Kraken2 — k25, k29, k35, sur merged (single-end) ET unmerged
    #     (paired-end), memes parametres --conf 0.2 que les scripts 06_*
    # ------------------------------------------------------------------
    for KRAKEN_SET in "k25:${KRAKEN2_DB_K25}:${KRAKEN_K25_DIR}" \
                       "k29:${KRAKEN2_DB_K29}:${KRAKEN_K29_DIR}" \
                       "k35:${KRAKEN2_DB_K35}:${KRAKEN_K35_DIR}"; do

      IFS=":" read -r KLABEL KDB KOUTROOT <<< "${KRAKEN_SET}"
      KOUT="${KOUTROOT}/${RUN_LABEL}/${SAMPLE}"
      mkdir -p "${KOUT}"

      if [[ -s "${FASTP_MERGED}" ]]; then
        echo " [Kraken2 ${KLABEL}] ${UNIT} (merged)"
        kraken2 --conf 0.2 --db "${KDB}" --threads "${THREADS}" \
          --output "${KOUT}/${UNIT}_merged.kraken" \
          --report "${KOUT}/${UNIT}_merged.report" \
          "${FASTP_MERGED}"
      fi

      if [[ -s "${FASTP_UNMERGED_R1}" && -s "${FASTP_UNMERGED_R2}" ]]; then
        echo " [Kraken2 ${KLABEL}] ${UNIT} (unmerged)"
        kraken2 --conf 0.2 --paired --db "${KDB}" --threads "${THREADS}" \
          --output "${KOUT}/${UNIT}_unmerged.kraken" \
          --report "${KOUT}/${UNIT}_unmerged.report" \
          "${FASTP_UNMERGED_R1}" "${FASTP_UNMERGED_R2}"
      fi

    done

    echo " >>> Pipeline complet termine pour ${UNIT} <<<"

  done
done

echo ""
echo " STEP 1 done. Tous les runs/echantillons ont ete traites individuellement."

# ==============================================================================
# STEP 2 : Bracken — sur chaque report Kraken2 (k25/k29/k35), merged/unmerged
#          separement, run par run (memes parametres que 08_bracken_krona_mpa-10.sh)
# ==============================================================================

echo "======================================================================="
echo " STEP 2: Bracken par run x echantillon x base k-mer"
echo "======================================================================="

run_bracken_for_kmer() {
  local KRAKEN_ROOT="$1"
  local DB="$2"
  local OUT_ROOT_K="$3"
  local KLABEL="$4"

  for i in "${!ALL_RUNS[@]}"; do
    RUN_LABEL="${RUN_LABELS[$i]}"

    for SAMPLE in "${SAMPLES[@]}"; do
      UNIT="${SAMPLE}_${RUN_LABEL}"
      KRAKEN_DIR="${KRAKEN_ROOT}/${RUN_LABEL}/${SAMPLE}"
      OUT_DIR="${OUT_ROOT_K}/${RUN_LABEL}/${SAMPLE}"
      mkdir -p "${OUT_DIR}"

      MERGED_REPORT="${KRAKEN_DIR}/${UNIT}_merged.report"
      UNMERGED_REPORT="${KRAKEN_DIR}/${UNIT}_unmerged.report"

      if [[ -f "${MERGED_REPORT}" ]]; then
        echo " [Bracken ${KLABEL}] ${UNIT} (merged)"
        bracken \
          -d "${DB}" \
          -i "${MERGED_REPORT}" \
          -o "${OUT_DIR}/${UNIT}_merged.bracken" \
          -w "${OUT_DIR}/${UNIT}_merged.bracken.report" \
          -r "${BRACKEN_READ_LEN}" \
          -l "${BRACKEN_LEVEL}" \
          -t "${BRACKEN_THRESHOLD}"
      fi

      if [[ -f "${UNMERGED_REPORT}" ]]; then
        echo " [Bracken ${KLABEL}] ${UNIT} (unmerged)"
        bracken \
          -d "${DB}" \
          -i "${UNMERGED_REPORT}" \
          -o "${OUT_DIR}/${UNIT}_unmerged.bracken" \
          -w "${OUT_DIR}/${UNIT}_unmerged.bracken.report" \
          -r "${BRACKEN_READ_LEN}" \
          -l "${BRACKEN_LEVEL}" \
          -t "${BRACKEN_THRESHOLD}"
      fi
    done
  done
}

run_bracken_for_kmer "${KRAKEN_K25_DIR}" "${KRAKEN2_DB_K25}" "${BRACKEN_DIR}/k25" "k25"
run_bracken_for_kmer "${KRAKEN_K29_DIR}" "${KRAKEN2_DB_K29}" "${BRACKEN_DIR}/k29" "k29"
run_bracken_for_kmer "${KRAKEN_K35_DIR}" "${KRAKEN2_DB_K35}" "${BRACKEN_DIR}/k35" "k35"

echo ""
echo " STEP 2 done."

# ==============================================================================
# STEP 3 : Export MPA + tables combinees (par base k-mer, merged/unmerged
#          separement, run par run) — memes outils que 07_mpa_tables-9.sh
#
#          Les en-tetes de colonnes des tables combinees porteront le nom
#          des fichiers .mpa, donc : <sample>_<run>_<merged|unmerged>.mpa
# ==============================================================================

echo "======================================================================="
echo " STEP 3: Export MPA + tables combinees par base k-mer"
echo "======================================================================="

export_mpa_for_kmer() {
  local KRAKEN_ROOT="$1"
  local BRACKEN_ROOT="$2"
  local MPA_OUT_ROOT="$3"
  local KLABEL="$4"

  local MPA_MERGED_DIR="${MPA_OUT_ROOT}/merged"
  local MPA_UNMERGED_DIR="${MPA_OUT_ROOT}/unmerged"
  mkdir -p "${MPA_MERGED_DIR}" "${MPA_UNMERGED_DIR}"

  declare -a merged_mpa_files=()
  declare -a unmerged_mpa_files=()

  for i in "${!ALL_RUNS[@]}"; do
    RUN_LABEL="${RUN_LABELS[$i]}"

    for SAMPLE in "${SAMPLES[@]}"; do
      UNIT="${SAMPLE}_${RUN_LABEL}"
      BRACKEN_DIR_UNIT="${BRACKEN_ROOT}/${RUN_LABEL}/${SAMPLE}"

      MERGED_BRACKEN_REPORT="${BRACKEN_DIR_UNIT}/${UNIT}_merged.bracken.report"
      UNMERGED_BRACKEN_REPORT="${BRACKEN_DIR_UNIT}/${UNIT}_unmerged.bracken.report"

      if [[ -f "${MERGED_BRACKEN_REPORT}" ]]; then
        MPA_FILE="${MPA_MERGED_DIR}/${UNIT}_merged.mpa"
        echo " [MPA ${KLABEL}] ${UNIT} (merged)"
        python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" -r "${MERGED_BRACKEN_REPORT}" -o "${MPA_FILE}"
        merged_mpa_files+=("${MPA_FILE}")
      fi

      if [[ -f "${UNMERGED_BRACKEN_REPORT}" ]]; then
        MPA_FILE="${MPA_UNMERGED_DIR}/${UNIT}_unmerged.mpa"
        echo " [MPA ${KLABEL}] ${UNIT} (unmerged)"
        python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" -r "${UNMERGED_BRACKEN_REPORT}" -o "${MPA_FILE}"
        unmerged_mpa_files+=("${MPA_FILE}")
      fi
    done
  done

  if [[ ${#merged_mpa_files[@]} -gt 0 ]]; then
    echo " -> Combinaison table MPA merged pour ${KLABEL}"
    python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" -i "${merged_mpa_files[@]}" -o "${MPA_MERGED_DIR}/combined_per_run_merged.tsv"
  fi

  if [[ ${#unmerged_mpa_files[@]} -gt 0 ]]; then
    echo " -> Combinaison table MPA unmerged pour ${KLABEL}"
    python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" -i "${unmerged_mpa_files[@]}" -o "${MPA_UNMERGED_DIR}/combined_per_run_unmerged.tsv"
  fi
}

export_mpa_for_kmer "${KRAKEN_K25_DIR}" "${BRACKEN_DIR}/k25" "${MPA_DIR}/k25" "k25"
export_mpa_for_kmer "${KRAKEN_K29_DIR}" "${BRACKEN_DIR}/k29" "${MPA_DIR}/k29" "k29"
export_mpa_for_kmer "${KRAKEN_K35_DIR}" "${BRACKEN_DIR}/k35" "${MPA_DIR}/k35" "k35"

# ==============================================================================
# RESUME FINAL
# ==============================================================================

echo ""
echo "======================================================================="
echo " ALL STEPS COMPLETED — Pipeline par run, sans concatenation"
echo "======================================================================="
echo " Racine des resultats     : ${OUT_ROOT}"
echo " BBDuk                    : ${BBDUK_DIR}/<run>/<sample>/"
echo " FastUniq                 : ${FASTUNIQ_DIR}/<run>/<sample>/"
echo " Clumpify                 : ${CLUMPIFY_DIR}/<run>/<sample>/"
echo " fastp                    : ${FASTP_DIR}/<run>/<sample>/"
echo " Kraken2 k25/k29/k35      : ${KRAKEN_K25_DIR}, ${KRAKEN_K29_DIR}, ${KRAKEN_K35_DIR}"
echo " Bracken                  : ${BRACKEN_DIR}/k25|k29|k35/<run>/<sample>/"
echo " Tables MPA combinees     : ${MPA_DIR}/k25|k29|k35/merged|unmerged/combined_per_run_*.tsv"
echo ""
echo " Convention de nommage    : <sample>_<run> puis suffixe merged/unmerged"
echo "   Exemple : cop408_run1_merged.report / cop408_run3_unmerged.bracken.report"
echo "======================================================================="
