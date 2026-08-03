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
# Script : 50_per_run_pipeline_fixed.sh
# Author : Pierre-Louis Stenger
# Purpose :
#   Pipeline complet BBDuk -> FastUniq -> Clumpify -> fastp -> Kraken2 -> Bracken
#   -> export MPA + table combinee, execute independamment pour chaque RUN
#   (RUN1..RUN5) x SAMPLE (cop408, cop410, cop412, cop414, cop417), avec QUATRE
#   bases de classification en parallele :
#     - k25   : core_nt_k25
#     - k29   : core_nt_k29
#     - k35   : k2_core_nt_20250609
#     - pracken : k2_NCBI_reference_20251007 (base Pracken pre-construite,
#                 fichiers kmer_distrib deja fournis pour Bracken)
#
#   Corrections apportees par rapport a 50_per_run_pipeline-3.sh :
#     - fonction find_fastq() correctement fermee (bug d'origine)
#     - STEP 0 correctement isole de la boucle principale (accolades fermees)
#     - garde-fous Bracken repris de 50b_fixed_bug.sh (skip proprement les
#       reports absents au lieu de planter tout le job)
#     - set -uo pipefail actif, avec set +e localise autour des appels externes
#       potentiellement absents (bracken, kraken2) pour ne jamais perdre tout
#       le run sur un seul echec
#     - integration de l'etape 51_mpa_all (export MPA global + table combinee
#       avec header informatif) directement a la fin du meme script
#     - ajout de la base Pracken comme quatrieme "kmer set" complet
#       (Kraken2 + Bracken + MPA), au meme titre que k25/k29/k35
#
#   Arborescence de sortie (isolee des anciens resultats "grouped") :
#     11_per_run_analysis/
#       01_bbduk/<run>/<sample>/...
#       02_fastuniq/<run>/<sample>/...
#       03_clumpify/<run>/<sample>/...
#       04_fastp/<run>/<sample>/...
#       05_kraken2_k25/<run>/<sample>/...
#       05_kraken2_k29/<run>/<sample>/...
#       05_kraken2_k35/<run>/<sample>/...
#       05_kraken2_pracken/<run>/<sample>/...
#       06_bracken/k25|k29|k35|pracken/<run>/<sample>/...
#       07_mpa_tables/k25|k29|k35|pracken/merged|unmerged/...
#       07_mpa_tables/combined_all_k25_k29_k35_pracken_merged_unmerged_bracken_kraken.tsv
# ==============================================================================

set -uo pipefail

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
# OUTPUT ROOT
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
KRAKEN_PRACKEN_DIR="${OUT_ROOT}/05_kraken2_pracken"

BRACKEN_DIR="${OUT_ROOT}/06_bracken"
MPA_DIR="${OUT_ROOT}/07_mpa_tables"

KRAKENTOOLS_DIR="${WORKDIR}/08_krakentools/KrakenTools"

THREADS="${SLURM_CPUS_PER_TASK:-36}"

# ==============================================================================
# TOOLS / DB
# ==============================================================================

BBDUK="/home/plstenge/bbmap/bbduk.sh"
PHIX="/home/plstenge/bbmap/resources/phix174_ill.ref.fa.gz"
CLUMPIFY="clumpify.sh"

KRAKEN2_DB_K25="/storage/groups/gdec/shared/Kraken_database/core_nt_k25"
KRAKEN2_DB_K29="/storage/groups/gdec/shared/Kraken_database/core_nt_k29"
KRAKEN2_DB_K35="/storage/groups/gdec/shared/Kraken_database/k2_core_nt_20250609"
KRAKEN2_DB_PRACKEN="/storage/groups/gdec/shared/Kraken_database/k2_NCBI_reference_20251007"

# Parametres Bracken classiques (k25/k29/k35), read length fixe a 50
BRACKEN_READ_LEN=50
BRACKEN_LEVEL="S"
BRACKEN_THRESHOLD=10

# Parametres Bracken pour la base Pracken : celle-ci fournit deja plusieurs
# fichiers database<N>mers.kmer_distrib (50/75/100/150/200/250/300). On choisit
# ici 50 par defaut pour rester coherent avec k25/k29/k35, mais la base
# permet de recalculer facilement pour d'autres longueurs de reads si besoin
# (il suffit d'ajouter la longueur a PRACKEN_READ_LENS).
PRACKEN_READ_LENS=(50)
PRACKEN_LEVEL="S"
PRACKEN_THRESHOLD=10

# ==============================================================================
# HELPER : find_fastq
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
mkdir -p "${KRAKEN_K25_DIR}" "${KRAKEN_K29_DIR}" "${KRAKEN_K35_DIR}" "${KRAKEN_PRACKEN_DIR}"
mkdir -p "${BRACKEN_DIR}/k25" "${BRACKEN_DIR}/k29" "${BRACKEN_DIR}/k35" "${BRACKEN_DIR}/pracken"
mkdir -p "${MPA_DIR}/k25/merged" "${MPA_DIR}/k25/unmerged"
mkdir -p "${MPA_DIR}/k29/merged" "${MPA_DIR}/k29/unmerged"
mkdir -p "${MPA_DIR}/k35/merged" "${MPA_DIR}/k35/unmerged"
mkdir -p "${MPA_DIR}/pracken/merged" "${MPA_DIR}/pracken/unmerged"

if [[ ! -d "${KRAKENTOOLS_DIR}" ]]; then
    echo "  KrakenTools non trouve, installation..."
    mkdir -p "${WORKDIR}/08_krakentools"
    ( cd "${WORKDIR}/08_krakentools" && git clone https://github.com/jenniferlu717/KrakenTools.git )
else
    echo "  KrakenTools deja present dans : ${KRAKENTOOLS_DIR}"
fi

# Verification que les bases existent avant de lancer quoi que ce soit
for DBCHECK in "${KRAKEN2_DB_K25}" "${KRAKEN2_DB_K29}" "${KRAKEN2_DB_K35}" "${KRAKEN2_DB_PRACKEN}"; do
    if [[ ! -d "${DBCHECK}" ]]; then
        echo "  ATTENTION: base introuvable -> ${DBCHECK}"
    fi
done

echo "  Done."

# ==============================================================================
# STEP 1 : Boucle principale — pour chaque RUN x SAMPLE, pipeline complet
#          BBDuk -> FastUniq -> Clumpify -> fastp -> Kraken2 (k25/k29/k35/pracken)
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
            echo "  [SKIP] ${SAMPLE} / ${RUN_LABEL} : fichier(s) R1/R2 non trouve(s) dans ${RUN_DIR}"
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
        # 1a. BBDuk — decontamination phiX + trimming
        # ------------------------------------------------------------------
        BBDUK_OUT="${BBDUK_DIR}/${RUN_LABEL}/${SAMPLE}"
        mkdir -p "${BBDUK_OUT}"

        CLEAN_R1="${BBDUK_OUT}/clean_${UNIT}_R1.fastq.gz"
        CLEAN_R2="${BBDUK_OUT}/clean_${UNIT}_R2.fastq.gz"

        echo "  [BBDuk] ${UNIT}"
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
            echo "  ERREUR: BBDuk n'a pas produit de sortie valide pour ${UNIT}, on saute."
            continue
        fi

        # ------------------------------------------------------------------
        # 1b. FastUniq — deduplication
        # ------------------------------------------------------------------
        FASTUNIQ_OUT="${FASTUNIQ_DIR}/${RUN_LABEL}/${SAMPLE}"
        FASTUNIQ_TMP="${FASTUNIQ_OUT}/tmp"
        mkdir -p "${FASTUNIQ_OUT}" "${FASTUNIQ_TMP}"

        R1_TMP="${FASTUNIQ_TMP}/${UNIT}_R1.fastq"
        R2_TMP="${FASTUNIQ_TMP}/${UNIT}_R2.fastq"
        LISTFILE="${FASTUNIQ_TMP}/${UNIT}.list"

        echo "  [FastUniq] decompression temporaire pour ${UNIT}"
        gzip -dc "${CLEAN_R1}" > "${R1_TMP}"
        gzip -dc "${CLEAN_R2}" > "${R2_TMP}"
        echo -e "${R1_TMP}\n${R2_TMP}" > "${LISTFILE}"

        DEDUP_R1="${FASTUNIQ_OUT}/clean_${UNIT}_dedup_R1.fastq"
        DEDUP_R2="${FASTUNIQ_OUT}/clean_${UNIT}_dedup_R2.fastq"

        echo "  [FastUniq] deduplication pour ${UNIT}"
        fastuniq \
            -i "${LISTFILE}" \
            -t q \
            -o "${DEDUP_R1}" \
            -p "${DEDUP_R2}"

        rm -f "${R1_TMP}" "${R2_TMP}" "${LISTFILE}"
        rmdir "${FASTUNIQ_TMP}" 2>/dev/null || true

        if [[ ! -s "${DEDUP_R1}" || ! -s "${DEDUP_R2}" ]]; then
            echo "  ERREUR: FastUniq n'a pas produit de sortie valide pour ${UNIT}, on saute."
            continue
        fi

        # ------------------------------------------------------------------
        # 1c. Clumpify — dedupe=t
        # ------------------------------------------------------------------
        CLUMPIFY_OUT="${CLUMPIFY_DIR}/${RUN_LABEL}/${SAMPLE}"
        mkdir -p "${CLUMPIFY_OUT}"

        CLUMP_R1="${CLUMPIFY_OUT}/clean_${UNIT}_dedup_clumpify_R1.fastq.gz"
        CLUMP_R2="${CLUMPIFY_OUT}/clean_${UNIT}_dedup_clumpify_R2.fastq.gz"

        echo "  [Clumpify] ${UNIT}"
        "${CLUMPIFY}" \
            in="${DEDUP_R1}" \
            in2="${DEDUP_R2}" \
            out="${CLUMP_R1}" \
            out2="${CLUMP_R2}" \
            dedupe=t \
            threads="${THREADS}"

        if [[ ! -s "${CLUMP_R1}" || ! -s "${CLUMP_R2}" ]]; then
            echo "  ERREUR: Clumpify n'a pas produit de sortie valide pour ${UNIT}, on saute."
            continue
        fi

        # ------------------------------------------------------------------
        # 1d. fastp — trimming + merge
        # ------------------------------------------------------------------
        FASTP_OUT="${FASTP_DIR}/${RUN_LABEL}/${SAMPLE}"
        mkdir -p "${FASTP_OUT}"

        FASTP_UNMERGED_R1="${FASTP_OUT}/${UNIT}_fastp_unmerged_R1.fastq.gz"
        FASTP_UNMERGED_R2="${FASTP_OUT}/${UNIT}_fastp_unmerged_R2.fastq.gz"
        FASTP_MERGED="${FASTP_OUT}/${UNIT}_fastp_merged.fastq.gz"
        FASTP_HTML="${FASTP_OUT}/${UNIT}_fastp_report.html"
        FASTP_JSON="${FASTP_OUT}/${UNIT}_fastp_report.json"

        echo "  [fastp] ${UNIT}"
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
        # 1e. Kraken2 — k25, k29, k35, pracken, sur merged (single-end) ET
        #     unmerged (paired-end), --conf 0.2
        # ------------------------------------------------------------------
        for KRAKEN_SET in "k25:${KRAKEN2_DB_K25}:${KRAKEN_K25_DIR}" \
                          "k29:${KRAKEN2_DB_K29}:${KRAKEN_K29_DIR}" \
                          "k35:${KRAKEN2_DB_K35}:${KRAKEN_K35_DIR}" \
                          "pracken:${KRAKEN2_DB_PRACKEN}:${KRAKEN_PRACKEN_DIR}"; do

            IFS=":" read -r KLABEL KDB KOUTROOT <<< "${KRAKEN_SET}"

            if [[ ! -d "${KDB}" ]]; then
                echo "  [SKIP Kraken2 ${KLABEL}] base introuvable: ${KDB}"
                continue
            fi

            KOUT="${KOUTROOT}/${RUN_LABEL}/${SAMPLE}"
            mkdir -p "${KOUT}"

            if [[ -s "${FASTP_MERGED}" ]]; then
                echo "  [Kraken2 ${KLABEL}] ${UNIT} (merged)"
                set +e
                kraken2 --conf 0.2 --db "${KDB}" --threads "${THREADS}" \
                    --output "${KOUT}/${UNIT}_merged.kraken" \
                    --report "${KOUT}/${UNIT}_merged.report" \
                    "${FASTP_MERGED}"
                set -e
            fi

            if [[ -s "${FASTP_UNMERGED_R1}" && -s "${FASTP_UNMERGED_R2}" ]]; then
                echo "  [Kraken2 ${KLABEL}] ${UNIT} (unmerged)"
                set +e
                kraken2 --conf 0.2 --paired --db "${KDB}" --threads "${THREADS}" \
                    --output "${KOUT}/${UNIT}_unmerged.kraken" \
                    --report "${KOUT}/${UNIT}_unmerged.report" \
                    "${FASTP_UNMERGED_R1}" "${FASTP_UNMERGED_R2}"
                set -e
            fi

        done

        echo "  >>> Pipeline complet termine pour ${UNIT} <<<"

    done
done

echo ""
echo "  STEP 1 done. Tous les runs/echantillons ont ete traites individuellement."

# ==============================================================================
# STEP 2 : Bracken — sur chaque report Kraken2 (k25/k29/k35/pracken),
#          merged/unmerged separement, run par run
# ==============================================================================

echo "======================================================================="
echo " STEP 2: Bracken par run x echantillon x base k-mer (incl. Pracken)"
echo "======================================================================="

run_bracken_for_kmer() {
    local KRAKEN_ROOT="$1"
    local DB="$2"
    local OUT_ROOT_K="$3"
    local KLABEL="$4"
    local READ_LEN="$5"
    local LEVEL="$6"
    local THRESHOLD="$7"

    if [[ ! -d "${DB}" ]]; then
        echo "  [SKIP Bracken ${KLABEL}] base introuvable: ${DB}"
        return
    fi

    echo "-----------------------------------------------------------------"
    echo " Bracken pour ${KLABEL} (read-len=${READ_LEN}, level=${LEVEL}, threshold=${THRESHOLD})"
    echo " Kraken root : ${KRAKEN_ROOT}"
    echo " Output root : ${OUT_ROOT_K}"
    echo "-----------------------------------------------------------------"

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
                echo "  [Bracken ${KLABEL}] ${UNIT} (merged)"
                set +e
                bracken \
                    -d "${DB}" \
                    -i "${MERGED_REPORT}" \
                    -o "${OUT_DIR}/${UNIT}_merged.bracken" \
                    -w "${OUT_DIR}/${UNIT}_merged.bracken.report" \
                    -r "${READ_LEN}" \
                    -l "${LEVEL}" \
                    -t "${THRESHOLD}"
                set -e
            else
                echo "  [SKIP ${KLABEL}] ${UNIT} (merged) : report absent -> ${MERGED_REPORT}"
            fi

            if [[ -f "${UNMERGED_REPORT}" ]]; then
                echo "  [Bracken ${KLABEL}] ${UNIT} (unmerged)"
                set +e
                bracken \
                    -d "${DB}" \
                    -i "${UNMERGED_REPORT}" \
                    -o "${OUT_DIR}/${UNIT}_unmerged.bracken" \
                    -w "${OUT_DIR}/${UNIT}_unmerged.bracken.report" \
                    -r "${READ_LEN}" \
                    -l "${LEVEL}" \
                    -t "${THRESHOLD}"
                set -e
            else
                echo "  [SKIP ${KLABEL}] ${UNIT} (unmerged) : report absent -> ${UNMERGED_REPORT}"
            fi
        done
    done
}

run_bracken_for_kmer "${KRAKEN_K25_DIR}" "${KRAKEN2_DB_K25}" "${BRACKEN_DIR}/k25" "k25" "${BRACKEN_READ_LEN}" "${BRACKEN_LEVEL}" "${BRACKEN_THRESHOLD}"
run_bracken_for_kmer "${KRAKEN_K29_DIR}" "${KRAKEN2_DB_K29}" "${BRACKEN_DIR}/k29" "k29" "${BRACKEN_READ_LEN}" "${BRACKEN_LEVEL}" "${BRACKEN_THRESHOLD}"
run_bracken_for_kmer "${KRAKEN_K35_DIR}" "${KRAKEN2_DB_K35}" "${BRACKEN_DIR}/k35" "k35" "${BRACKEN_READ_LEN}" "${BRACKEN_LEVEL}" "${BRACKEN_THRESHOLD}"

# Base Pracken : les fichiers database<N>mers.kmer_distrib sont deja fournis
# (50/75/100/150/200/250/300). On boucle sur PRACKEN_READ_LENS pour ne generer
# que les longueurs demandees, en verifiant que le fichier kmer_distrib
# correspondant existe bien avant d'appeler bracken.
for RLEN in "${PRACKEN_READ_LENS[@]}"; do
    KMER_DISTRIB_FILE="${KRAKEN2_DB_PRACKEN}/database${RLEN}mers.kmer_distrib"
    if [[ ! -f "${KMER_DISTRIB_FILE}" ]]; then
        echo "  [SKIP Bracken pracken] Aucun fichier kmer_distrib pour read-len=${RLEN} (${KMER_DISTRIB_FILE} absent)"
        continue
    fi
    run_bracken_for_kmer "${KRAKEN_PRACKEN_DIR}" "${KRAKEN2_DB_PRACKEN}" "${BRACKEN_DIR}/pracken" "pracken" "${RLEN}" "${PRACKEN_LEVEL}" "${PRACKEN_THRESHOLD}"
done

echo ""
echo "  STEP 2 done."

# ==============================================================================
# STEP 3 : Export MPA (Bracken + Kraken) par base k-mer, merged/unmerged
#          separement, run par run — inclut desormais k25/k29/k35/pracken
# ==============================================================================

echo "======================================================================="
echo " STEP 3: Export MPA (Bracken + Kraken) par base k-mer (incl. Pracken)"
echo "======================================================================="

gen_mpa_for_source() {
    local SOURCE="$1"        # "bracken" ou "kraken"
    local KLABEL="$2"        # k25/k29/k35/pracken
    local REPORT_ROOT="$3"   # racine des reports
    local MPA_OUT_ROOT="$4"  # racine des .mpa

    for i in "${!ALL_RUNS[@]}"; do
        RUN_LABEL="${RUN_LABELS[$i]}"

        for SAMPLE in "${SAMPLES[@]}"; do
            UNIT="${SAMPLE}_${RUN_LABEL}"

            if [[ "${SOURCE}" == "bracken" ]]; then
                REPORT_MERGED="${REPORT_ROOT}/${RUN_LABEL}/${SAMPLE}/${UNIT}_merged.bracken.report"
                REPORT_UNMERGED="${REPORT_ROOT}/${RUN_LABEL}/${SAMPLE}/${UNIT}_unmerged.bracken.report"
            else
                REPORT_MERGED="${REPORT_ROOT}/${RUN_LABEL}/${SAMPLE}/${UNIT}_merged.report"
                REPORT_UNMERGED="${REPORT_ROOT}/${RUN_LABEL}/${SAMPLE}/${UNIT}_unmerged.report"
            fi

            if [[ -f "${REPORT_MERGED}" ]]; then
                MPA_MERGED="${MPA_OUT_ROOT}/merged/${UNIT}_merged.${SOURCE}.mpa"
                echo "  [${SOURCE} ${KLABEL}] ${UNIT} (merged) -> ${MPA_MERGED}"
                python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" -r "${REPORT_MERGED}" -o "${MPA_MERGED}"
            fi

            if [[ -f "${REPORT_UNMERGED}" ]]; then
                MPA_UNMERGED="${MPA_OUT_ROOT}/unmerged/${UNIT}_unmerged.${SOURCE}.mpa"
                echo "  [${SOURCE} ${KLABEL}] ${UNIT} (unmerged) -> ${MPA_UNMERGED}"
                python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" -r "${REPORT_UNMERGED}" -o "${MPA_UNMERGED}"
            fi
        done
    done
}

# Bracken : reports dans 06_bracken/k*/run*/cop*/
gen_mpa_for_source "bracken" "k25"     "${BRACKEN_DIR}/k25"     "${MPA_DIR}/k25"
gen_mpa_for_source "bracken" "k29"     "${BRACKEN_DIR}/k29"     "${MPA_DIR}/k29"
gen_mpa_for_source "bracken" "k35"     "${BRACKEN_DIR}/k35"     "${MPA_DIR}/k35"
gen_mpa_for_source "bracken" "pracken" "${BRACKEN_DIR}/pracken" "${MPA_DIR}/pracken"

# Kraken : reports dans 05_kraken2_k*/run*/cop*/
gen_mpa_for_source "kraken" "k25"     "${KRAKEN_K25_DIR}"     "${MPA_DIR}/k25"
gen_mpa_for_source "kraken" "k29"     "${KRAKEN_K29_DIR}"     "${MPA_DIR}/k29"
gen_mpa_for_source "kraken" "k35"     "${KRAKEN_K35_DIR}"     "${MPA_DIR}/k35"
gen_mpa_for_source "kraken" "pracken" "${KRAKEN_PRACKEN_DIR}" "${MPA_DIR}/pracken"

echo ""
echo "  STEP 3 done."

# ==============================================================================
# STEP 4 : Combinaison intra-base k-mer (equivalent de l'ancien STEP 3 de
#          50_per_run_pipeline-3.sh) — une table combinee par base et par mode
# ==============================================================================

echo "======================================================================="
echo " STEP 4: Tables combinees intra-base (k25/k29/k35/pracken, merged/unmerged)"
echo "======================================================================="

combine_intra_kmer() {
    local KLABEL="$1"
    local MPA_ROOT="$2"

    for MODE in merged unmerged; do
        local SUBDIR="${MPA_ROOT}/${MODE}"
        [[ -d "${SUBDIR}" ]] || continue

        declare -a files=()
        while IFS= read -r f; do
            files+=("$f")
        done < <(find "${SUBDIR}" -maxdepth 1 -type f -name "*.bracken.mpa" | sort)

        if [[ ${#files[@]} -gt 0 ]]; then
            echo "  -> Combinaison table MPA ${MODE} (bracken) pour ${KLABEL}"
            python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" -i "${files[@]}" -o "${SUBDIR}/combined_per_run_${MODE}.tsv"
        fi
    done
}

combine_intra_kmer "k25"     "${MPA_DIR}/k25"
combine_intra_kmer "k29"     "${MPA_DIR}/k29"
combine_intra_kmer "k35"     "${MPA_DIR}/k35"
combine_intra_kmer "pracken" "${MPA_DIR}/pracken"

echo ""
echo "  STEP 4 done."

# ==============================================================================
# STEP 5 : Table finale globale — equivalent de 51_mpa_all-2.sh, mais integrant
#          desormais k25/k29/k35/pracken, avec header informatif
#          Ordre colonnes : k25/k29/k35/pracken -> merged/unmerged -> cop4* ->
#          run* -> bracken/kraken
# ==============================================================================

echo "======================================================================="
echo " STEP 5: Table finale globale (k25/k29/k35/pracken x merged/unmerged x bracken/kraken)"
echo "======================================================================="

FINAL_MPA="${MPA_DIR}/combined_all_k25_k29_k35_pracken_merged_unmerged_bracken_kraken.tsv"

declare -a MPA_FILES=()

for KLABEL in k25 k29 k35 pracken; do
    for MODE in merged unmerged; do
        MPA_SUBDIR="${MPA_DIR}/${KLABEL}/${MODE}"
        [[ -d "${MPA_SUBDIR}" ]] || continue

        # Bracken d'abord
        while IFS= read -r f; do
            MPA_FILES+=("$f")
        done < <(find "${MPA_SUBDIR}" -maxdepth 1 -type f -name "*.bracken.mpa" | sort)

        # Kraken ensuite
        while IFS= read -r f; do
            MPA_FILES+=("$f")
        done < <(find "${MPA_SUBDIR}" -maxdepth 1 -type f -name "*.kraken.mpa" | sort)
    done
done

if [[ ${#MPA_FILES[@]} -eq 0 ]]; then
    echo "[ERREUR] Aucun fichier .mpa genere dans ${MPA_DIR}/k*|pracken/merged|unmerged."
    echo "         La table finale combinee ne sera pas produite."
else
    echo "[INFO] Nombre de fichiers MPA detectes : ${#MPA_FILES[@]}"
    printf '%s\n' "${MPA_FILES[@]}"

    # Construction du header informatif
    declare -a HEADER_COLS=()

    for f in "${MPA_FILES[@]}"; do
        basename_f=$(basename "$f")

        klabel=$(echo "$f" | awk -F'/' '{for(i=1;i<=NF;i++) if($i=="k25"||$i=="k29"||$i=="k35"||$i=="pracken"){print $i; break}}')
        mode=$(echo "$f" | awk -F'/' '{for(i=1;i<=NF;i++) if($i=="merged"||$i=="unmerged"){print $i; break}}')
        sample=$(echo "$basename_f" | sed -E 's/^([^_]+)_.*$/\1/')
        run=$(echo "$basename_f" | sed -E 's/^[^_]+_([^_]+)_.*$/\1/')
        source=$(echo "$basename_f" | sed -E 's/^.*_([^_]+)\.mpa$/\1/')

        colname="${klabel}_${mode}_${sample}_${run}_${source}"
        HEADER_COLS+=("$colname")
    done

    echo "[INFO] Nombre de colonnes dans le header : ${#HEADER_COLS[@]}"

    echo "== Combinaison des fichiers MPA avec combine_mpa.py =="
    python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" \
        -i "${MPA_FILES[@]}" \
        -o "${FINAL_MPA}"

    if [[ ! -f "${FINAL_MPA}" ]]; then
        echo "[ERREUR] Le fichier final ${FINAL_MPA} n'a pas ete cree."
    else
        echo "== Reecriture du header de la table finale =="

        tmp_out="${FINAL_MPA}.tmp"

        {
            IFS=$'\n' read -r header_line
            IFS=$'\t' read -r -a header_fields <<< "$header_line"

            new_header=("${header_fields[0]}")
            for col in "${HEADER_COLS[@]}"; do
                new_header+=("$col")
            done

            printf '%s' "${new_header[0]}"
            for ((idx=1; idx<${#new_header[@]}; idx++)); do
                printf '\t%s' "${new_header[idx]}"
            done
            printf '\n'

            cat
        } < "${FINAL_MPA}" > "${tmp_out}"

        mv "${tmp_out}" "${FINAL_MPA}"

        echo "[INFO] Table finale ecrite dans : ${FINAL_MPA}"
        echo "[INFO] Ordre des colonnes : k25/k29/k35/pracken -> merged/unmerged -> cop4* -> run* -> bracken/kraken"
        echo "[INFO] Exemple de nom de colonne : pracken_merged_cop408_run1_bracken"
    fi
fi

# ==============================================================================
# RESUME FINAL
# ==============================================================================

echo ""
echo "======================================================================="
echo " ALL STEPS COMPLETED — Pipeline par run, avec k25/k29/k35 + Pracken"
echo "======================================================================="
echo " Racine des resultats  : ${OUT_ROOT}"
echo " BBDuk                 : ${BBDUK_DIR}/<run>/<sample>/"
echo " FastUniq               : ${FASTUNIQ_DIR}/<run>/<sample>/"
echo " Clumpify               : ${CLUMPIFY_DIR}/<run>/<sample>/"
echo " fastp                  : ${FASTP_DIR}/<run>/<sample>/"
echo " Kraken2 k25/k29/k35    : ${KRAKEN_K25_DIR}, ${KRAKEN_K29_DIR}, ${KRAKEN_K35_DIR}"
echo " Kraken2 Pracken        : ${KRAKEN_PRACKEN_DIR}"
echo " Bracken                : ${BRACKEN_DIR}/k25|k29|k35|pracken/<run>/<sample>/"
echo " Tables MPA par base    : ${MPA_DIR}/k25|k29|k35|pracken/merged|unmerged/combined_per_run_*.tsv"
echo " Table finale globale   : ${FINAL_MPA}"
echo ""
echo " Convention de nommage : <sample>_<run> puis suffixe merged/unmerged"
echo " Exemple : cop408_run1_merged.report / cop408_run3_unmerged.bracken.report"
echo "======================================================================="
