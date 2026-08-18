#!/usr/bin/env bash
#SBATCH --job-name=13_extract_rangifer_kraken
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=/home/plstenge/coprolites_comparison/00_scripts/13_extract_rangifer_kraken_reads.err
#SBATCH --output=/home/plstenge/coprolites_comparison/00_scripts/13_extract_rangifer_kraken_reads.out

# ==============================================================================
# Extrait les sequences des reads attribues par Kraken2 au genre Rangifer
# (TaxID 9870), en incluant ses descendants, pour chaque echantillon Cop*.
#
# Sortie principale : FASTA avec les sequences ATGC... :
#   /home/plstenge/coprolites_comparison/13_extract_rangifer_kraken_reads/
#       <copXXX>/<copXXX>_Rangifer_kraken_reads.fasta.gz
#
# Le script traite la cascade Kraken k35 -> k29 -> k25 de la meme maniere que
# 10_mapdamage_following.sh : les FASTQ non assignes sont reconstruits entre
# les niveaux afin que chaque fichier .kraken corresponde au bon FASTQ.
# ==============================================================================

module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics

WORKDIR="/home/plstenge/coprolites_comparison"
FASTP_DIR="${WORKDIR}/06_fastp"
FOLLOWING_DIR="${WORKDIR}/09_kraken_following"
K35_DIR="${FOLLOWING_DIR}/k35"
K29_DIR="${FOLLOWING_DIR}/k29"
K25_DIR="${FOLLOWING_DIR}/k25"
KRAKENTOOLS_DIR="${WORKDIR}/08_krakentools/KrakenTools"

# TaxID du GENRE Rangifer. Avec --include-children : toutes les especes
# descendants de Rangifer dans le report Kraken sont incluses.
RANGIFER_TAXID="9870"

SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")
OUT_ROOT="${WORKDIR}/13_extract_rangifer_kraken_reads"
TMP_ROOT="${OUT_ROOT}/tmp_cascade_fastq"
SUMMARY_DIR="${OUT_ROOT}/summary"
mkdir -p "${OUT_ROOT}" "${TMP_ROOT}" "${SUMMARY_DIR}"

LOGFILE="${WORKDIR}/00_scripts/extract_rangifer_kraken_reads_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOGFILE}") 2>&1

for CMD in python3 seqkit; do
    command -v "${CMD}" >/dev/null || {
        echo "ERREUR: ${CMD} introuvable dans l'environnement metagenomics." >&2
        exit 1
    }
done

[[ -f "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" ]] || {
    echo "ERREUR: script KrakenTools absent : ${KRAKENTOOLS_DIR}/extract_kraken_reads.py" >&2
    exit 1
}

count_fastq_reads() {
    [[ -s "$1" ]] && awk 'END { print int(NR / 4) }' "$1" || echo 0
}

extract_unclassified() {
    local kraken="$1" fastq="$2" output_fastq="$3"
    python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
        -k "${kraken}" -s "${fastq}" -t 0 -o "${output_fastq}" --fastq-output
}

extract_rangifer() {
    local kraken="$1" report="$2" fastq="$3" output_fastq="$4"
    python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
        -k "${kraken}" -r "${report}" -s "${fastq}" -t "${RANGIFER_TAXID}" \
        --include-children -o "${output_fastq}" --fastq-output
}

printf 'sample\treads_k35\treads_k29\treads_k25\treads_total\n' \
    > "${SUMMARY_DIR}/rangifer_kraken_extraction_counts.tsv"

for SAMPLE in "${SAMPLES[@]}"; do
    ORIGINAL_FASTQ="${FASTP_DIR}/clean_${SAMPLE}_grouped_dedup_fastp_merged.fastq.gz"
    K35_KRAKEN="${K35_DIR}/${SAMPLE}_merged.kraken"
    K35_REPORT="${K35_DIR}/${SAMPLE}_merged.report"
    K29_KRAKEN="${K29_DIR}/${SAMPLE}_merged_unassigned_k35.kraken"
    K29_REPORT="${K29_DIR}/${SAMPLE}_merged_unassigned_k35.report"
    K25_KRAKEN="${K25_DIR}/${SAMPLE}_merged_unassigned_k29.kraken"
    K25_REPORT="${K25_DIR}/${SAMPLE}_merged_unassigned_k29.report"

    if [[ ! -s "${ORIGINAL_FASTQ}" || ! -s "${K35_KRAKEN}" || ! -s "${K35_REPORT}" ]]; then
        echo "ATTENTION: FASTQ ou resultats k35 manquants pour ${SAMPLE}; ignore."
        continue
    fi

    echo "=============================================================="
    echo "Extraction Rangifer : ${SAMPLE}"
    echo "=============================================================="

    SAMPLE_TMP="${TMP_ROOT}/${SAMPLE}"
    OUTDIR="${OUT_ROOT}/${SAMPLE}"
    mkdir -p "${SAMPLE_TMP}" "${OUTDIR}"

    UN_K35="${SAMPLE_TMP}/${SAMPLE}_unassigned_after_k35.fastq"
    UN_K29="${SAMPLE_TMP}/${SAMPLE}_unassigned_after_k29.fastq"
    FQ_K35="${OUTDIR}/${SAMPLE}_Rangifer_k35.fastq"
    FQ_K29="${OUTDIR}/${SAMPLE}_Rangifer_k29.fastq"
    FQ_K25="${OUTDIR}/${SAMPLE}_Rangifer_k25.fastq"
    FQ_ALL="${OUTDIR}/${SAMPLE}_Rangifer_kraken_reads.fastq"
    FASTA_ALL="${OUTDIR}/${SAMPLE}_Rangifer_kraken_reads.fasta"

    # k29 est aligne sur les reads U de k35; k25 sur les reads U de k29.
    extract_unclassified "${K35_KRAKEN}" "${ORIGINAL_FASTQ}" "${UN_K35}"
    if [[ -s "${K29_KRAKEN}" && -s "${UN_K35}" ]]; then
        extract_unclassified "${K29_KRAKEN}" "${UN_K35}" "${UN_K29}"
    else
        : > "${UN_K29}"
    fi

    # Rangifer a chaque niveau de la cascade, genre + descendants.
    extract_rangifer "${K35_KRAKEN}" "${K35_REPORT}" "${ORIGINAL_FASTQ}" "${FQ_K35}"
    if [[ -s "${K29_KRAKEN}" && -s "${K29_REPORT}" && -s "${UN_K35}" ]]; then
        extract_rangifer "${K29_KRAKEN}" "${K29_REPORT}" "${UN_K35}" "${FQ_K29}"
    else
        : > "${FQ_K29}"
    fi
    if [[ -s "${K25_KRAKEN}" && -s "${K25_REPORT}" && -s "${UN_K29}" ]]; then
        extract_rangifer "${K25_KRAKEN}" "${K25_REPORT}" "${UN_K29}" "${FQ_K25}"
    else
        : > "${FQ_K25}"
    fi

    # FASTQ conserve avec qualites; FASTA contient uniquement les identifiants
    # et les sequences ATGC..., sans qualites, et constitue la sortie demandee.
    cat "${FQ_K35}" "${FQ_K29}" "${FQ_K25}" > "${FQ_ALL}"
    seqkit seq -w 0 "${FQ_ALL}" > "${FASTA_ALL}"
    gzip -f "${FASTA_ALL}"

    N35=$(count_fastq_reads "${FQ_K35}")
    N29=$(count_fastq_reads "${FQ_K29}")
    N25=$(count_fastq_reads "${FQ_K25}")
    NTOTAL=$(count_fastq_reads "${FQ_ALL}")
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "${SAMPLE}" "${N35}" "${N29}" "${N25}" "${NTOTAL}" \
        >> "${SUMMARY_DIR}/rangifer_kraken_extraction_counts.tsv"

    echo "${SAMPLE}: ${NTOTAL} reads Rangifer ecrits dans ${FASTA_ALL}.gz"

    # Les FASTQ individuels et concatene sont conserves. Le FASTQ temporaire
    # de reconstruction de la cascade est supprime.
    rm -rf "${SAMPLE_TMP}"
done

echo "=============================================================="
echo "EXTRACTION TERMINEE"
echo "FASTA : ${OUT_ROOT}/<copXXX>/<copXXX>_Rangifer_kraken_reads.fasta.gz"
echo "Table : ${SUMMARY_DIR}/rangifer_kraken_extraction_counts.tsv"
echo "Log   : ${LOGFILE}"
echo "=============================================================="

conda deactivate
