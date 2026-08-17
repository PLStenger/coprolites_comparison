#!/bin/bash

#SBATCH --job-name=10_mapdamage_following
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=36
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/10_mapdamage_following.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/10_mapdamage_following.out"

# ==============================================================================
# MapDamage sur la cascade Kraken2 k35 -> k29 -> k25, reads merged uniquement, only on "goat" and "sheep" sample
#
# IMPORTANT
# Les .kraken k29 et k25 n'ont pas le meme nombre ni le meme ordre de reads que
# le FASTQ merged original : ils correspondent respectivement aux reads non
# classes par k35, puis par k35+k29. On reconstruit donc ces deux FASTQ
# intermediaires et on extrait les reads de chaque taxon dans chaque niveau
# separement, avec SON report Kraken2. Les trois FASTQ extraits sont ensuite
# concatenes avant le mapping et MapDamage.
# ==============================================================================

module load conda/4.12.0
source ~/.bashrc
conda activate mapdamage_py39

WORKDIR="/home/plstenge/coprolites_comparison"
FASTP_DIR="${WORKDIR}/06_fastp"
FOLLOWING_DIR="${WORKDIR}/09_kraken_following"
K35_DIR="${FOLLOWING_DIR}/k35"
K29_DIR="${FOLLOWING_DIR}/k29"
K25_DIR="${FOLLOWING_DIR}/k25"
KRAKENTOOLS_DIR="${WORKDIR}/08_krakentools/KrakenTools"

THREADS="${SLURM_CPUS_PER_TASK:-36}"
SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")

OUT_ROOT="${WORKDIR}/10_mapdamage_following"
TMP_ROOT="${OUT_ROOT}/tmp_cascade_fastq"
SUMMARY_DIR="${WORKDIR}/10_summary_tables"
mkdir -p "${OUT_ROOT}" "${TMP_ROOT}" "${SUMMARY_DIR}"

LOGFILE="${WORKDIR}/00_scripts/mapdamage_following_$(date +%Y%m%d_%H%M%S).txt"
MAPPING_INFO="${SUMMARY_DIR}/mapping_bwa_info_following.tsv"
echo -e "Sample\tTaxon\tExtracted_k35\tExtracted_k29\tExtracted_k25\tTotal_extracted\tMapped_reads\tMapping_rate_percent" > "${MAPPING_INFO}"

bwa index /home/plstenge/genomes/Alces_alces/GCA_059051365.1_mAlcAlc2_p1.1_genomic.fna
bwa index /home/plstenge/genomes/Rangifer_tarandus/GCA_949782905.1_mRanTar1.h1.1_genomic.fna

declare -A TAXONS=(
   # ["Ovis_aries"]="9940:/home/plstenge/genomes/Ovis_aries/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa"
   # ["Capra_hircus"]="9925:/home/plstenge/genomes/Capra_hircus/Capra_hircus.ARS1.dna.toplevel.fa"
   # ["Corylus_avellana"]="13451:/home/plstenge/genomes/Corylus_avellana/Corylus_avellana_CavTom2PMs_1_0.fasta"
   ["Mus_musculus"]="10090:/home/plstenge/genomes/Mus_musculus/Mus_musculus.GRCm39.dna.toplevel.fa"
   ["Alces_alces"]="9852:/home/plstenge/genomes/Alces_alces/GCA_059051365.1_mAlcAlc2_p1.1_genomic.fna"
   ["Rangifer_tarandus"]="9870:/home/plstenge/genomes/Rangifer_tarandus/GCA_949782905.1_mRanTar1.h1.1_genomic.fna"
   ["Cannabis_sativa"]="3483:/home/plstenge/genomes/Cannabis_sativa/GCF_029168945.1_ASM2916894v1_genomic.fna"
   ["Homo_sapiens"]="9606:/home/plstenge/genomes/Homo_sapiens/GCF_000001405.40_GRCh38.p14_genomic.fna"
   ["Rubus_caesius"]="75065:/home/plstenge/genomes/Rubus_caesius/GCA_964235055.1_drRubCaes1.hap1.1_genomic.fna"
)

count_fastq_reads() {
    [[ -s "$1" ]] && grep -c '^@' "$1" || echo 0
}

check_bwa_index() {
    local ref="$1"
    [[ -f "${ref}.amb" && -f "${ref}.ann" && -f "${ref}.bwt" && -f "${ref}.pac" && -f "${ref}.sa" ]]
}

extract_unclassified() {
    # TaxID 0 correspond exactement aux lignes U : pas de --include-children,
    # et donc pas besoin de --report.
    local kraken="$1" fastq="$2" out_fastq="$3"
    python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
        -k "${kraken}" -s "${fastq}" -t 0 -o "${out_fastq}" --fastq-output
}

extract_taxon() {
    # --include-children exige obligatoirement le report correspondant au
    # fichier .kraken. Cela recupere le taxon et ses descendants.
    local kraken="$1" report="$2" fastq="$3" taxid="$4" out_fastq="$5"
    python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
        -k "${kraken}" -r "${report}" -s "${fastq}" -t "${taxid}" \
        --include-children -o "${out_fastq}" --fastq-output
}

# Verification et indexation des references.
declare -A VALID_REF
for TAXON in "${!TAXONS[@]}"; do
    REF="${TAXONS[${TAXON}]#*:}"
    if [[ ! -f "${REF}" ]]; then
        echo "ERREUR: genome absent pour ${TAXON}: ${REF}" | tee -a "${LOGFILE}"
        continue
    fi
    if ! check_bwa_index "${REF}"; then
        echo "Indexation BWA: ${TAXON}" | tee -a "${LOGFILE}"
        bwa index "${REF}" 2>>"${LOGFILE}" || continue
    fi
    [[ -f "${REF}.fai" ]] || samtools faidx "${REF}" 2>>"${LOGFILE}" || continue
    VALID_REF[${TAXON}]="${REF}"
done

[[ ${#VALID_REF[@]} -gt 0 ]] || { echo "ERREUR: aucune reference indexable."; exit 1; }

echo "=========================================="
echo "PHASE 2: extraction cascade + mapping BWA + MapDamage"
echo "=========================================="

for SAMPLE in "${SAMPLES[@]}"; do
    ORIGINAL_FASTQ="${FASTP_DIR}/clean_${SAMPLE}_grouped_dedup_fastp_merged.fastq.gz"
    K35_KRAKEN="${K35_DIR}/${SAMPLE}_merged.kraken"
    K35_REPORT="${K35_DIR}/${SAMPLE}_merged.report"
    K29_KRAKEN="${K29_DIR}/${SAMPLE}_merged_unassigned_k35.kraken"
    K29_REPORT="${K29_DIR}/${SAMPLE}_merged_unassigned_k35.report"
    K25_KRAKEN="${K25_DIR}/${SAMPLE}_merged_unassigned_k29.kraken"
    K25_REPORT="${K25_DIR}/${SAMPLE}_merged_unassigned_k29.report"

    if [[ ! -s "${ORIGINAL_FASTQ}" || ! -s "${K35_KRAKEN}" || ! -s "${K35_REPORT}" ]]; then
        echo "ATTENTION: entree k35/FASTQ manquante pour ${SAMPLE}; echantillon ignore." | tee -a "${LOGFILE}"
        continue
    fi

    SAMPLE_TMP="${TMP_ROOT}/${SAMPLE}"
    mkdir -p "${SAMPLE_TMP}"
    UN_K35="${SAMPLE_TMP}/${SAMPLE}_unassigned_after_k35.fastq"
    UN_K29="${SAMPLE_TMP}/${SAMPLE}_unassigned_after_k29.fastq"

    echo "--- ${SAMPLE}: reconstruction des FASTQ de la cascade ---"
    extract_unclassified "${K35_KRAKEN}" "${ORIGINAL_FASTQ}" "${UN_K35}" 2>>"${LOGFILE}"
    if [[ -s "${K29_KRAKEN}" && -s "${UN_K35}" ]]; then
        extract_unclassified "${K29_KRAKEN}" "${UN_K35}" "${UN_K29}" 2>>"${LOGFILE}"
    else
        : > "${UN_K29}"
    fi

    for TAXON in "${!VALID_REF[@]}"; do
        TAXID="${TAXONS[${TAXON}]%%:*}"
        REF="${VALID_REF[${TAXON}]}"
        OUTDIR="${OUT_ROOT}/${TAXON}/${SAMPLE}"
        mkdir -p "${OUTDIR}"

        FQ_K35="${OUTDIR}/${SAMPLE}_${TAXON}_k35.fastq"
        FQ_K29="${OUTDIR}/${SAMPLE}_${TAXON}_k29.fastq"
        FQ_K25="${OUTDIR}/${SAMPLE}_${TAXON}_k25.fastq"
        FQ_ALL="${OUTDIR}/${SAMPLE}_${TAXON}_following.fastq"

        echo "  -> ${TAXON} (TaxID ${TAXID})"
        extract_taxon "${K35_KRAKEN}" "${K35_REPORT}" "${ORIGINAL_FASTQ}" "${TAXID}" "${FQ_K35}" 2>>"${LOGFILE}"

        if [[ -s "${K29_KRAKEN}" && -s "${K29_REPORT}" && -s "${UN_K35}" ]]; then
            extract_taxon "${K29_KRAKEN}" "${K29_REPORT}" "${UN_K35}" "${TAXID}" "${FQ_K29}" 2>>"${LOGFILE}"
        else
            : > "${FQ_K29}"
        fi

        if [[ -s "${K25_KRAKEN}" && -s "${K25_REPORT}" && -s "${UN_K29}" ]]; then
            extract_taxon "${K25_KRAKEN}" "${K25_REPORT}" "${UN_K29}" "${TAXID}" "${FQ_K25}" 2>>"${LOGFILE}"
        else
            : > "${FQ_K25}"
        fi

        N35=$(count_fastq_reads "${FQ_K35}")
        N29=$(count_fastq_reads "${FQ_K29}")
        N25=$(count_fastq_reads "${FQ_K25}")
        cat "${FQ_K35}" "${FQ_K29}" "${FQ_K25}" > "${FQ_ALL}"
        NTOT=$(count_fastq_reads "${FQ_ALL}")
        echo "     reads extraits : k35=${N35}; k29=${N29}; k25=${N25}; total=${NTOT}"

        if [[ "${NTOT}" -eq 0 ]]; then
            echo -e "${SAMPLE}\t${TAXON}\t${N35}\t${N29}\t${N25}\t0\t0\t0" >> "${MAPPING_INFO}"
            rm -f "${FQ_K35}" "${FQ_K29}" "${FQ_K25}" "${FQ_ALL}"
            continue
        fi

        SAI="${OUTDIR}/${SAMPLE}_${TAXON}_following.sai"
        BAM="${OUTDIR}/${SAMPLE}_${TAXON}_following.sorted.bam"
        bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "${REF}" "${FQ_ALL}" > "${SAI}" 2>>"${LOGFILE}"
        bwa samse "${REF}" "${SAI}" "${FQ_ALL}" 2>>"${LOGFILE}" | samtools view -bS - | samtools sort -o "${BAM}" -
        samtools index "${BAM}" 2>>"${LOGFILE}"

        TOTAL=$(samtools view -c "${BAM}" 2>/dev/null || echo 0)
        MAPPED=$(samtools view -c -F 4 "${BAM}" 2>/dev/null || echo 0)
        RATE=$(awk -v m="${MAPPED}" -v t="${TOTAL}" 'BEGIN {if (t>0) printf "%.1f", 100*m/t; else print 0}')
        echo -e "${SAMPLE}\t${TAXON}\t${N35}\t${N29}\t${N25}\t${NTOT}\t${MAPPED}\t${RATE}" >> "${MAPPING_INFO}"

        MD_DIR="${OUTDIR}/${SAMPLE}_${TAXON}_following_mapDamage"
        mapDamage -i "${BAM}" -r "${REF}" --folder "${MD_DIR}" --no-stats 2>>"${LOGFILE}" || echo "MapDamage echoue pour ${SAMPLE}/${TAXON}; voir ${LOGFILE}" >&2

        # Les BAM/.bai et resultats MapDamage sont conserves; seuls les FASTQ/SAI
        # d'extraction temporaires sont supprimes.
        rm -f "${SAI}" "${FQ_K35}" "${FQ_K29}" "${FQ_K25}" "${FQ_ALL}"
    done

    rm -rf "${SAMPLE_TMP}"
done

echo "=========================================="
echo "PIPELINE MAPDAMAGE FOLLOWING TERMINE"
echo "Table de mapping : ${MAPPING_INFO}"
echo "Log : ${LOGFILE}"
echo "=========================================="
conda deactivate
