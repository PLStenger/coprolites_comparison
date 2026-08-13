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
# Script : 10_mapdamage_following.sh
# Author : Pierre-Louis Stenger
# Purpose :
#   Genere les profils MapDamage (dommages aDNA) pour chaque echantillon
#   "grouped" (cop408, cop410, cop412, cop414, cop417), a partir des reads
#   merged assignes par la cascade Kraken2 "following" (k35 -> k29 -> k25)
#   produite par 09_kraken_following_k35_29_25_mpa.sh.
#
#   Cible UNIQUEMENT 3 taxons :
#     - Ovis_aries       (mouton)   TaxID 9940
#     - Capra_hircus     (chevre)   TaxID 9925
#     - Corylus_avellana (noisetier) TaxID 13451
#
#   Source des reads :
#     - Fichier kraken combine (k35+k29+k25) :
#       09_kraken_following/combined_kraken_for_mapdamage/<sample>_following.kraken
#     - Fastq merged original :
#       06_fastp/clean_<sample>_grouped_dedup_fastp_merged.fastq.gz
#
#   Pour chaque (echantillon x taxon) :
#     1. Extraction des reads assignes au taxon via extract_kraken_reads.py
#        (--include-children pour recuperer aussi les assignations a des
#        rangs infra-specifiques)
#     2. Mapping BWA aln + samse (single-end, reads merged)
#     3. Calcul du taux de mapping (samtools)
#     4. MapDamage (graphiques de dommages, --no-stats)
#
#   Arborescence de sortie :
#     10_mapdamage_following/<Taxon>/<sample>/
#       <sample>_<Taxon>_following.fastq        (temporaire, supprime apres usage)
#       <sample>_<Taxon>_following.sorted.bam(.bai)
#       <sample>_<Taxon>_following_mapDamage/    (graphiques MapDamage)
#     10_summary_tables/mapping_bwa_info_following.tsv (taux de mapping global)
# ==============================================================================

echo ""
echo "=========================================="
echo "MapDamage sur le pipeline kraken_following (k35+k29+k25, merged only)"
echo "3 taxons cibles : mouton, chevre, noisetier"
echo "=========================================="
echo ""

# ==============================================================================
# INITIALISATION CONDA
# ==============================================================================

module load conda/4.12.0
source ~/.bashrc
conda activate mapdamage_py39
echo "Environnement mapdamage_py39 active"

# ==============================================================================
# CONFIGURATION GLOBALE
# ==============================================================================

WORKDIR="/home/plstenge/coprolites_comparison"

FASTP_DIR="${WORKDIR}/06_fastp"
FOLLOWING_DIR="${WORKDIR}/09_kraken_following"
COMBINED_KRAKEN_DIR="${FOLLOWING_DIR}/combined_kraken_for_mapdamage"

KRAKENTOOLS_DIR="${WORKDIR}/08_krakentools/KrakenTools"

THREADS="${SLURM_CPUS_PER_TASK:-36}"

SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")

MAPDAMAGE_OUT_ROOT="${WORKDIR}/10_mapdamage_following"
mkdir -p "${MAPDAMAGE_OUT_ROOT}"

LOGFILE="${WORKDIR}/00_scripts/mapdamage_following_$(date +%Y%m%d_%H%M%S).txt"
touch "${LOGFILE}"

SUMMARY_DIR="${WORKDIR}/10_summary_tables"
mkdir -p "${SUMMARY_DIR}"
MAPPING_INFO="${SUMMARY_DIR}/mapping_bwa_info_following.tsv"
echo -e "Sample\tTaxon\tTotal_Reads\tMapped_Reads\tMapping_Rate(%)" > "${MAPPING_INFO}"

echo "$(date): Script MapDamage following demarre" | tee -a "${LOGFILE}"

# ==============================================================================
# GENOMES DE REFERENCE — uniquement les 3 taxons demandes
# ==============================================================================
# Format : ["Nom_espece"]="TaxID:/chemin/vers/genome.fa"

declare -A TAXONS=(
    ["Ovis_aries"]="9940:/home/plstenge/genomes/Ovis_aries/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa"
    ["Capra_hircus"]="9925:/home/plstenge/genomes/Capra_hircus/Capra_hircus.ARS1.dna.toplevel.fa"
    ["Corylus_avellana"]="13451:/home/plstenge/genomes/Corylus_avellana/Corylus_avellana_CavTom2PMs_1_0.fasta"
)

echo ""
echo "Genomes de reference configures: ${#TAXONS[@]}"
for sp in "${!TAXONS[@]}"; do
    echo " - ${sp}"
done
echo ""

# ==============================================================================
# FONCTIONS
# ==============================================================================

check_bwa_index_complete() {
    local ref="$1"
    [[ -f "${ref}.bwt" && -f "${ref}.amb" && -f "${ref}.ann" && -f "${ref}.pac" && -f "${ref}.sa" ]]
}

calculate_mapping_rate() {
    local bam_file="$1"
    local sample_name="$2"
    local taxon="$3"
    local outfile="$4"

    if [[ -f "${bam_file}" ]]; then
        local total mapped rate
        total=$(samtools view -c "${bam_file}" 2>/dev/null || echo 0)
        mapped=$(samtools view -c -F 4 "${bam_file}" 2>/dev/null || echo 0)
        rate=0
        if [[ "${total}" -gt 0 ]]; then
            rate=$(echo "scale=1; ${mapped} * 100 / ${total}" | bc 2>/dev/null || echo "0")
        fi
        echo -e "${sample_name}\t${taxon}\t${total}\t${mapped}\t${rate}" >> "${outfile}"
        echo "    Stats: ${mapped}/${total} mappes (${rate}%)" | tee -a "${LOGFILE}"
    fi
}

run_mapdamage() {
    local bam_file="$1"
    local ref_fasta="$2"
    local output_dir="$3"

    if [[ ! -s "${bam_file}" ]]; then
        echo "    (BAM vide ou inexistant, MapDamage ignore)"
        return 0
    fi

    mkdir -p "${output_dir}"
    echo -n "    -> MapDamage en cours..."

    if mapDamage -i "${bam_file}" -r "${ref_fasta}" --folder "${output_dir}" --no-stats 2>>"${LOGFILE}"; then
        echo " OK (graphiques generes)"
    else
        echo " echoue (non bloquant)"
    fi
}

# ==============================================================================
# PHASE 1 : INDEXATION DES GENOMES DE REFERENCE
# ==============================================================================

echo ""
echo "=========================================="
echo "PHASE 1: Indexation des genomes de reference (3 taxons)"
echo "=========================================="
echo ""

declare -A VALID_GENOMES
INDEXATION_ERRORS=0

for GROUP in "${!TAXONS[@]}"; do
    entry="${TAXONS[${GROUP}]}"
    ORIGINAL_REF="${entry#*:}"

    echo -n "-> ${GROUP}: "

    if [[ ! -f "${ORIGINAL_REF}" ]]; then
        echo "GENOME MANQUANT: ${ORIGINAL_REF}"
        ((INDEXATION_ERRORS++))
        continue
    fi

    echo -n "($(du -h "${ORIGINAL_REF}" | cut -f1)) ... "

    if check_bwa_index_complete "${ORIGINAL_REF}"; then
        echo -n "[BWA:OK] "
    else
        echo -n "[BWA: indexation...] "
        if timeout 3600 bwa index "${ORIGINAL_REF}" 2>>"${LOGFILE}"; then
            check_bwa_index_complete "${ORIGINAL_REF}" || { echo "ERREUR: index BWA incomplet"; ((INDEXATION_ERRORS++)); continue; }
            echo -n "[OK] "
        else
            echo "ERREUR: bwa index timeout/echec"
            ((INDEXATION_ERRORS++))
            continue
        fi
    fi

    if [[ ! -f "${ORIGINAL_REF}.fai" ]]; then
        if samtools faidx "${ORIGINAL_REF}" 2>>"${LOGFILE}"; then
            echo "[SAM:OK]"
        else
            echo "ERREUR: samtools faidx echoue"
            ((INDEXATION_ERRORS++))
            continue
        fi
    else
        echo "[SAM:OK]"
    fi

    VALID_GENOMES[${GROUP}]="${ORIGINAL_REF}"
done

echo ""
[[ ${INDEXATION_ERRORS} -gt 0 ]] && echo "ATTENTION: ${INDEXATION_ERRORS} genome(s) ont echoue l'indexation"
echo "Genomes prets pour le mapping: ${#VALID_GENOMES[@]} / ${#TAXONS[@]}"

if [[ ${#VALID_GENOMES[@]} -eq 0 ]]; then
    echo "ERREUR CRITIQUE: Aucun genome valide disponible!"
    exit 1
fi

# ==============================================================================
# PHASE 2 : EXTRACTION + MAPPING + MAPDAMAGE PAR SAMPLE x TAXON
# ==============================================================================

echo ""
echo "=========================================="
echo "PHASE 2: Extraction reads (following k35+k29+k25) + mapping BWA + MapDamage"
echo "=========================================="
echo ""

for SAMPLE in "${SAMPLES[@]}"; do

    COMBINED_KRAKEN="${COMBINED_KRAKEN_DIR}/${SAMPLE}_following.kraken"
    MERGED_FASTQ="${FASTP_DIR}/clean_${SAMPLE}_grouped_dedup_fastp_merged.fastq.gz"

    if [[ ! -s "${COMBINED_KRAKEN}" ]]; then
        echo "ATTENTION: fichier kraken combine introuvable pour ${SAMPLE} (${COMBINED_KRAKEN}), on saute."
        continue
    fi

    if [[ ! -s "${MERGED_FASTQ}" ]]; then
        echo "ATTENTION: fastq merged introuvable pour ${SAMPLE} (${MERGED_FASTQ}), on saute."
        continue
    fi

    echo ""
    echo "-----------------------------------------------------------------"
    echo " Echantillon : ${SAMPLE} (grouped, merged, following k35+k29+k25)"
    echo "-----------------------------------------------------------------"

    for GROUP in "${!VALID_GENOMES[@]}"; do
        REF="${VALID_GENOMES[${GROUP}]}"
        TAX_ID="${TAXONS[${GROUP}]%%:*}"

        OUTDIR="${MAPDAMAGE_OUT_ROOT}/${GROUP}/${SAMPLE}"
        mkdir -p "${OUTDIR}"

        echo " -> ${GROUP} (TaxID: ${TAX_ID})"

        OUT_FQ="${OUTDIR}/${SAMPLE}_${GROUP}_following.fastq"

        set +e
        python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
            -k "${COMBINED_KRAKEN}" \
            -s "${MERGED_FASTQ}" \
            -t "${TAX_ID}" \
            --include-children \
            -o "${OUT_FQ}" \
            --fastq-output 2>>"${LOGFILE}"
        set -e

        if [[ ! -s "${OUT_FQ}" ]]; then
            echo "    (aucun read extrait, on saute)"
            rm -f "${OUT_FQ}" 2>/dev/null
            continue
        fi

        READ_COUNT=$(grep -c "^@" "${OUT_FQ}" 2>/dev/null || echo 0)
        echo "    ${READ_COUNT} reads extraits"

        SAI="${OUTDIR}/${SAMPLE}_${GROUP}_following.sai"
        BAM="${OUTDIR}/${SAMPLE}_${GROUP}_following.sorted.bam"

        set +e
        bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "${REF}" "${OUT_FQ}" > "${SAI}" 2>>"${LOGFILE}"
        bwa samse "${REF}" "${SAI}" "${OUT_FQ}" 2>>"${LOGFILE}" | \
            samtools view -bS - 2>>"${LOGFILE}" | \
            samtools sort -o "${BAM}" - 2>>"${LOGFILE}"
        set -e

        [[ -f "${BAM}" ]] || { echo "    ERREUR mapping"; rm -f "${SAI}" "${OUT_FQ}"; continue; }

        samtools index "${BAM}" 2>>"${LOGFILE}"

        calculate_mapping_rate "${BAM}" "${SAMPLE}" "${GROUP}" "${MAPPING_INFO}"

        run_mapdamage "${BAM}" "${REF}" "${OUTDIR}/${SAMPLE}_${GROUP}_following_mapDamage"

        rm -f "${SAI}" "${OUT_FQ}"
    done
done

# ==============================================================================
# RESUME ET FIN
# ==============================================================================

echo ""
echo "=========================================="
echo "PIPELINE MAPDAMAGE FOLLOWING TERMINE"
echo "Date: $(date)"
echo "=========================================="
echo ""
echo "Statistiques de mapping : ${MAPPING_INFO}"
echo "Log complet : ${LOGFILE}"
echo "Resultats BAM+MapDamage : ${MAPDAMAGE_OUT_ROOT}/<Taxon>/<sample>/"
echo ""

if [[ -f "${MAPPING_INFO}" ]]; then
    echo "Apercu des resultats :"
    column -t "${MAPPING_INFO}" | head -50
fi

conda deactivate

echo "Script termine."
echo "Pipeline MapDamage following termine le $(date +%Y-%m-%d\ %H:%M:%S). Resultats: ${MAPDAMAGE_OUT_ROOT}" | \
    mail -s "Pipeline Coprolites MapDamage following - Termine" pierrelouis.stenger@gmail.com 2>/dev/null || true

exit 0
