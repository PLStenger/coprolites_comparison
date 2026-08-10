#!/bin/bash

#SBATCH --job-name=52_mapdamage
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=36
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/52_mapdamage.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/52_mapdamage.out"

# ==============================================================================
# Script : 52_mapdamage.sh
# Author : Pierre-Louis Stenger
# Purpose :
#   Genere les profils MapDamage (graphiques de dommages aDNA + taux de mapping)
#   pour chaque echantillon (cop408, cop410, cop412, cop414, cop417) x chaque
#   run (run1..run5), a partir des reads Kraken2/Bracken/Pracken extraits pour
#   5 taxons cibles :
#     - Ovis_aries        (mouton)
#     - Capra_hircus      (chevre)
#     - Olea_europaea     (olivier)
#     - Fraxinus_excelsior (frene)
#     - Corylus_avellana  (noisetier)
#
#   Sources de reads couvertes (identiques a la structure de
#   50_per_run_pipeline_fixed.sh) :
#     - Kraken2 k25    : 05_kraken2_k25/<run>/<sample>/<unit>_merged.kraken + _unmerged.kraken
#     - Kraken2 k29    : 05_kraken2_k29/<run>/<sample>/...
#     - Kraken2 k35    : 05_kraken2_k35/<run>/<sample>/...
#     - Kraken2 pracken: 05_kraken2_pracken/<run>/<sample>/...
#     - Bracken (k35 uniquement, comme demande) : reclassification a partir
#       des reads Kraken2 k35 filtres par les taxa retenus par Bracken
#       (extract_kraken_reads.py utilise directement le fichier .kraken k35,
#       car Bracken ne fait que reestimer l'abondance sans regenerer de reads;
#       le sous-ensemble "bracken_k35" ici designe donc les memes reads bruts
#       Kraken2 k35, mais on ne les traite QUE si un report Bracken existe
#       pour ce taxon/unit, garantissant la coherence avec les abondances
#       Bracken validees)
#
#   Pour chaque combinaison (sample x run x base x merge_status x taxon) :
#     1. Extraction des reads assignes au taxon via extract_kraken_reads.py
#     2. Mapping BWA (aln + samse, single-end car merged ET unmerged sont
#        traites comme des lots de reads independants après fastp merge/R1+R2)
#     3. Calcul du taux de mapping (samtools)
#     4. MapDamage (graphiques de dommages, --no-stats pour eviter les soucis
#        Rcpp/GCC, comme dans l'ancien script)
#
#   Arborescence de sortie :
#     12_mapdamage_per_run/
#       <base>/<merge_status>/<taxon>/<run>/<sample>/
#         <unit>_<taxon>_<merge_status>.fastq (temporaire, supprime apres usage)
#         <unit>_<taxon>_<merge_status>.sorted.bam(.bai)
#         <unit>_<taxon>_<merge_status>_mapDamage/  (graphiques MapDamage)
#     11_summary_tables/mapping_bwa_info_per_run.tsv (taux de mapping global)
# ==============================================================================

echo ""
echo "=========================================="
echo "MapDamage sur le nouveau pipeline per-run (k25/k29/k35/pracken + Bracken k35)"
echo "5 taxons cibles : mouton, chevre, olivier, frene, noisetier"
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
# CONFIGURATION GLOBALE — chemins adaptes au nouveau pipeline
# ==============================================================================

WORKDIR="/home/plstenge/coprolites_comparison"
OUT_ROOT="${WORKDIR}/11_per_run_analysis"

KRAKEN_K25_DIR="${OUT_ROOT}/05_kraken2_k25"
KRAKEN_K29_DIR="${OUT_ROOT}/05_kraken2_k29"
KRAKEN_K35_DIR="${OUT_ROOT}/05_kraken2_k35"
KRAKEN_PRACKEN_DIR="${OUT_ROOT}/05_kraken2_pracken"

BRACKEN_K35_DIR="${OUT_ROOT}/06_bracken/k35"

FASTP_DIR="${OUT_ROOT}/04_fastp"

KRAKENTOOLS_DIR="${WORKDIR}/08_krakentools/KrakenTools"

THREADS="${SLURM_CPUS_PER_TASK:-36}"

RUN_LABELS=("run1" "run2" "run3" "run4" "run5")
SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")

MAPDAMAGE_OUT_ROOT="${WORKDIR}/12_mapdamage_per_run"
mkdir -p "${MAPDAMAGE_OUT_ROOT}"

LOGFILE="${WORKDIR}/00_scripts/mapdamage_per_run_$(date +%Y%m%d_%H%M%S).txt"
touch "${LOGFILE}"

SUMMARY_DIR="${WORKDIR}/11_summary_tables"
mkdir -p "${SUMMARY_DIR}"
MAPPING_INFO="${SUMMARY_DIR}/mapping_bwa_info_per_run.tsv"
echo -e "Sample\tRun\tBase\tMergeStatus\tTaxon\tTotal_Reads\tMapped_Reads\tMapping_Rate(%)" > "${MAPPING_INFO}"

echo "$(date): Script MapDamage per-run demarre" | tee -a "${LOGFILE}"

# ==============================================================================
# GENOMES DE REFERENCE — uniquement les 5 taxons demandes
# ==============================================================================
# Format : ["Nom_espece"]="TaxID:/chemin/vers/genome.fa"
# ATTENTION : verifier/adapter les TaxIDs et chemins de genomes ci-dessous
# selon vos fichiers reels sur /home/plstenge/genomes/

declare -A TAXONS=(
    ["Ovis_aries"]="9940:/home/plstenge/genomes/Ovis_aries/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa"
    ["Capra_hircus"]="9925:/home/plstenge/genomes/Capra_hircus/Capra_hircus.ARS1.dna.toplevel.fa"
    ["Olea_europaea"]="4146:/home/plstenge/genomes/Olea_europaea/CACTIH01.fasta"
    ["Fraxinus_excelsior"]="38873:/home/plstenge/genomes/Fraxinus_excelsior/GCA_965226085.2_daFraExce3.hap1.2_genomic.fna"
    ["Corylus_avellana"]="13451:/home/plstenge/genomes/Corylus_avellana/Corylus_avellana_CavTom2PMs_1_0.fasta"
)

echo ""
echo "Genomes de reference configures: ${#TAXONS[@]}"
for sp in "${!TAXONS[@]}"; do
    echo "  - ${sp}"
done
echo ""

# ==============================================================================
# BASES KRAKEN2 / BRACKEN A TRAITER
# ==============================================================================
# Format : "LABEL:KRAKEN_ROOT_DIR"
# Note : "bracken_k35" pointe vers le meme repertoire kraken2 k35 pour
# l'extraction des reads, mais on ne traite l'unite que si un report
# Bracken existe (verifie dans la boucle) pour rester coherent avec Bracken.

KRAKEN_BASES=(
    "k25:${KRAKEN_K25_DIR}"
    "k29:${KRAKEN_K29_DIR}"
    "k35:${KRAKEN_K35_DIR}"
    "pracken:${KRAKEN_PRACKEN_DIR}"
    "bracken_k35:${KRAKEN_K35_DIR}"
)

MERGE_STATUSES=("merged" "unmerged")

################################################################################
# FONCTIONS
################################################################################

calculate_mapping_rate() {
    local bam_file="$1"
    local sample_name="$2"
    local run_label="$3"
    local base_label="$4"
    local merge_status="$5"
    local taxon="$6"
    local outfile="$7"

    if [[ -f "${bam_file}" ]]; then
        local total mapped rate
        total=$(samtools view -c "${bam_file}" 2>/dev/null || echo 0)
        mapped=$(samtools view -c -F 4 "${bam_file}" 2>/dev/null || echo 0)
        rate=0
        if [[ "${total}" -gt 0 ]]; then
            rate=$(echo "scale=1; ${mapped} * 100 / ${total}" | bc 2>/dev/null || echo "0")
        fi
        echo -e "${sample_name}\t${run_label}\t${base_label}\t${merge_status}\t${taxon}\t${total}\t${mapped}\t${rate}" >> "${outfile}"
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

check_bwa_index_complete() {
    local ref="$1"
    [[ -f "${ref}.bwt" && -f "${ref}.amb" && -f "${ref}.ann" && -f "${ref}.pac" && -f "${ref}.sa" ]]
}

################################################################################
# PHASE 1 : INDEXATION DES GENOMES DE REFERENCE
################################################################################

echo ""
echo "=========================================="
echo "PHASE 1: Indexation des genomes de reference (5 taxons)"
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

################################################################################
# PHASE 2 : EXTRACTION + MAPPING + MAPDAMAGE PAR RUN x SAMPLE x BASE x TAXON
################################################################################

echo ""
echo "=========================================="
echo "PHASE 2: Extraction reads Kraken2/Bracken + mapping BWA + MapDamage"
echo "=========================================="
echo ""

for BASE_ENTRY in "${KRAKEN_BASES[@]}"; do
    IFS=":" read -r BASE_LABEL KRAKEN_ROOT <<< "${BASE_ENTRY}"

    echo ""
    echo "###########################################################"
    echo "# Base : ${BASE_LABEL}"
    echo "###########################################################"

    for RUN_LABEL in "${RUN_LABELS[@]}"; do
        for SAMPLE in "${SAMPLES[@]}"; do

            UNIT="${SAMPLE}_${RUN_LABEL}"
            KRAKEN_DIR="${KRAKEN_ROOT}/${RUN_LABEL}/${SAMPLE}"

            for MERGE_STATUS in "${MERGE_STATUSES[@]}"; do

                KRAKEN_FILE="${KRAKEN_DIR}/${UNIT}_${MERGE_STATUS}.kraken"

                if [[ ! -s "${KRAKEN_FILE}" ]]; then
                    continue
                fi

                # Pour bracken_k35 : ne traiter l'unite que si un report
                # Bracken existe reellement pour ce sample/run/merge_status,
                # afin de rester coherent avec les abondances validees par
                # Bracken (et pas simplement les assignations brutes Kraken2).
                if [[ "${BASE_LABEL}" == "bracken_k35" ]]; then
                    BRACKEN_REPORT="${BRACKEN_K35_DIR}/${RUN_LABEL}/${SAMPLE}/${UNIT}_${MERGE_STATUS}.bracken"
                    if [[ ! -s "${BRACKEN_REPORT}" ]]; then
                        continue
                    fi
                fi

                # Fichier(s) fastq source correspondant au merge_status
                if [[ "${MERGE_STATUS}" == "merged" ]]; then
                    SEQ_FILE="${FASTP_DIR}/${RUN_LABEL}/${SAMPLE}/${UNIT}_fastp_merged.fastq.gz"
                    [[ -s "${SEQ_FILE}" ]] || continue
                else
                    SEQ_R1="${FASTP_DIR}/${RUN_LABEL}/${SAMPLE}/${UNIT}_fastp_unmerged_R1.fastq.gz"
                    SEQ_R2="${FASTP_DIR}/${RUN_LABEL}/${SAMPLE}/${UNIT}_fastp_unmerged_R2.fastq.gz"
                    [[ -s "${SEQ_R1}" && -s "${SEQ_R2}" ]] || continue
                fi

                echo ""
                echo "-----------------------------------------------------------------"
                echo " ${BASE_LABEL} | ${UNIT} | ${MERGE_STATUS}"
                echo "-----------------------------------------------------------------"

                for GROUP in "${!VALID_GENOMES[@]}"; do
                    REF="${VALID_GENOMES[${GROUP}]}"
                    TAX_ID="${TAXONS[${GROUP}]%%:*}"

                    OUTDIR="${MAPDAMAGE_OUT_ROOT}/${BASE_LABEL}/${MERGE_STATUS}/${GROUP}/${RUN_LABEL}/${SAMPLE}"
                    mkdir -p "${OUTDIR}"

                    echo "  -> ${GROUP} (TaxID: ${TAX_ID})"

                    if [[ "${MERGE_STATUS}" == "merged" ]]; then
                        OUT_FQ="${OUTDIR}/${UNIT}_${GROUP}_merged.fastq"

                        set +e
                        python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
                            -k "${KRAKEN_FILE}" \
                            -s "${SEQ_FILE}" \
                            -t "${TAX_ID}" \
                            -o "${OUT_FQ}" \
                            --fastq-output 2>>"${LOGFILE}"
                        set -e

                        if [[ ! -s "${OUT_FQ}" ]]; then
                            echo "     (aucun read extrait, on saute)"
                            rm -f "${OUT_FQ}" 2>/dev/null
                            continue
                        fi

                        READ_COUNT=$(grep -c "^@" "${OUT_FQ}" 2>/dev/null || echo 0)
                        echo "     ${READ_COUNT} reads extraits"

                        SAI="${OUTDIR}/${UNIT}_${GROUP}_merged.sai"
                        BAM="${OUTDIR}/${UNIT}_${GROUP}_merged.sorted.bam"

                        set +e
                        bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "${REF}" "${OUT_FQ}" > "${SAI}" 2>>"${LOGFILE}"
                        bwa samse "${REF}" "${SAI}" "${OUT_FQ}" 2>>"${LOGFILE}" | \
                            samtools view -bS - 2>>"${LOGFILE}" | \
                            samtools sort -o "${BAM}" - 2>>"${LOGFILE}"
                        set -e

                        [[ -f "${BAM}" ]] || { echo "     ERREUR mapping"; rm -f "${SAI}" "${OUT_FQ}"; continue; }

                        samtools index "${BAM}" 2>>"${LOGFILE}"

                        calculate_mapping_rate "${BAM}" "${SAMPLE}" "${RUN_LABEL}" "${BASE_LABEL}" "${MERGE_STATUS}" "${GROUP}" "${MAPPING_INFO}"

                        run_mapdamage "${BAM}" "${REF}" "${OUTDIR}/${UNIT}_${GROUP}_merged_mapDamage"

                        rm -f "${SAI}" "${OUT_FQ}"

                    else
                        OUT_R1="${OUTDIR}/${UNIT}_${GROUP}_unmerged_R1.fastq"
                        OUT_R2="${OUTDIR}/${UNIT}_${GROUP}_unmerged_R2.fastq"

                        set +e
                        python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
                            -k "${KRAKEN_FILE}" \
                            -s "${SEQ_R1}" -s2 "${SEQ_R2}" \
                            -t "${TAX_ID}" \
                            -o "${OUT_R1}" -o2 "${OUT_R2}" \
                            --fastq-output 2>>"${LOGFILE}"
                        set -e

                        if [[ ! -s "${OUT_R1}" || ! -s "${OUT_R2}" ]]; then
                            echo "     (aucun read extrait, on saute)"
                            rm -f "${OUT_R1}" "${OUT_R2}" 2>/dev/null
                            continue
                        fi

                        READ_COUNT=$(grep -c "^@" "${OUT_R1}" 2>/dev/null || echo 0)
                        echo "     ${READ_COUNT} paires extraites"

                        SAI1="${OUTDIR}/${UNIT}_${GROUP}_unmerged_R1.sai"
                        SAI2="${OUTDIR}/${UNIT}_${GROUP}_unmerged_R2.sai"
                        BAM="${OUTDIR}/${UNIT}_${GROUP}_unmerged.sorted.bam"

                        set +e
                        bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "${REF}" "${OUT_R1}" > "${SAI1}" 2>>"${LOGFILE}"
                        bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "${REF}" "${OUT_R2}" > "${SAI2}" 2>>"${LOGFILE}"
                        bwa sampe "${REF}" "${SAI1}" "${SAI2}" "${OUT_R1}" "${OUT_R2}" 2>>"${LOGFILE}" | \
                            samtools view -bS - 2>>"${LOGFILE}" | \
                            samtools sort -o "${BAM}" - 2>>"${LOGFILE}"
                        set -e

                        [[ -f "${BAM}" ]] || { echo "     ERREUR mapping"; rm -f "${SAI1}" "${SAI2}" "${OUT_R1}" "${OUT_R2}"; continue; }

                        samtools index "${BAM}" 2>>"${LOGFILE}"

                        calculate_mapping_rate "${BAM}" "${SAMPLE}" "${RUN_LABEL}" "${BASE_LABEL}" "${MERGE_STATUS}" "${GROUP}" "${MAPPING_INFO}"

                        run_mapdamage "${BAM}" "${REF}" "${OUTDIR}/${UNIT}_${GROUP}_unmerged_mapDamage"

                        rm -f "${SAI1}" "${SAI2}" "${OUT_R1}" "${OUT_R2}"
                    fi
                done
            done
        done
    done
done

################################################################################
# RESUME ET FIN
################################################################################

echo ""
echo "=========================================="
echo "PIPELINE MAPDAMAGE PER-RUN TERMINE"
echo "Date: $(date)"
echo "=========================================="
echo ""
echo "Statistiques de mapping : ${MAPPING_INFO}"
echo "Log complet             : ${LOGFILE}"
echo "Resultats BAM+MapDamage : ${MAPDAMAGE_OUT_ROOT}/<base>/<merge_status>/<taxon>/<run>/<sample>/"
echo ""

if [[ -f "${MAPPING_INFO}" ]]; then
    echo "Apercu des resultats :"
    column -t "${MAPPING_INFO}" | head -50
fi

conda deactivate

echo "Script termine."
echo "Pipeline MapDamage per-run termine le $(date +%Y-%m-%d\ %H:%M:%S). Resultats: ${MAPDAMAGE_OUT_ROOT}" | \
    mail -s "Pipeline Coprolites MapDamage per-run - Termine" pierrelouis.stenger@gmail.com 2>/dev/null || true

exit 0
