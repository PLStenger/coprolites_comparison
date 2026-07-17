#!/bin/bash

#SBATCH --job-name=05_fastp_cop
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/05_fastp.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/05_fastp.out"

# ==============================================================================
# Script : 05_fastp.sh
# Author : Pierre-Louis Stenger
# Purpose :
# 1. Trimming qualité, filtrage et merge des paires avec fastp
#    sur les fichiers issus de Clumpify
# 2. FastQC sur les reads mergés (.fastq.gz)
# 3. MultiQC sur les résultats FastQC post-fastp
# ==============================================================================

# set -euo pipefail

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics

# ==============================================================================
# PATHS
# ==============================================================================

WORKDIR="/home/plstenge/coprolites_comparison"

# Répertoire d'entrée : fichiers issus de Clumpify
ENTREE="${WORKDIR}/05_clumpify"

# Répertoire de sortie fastp
SORTIE="${WORKDIR}/06_fastp"

# Répertoire pour le contrôle qualité post-fastp
QUALITE="${WORKDIR}/06_fastp/controle_qualite"

# Threads
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ==============================================================================
# STEP 0 : Créer les répertoires de sortie
# ==============================================================================

echo "======================================================================="
echo " STEP 0: Creating output directories"
echo "======================================================================="

mkdir -p "${SORTIE}"
mkdir -p "${QUALITE}"

echo " Done."

# ==============================================================================
# STEP 1 : fastp — trimming qualité + merge des paires
#
# Pattern d'entrée : clean_cop408_grouped_dedup_clumpify_R1.fastq.gz
# Sorties par échantillon :
#   - unmerged R1/R2 (paires qui n'ont pas pu être fusionnées)
#   - merged (reads fusionnés, principal produit utile en aval)
#   - rapport HTML/JSON
#
# Paramètres (identiques au script de référence) :
#   --length_required 20
#   --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 10
#   --n_base_limit 5
#   --unqualified_percent_limit 40
#   --complexity_threshold 30 --low_complexity_filter
#   --qualified_quality_phred 20
#   --trim_poly_x --poly_x_min_len 10
#   --merge --correction
#   --overlap_len_require 10 --overlap_diff_limit 5 --overlap_diff_percent_limit 20
#   --adapter_sequence / --adapter_sequence_r2 (adaptateurs Illumina standards)
#   --detect_adapter_for_pe
# ==============================================================================

echo "======================================================================="
echo " STEP 1: fastp — trimming qualité et merge des paires"
echo "======================================================================="

for R1 in "${ENTREE}"/*_clumpify_R1.fastq.gz; do

    base=$(basename "${R1}" _clumpify_R1.fastq.gz)
    R2="${ENTREE}/${base}_clumpify_R2.fastq.gz"

    if [[ ! -f "${R2}" ]]; then
        echo "  Fichier R2 manquant pour ${R1}, échantillon ignoré." >&2
        continue
    fi

    echo ""
    echo "  Traitement fastp pour : ${base}"
    echo "    R1 : ${R1}"
    echo "    R2 : ${R2}"

    OUT_R1="${SORTIE}/${base}_fastp_unmerged_R1.fastq.gz"
    OUT_R2="${SORTIE}/${base}_fastp_unmerged_R2.fastq.gz"
    MERGED="${SORTIE}/${base}_fastp_merged.fastq.gz"
    HTML="${SORTIE}/${base}_fastp_report.html"
    JSON="${SORTIE}/${base}_fastp_report.json"

    fastp \
        --in1 "${R1}" --in2 "${R2}" \
        --out1 "${OUT_R1}" --out2 "${OUT_R2}" \
        --merged_out "${MERGED}" \
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
        --html "${HTML}" \
        --json "${JSON}" \
        --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --detect_adapter_for_pe \
        --thread "${THREADS}"

    echo "  >>> fastp terminé pour l'échantillon : ${base} <<<"

done

echo ""
echo " STEP 1 done."
echo "Tous les traitements fastp sont terminés."

# ==============================================================================
# STEP 2 : FastQC sur les reads mergés (.fastq.gz), fichier par fichier
# ==============================================================================

echo "======================================================================="
echo " STEP 2: FastQC sur les reads mergés (post-fastp)"
echo "======================================================================="

echo "Lancement de FastQC (fichier par fichier)..."

for fq in "${SORTIE}"/*_merged.fastq.gz; do
    if [[ -f "${fq}" ]]; then
        echo "  --> FastQC sur : $(basename "${fq}")"
        fastqc \
            --threads "${THREADS}" \
            --outdir "${QUALITE}" \
            "${fq}"
    fi
done

echo ""
echo " STEP 2 done. Rapports FastQC dans : ${QUALITE}"

# ==============================================================================
# STEP 3 : MultiQC sur les résultats FastQC + rapports fastp (JSON)
# ==============================================================================

echo "======================================================================="
echo " STEP 3: MultiQC sur les résultats post-fastp"
echo "======================================================================="

echo "Lancement de MultiQC..."

multiqc \
    "${QUALITE}" \
    "${SORTIE}" \
    --outdir "${QUALITE}" \
    --filename "multiqc_report_post_fastp" \
    --title "MultiQC — Post-fastp trimming and merging (coprolites)" \
    --force

echo " STEP 3 done. Rapport MultiQC dans : ${QUALITE}/multiqc_report_post_fastp.html"

# ==============================================================================
# RÉSUMÉ FINAL
# ==============================================================================

echo ""
echo "======================================================================="
echo " ALL STEPS COMPLETED SUCCESSFULLY"
echo "======================================================================="
echo ""
echo "  Reads mergés fastp        : ${SORTIE}/*_fastp_merged.fastq.gz"
echo "  Reads non-mergés (paires) : ${SORTIE}/*_fastp_unmerged_R1/R2.fastq.gz"
echo "  Rapports fastp (HTML/JSON): ${SORTIE}/*_fastp_report.*"
echo "  FastQC post-fastp         : ${QUALITE}"
echo "  MultiQC post-fastp        : ${QUALITE}/multiqc_report_post_fastp.html"
echo "======================================================================="
