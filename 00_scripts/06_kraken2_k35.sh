#!/bin/bash

#SBATCH --job-name=07_kraken2_k35_cop
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/07_kraken2_k35.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/07_kraken2_k35.out"

# ==============================================================================
# Script : 07_kraken2_k35.sh
# Author : Pierre-Louis Stenger
# Purpose :
# 1. Classification taxonomique avec Kraken2 (base de données k35, standard)
#    sur les reads mergés (single-end) et non-mergés (paired-end) issus de fastp
# 2. Visualisation Krona (merged / unmerged séparément)
# ==============================================================================

# set -euo pipefail

module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics

# ==============================================================================
# PATHS
# ==============================================================================

WORKDIR="/home/plstenge/coprolites_comparison"

ENTREE="${WORKDIR}/06_fastp"
KRAKEN2_DB="/storage/groups/gdec/shared/Kraken_database/k2_core_nt_20250609"
SORTIE="${WORKDIR}/07_kraken2_k35"
SORTIE_KRONA="${SORTIE}/krona"

THREADS="${SLURM_CPUS_PER_TASK:-36}"

mkdir -p "${SORTIE}"
mkdir -p "${SORTIE_KRONA}"

# ==============================================================================
# STEP 1 : Kraken2 sur les reads mergés (single-end)
# ==============================================================================

echo "======================================================================="
echo " STEP 1: Kraken2 (k35) — reads merged (single-end)"
echo "======================================================================="

for MERGED in "${ENTREE}"/*_fastp_merged.fastq.gz; do

    SAMPLE=$(basename "${MERGED}" _fastp_merged.fastq.gz)
    SORTIE_KRAKEN="${SORTIE}/${SAMPLE}_merged.kraken"
    SORTIE_REPORT="${SORTIE}/${SAMPLE}_merged.report"

    echo "  Traitement : ${SAMPLE} (merged)"

    kraken2 --conf 0.2 --db "${KRAKEN2_DB}" --threads "${THREADS}" \
        --output "${SORTIE_KRAKEN}" --report "${SORTIE_REPORT}" "${MERGED}"

    echo "  >>> Termine: ${SAMPLE} (merged) <<<"

done

# ==============================================================================
# STEP 2 : Kraken2 sur les reads non-mergés (paired-end)
# ==============================================================================

echo "======================================================================="
echo " STEP 2: Kraken2 (k35) — reads unmerged (paired-end)"
echo "======================================================================="

for R1 in "${ENTREE}"/*_fastp_unmerged_R1.fastq.gz; do

    SAMPLE=$(basename "${R1}" _fastp_unmerged_R1.fastq.gz)
    R2="${ENTREE}/${SAMPLE}_fastp_unmerged_R2.fastq.gz"

    if [[ -f "${R2}" ]]; then
        SORTIE_KRAKEN="${SORTIE}/${SAMPLE}_unmerged.kraken"
        SORTIE_REPORT="${SORTIE}/${SAMPLE}_unmerged.report"

        echo "  Traitement : ${SAMPLE} (unmerged)"

        kraken2 --conf 0.2 --paired --db "${KRAKEN2_DB}" --threads "${THREADS}" \
            --output "${SORTIE_KRAKEN}" --report "${SORTIE_REPORT}" "${R1}" "${R2}"

        echo "  >>> Termine : ${SAMPLE} (unmerged) <<<"
    else
        echo "  ATTENTION: R2 manquant pour ${SAMPLE}, ignore."
    fi

done

echo ""
echo "Analyse Kraken2 (k35) terminee pour tous les echantillons."

# ==============================================================================
# STEP 3 : Visualisation Krona
# ==============================================================================

echo "======================================================================="
echo " STEP 3: Generation des graphiques Krona (k35)"
echo "======================================================================="

echo "  Generation du Krona pour les reads merged"
ktImportTaxonomy \
    -t 5 \
    -m 3 \
    -o "${SORTIE_KRONA}/krona_merged.html" \
    "${SORTIE}"/*_merged.report

echo "  Generation du Krona pour les reads unmerged"
ktImportTaxonomy \
    -t 5 \
    -m 3 \
    -o "${SORTIE_KRONA}/krona_unmerged.html" \
    "${SORTIE}"/*_unmerged.report

echo ""
echo "======================================================================="
echo " ALL STEPS COMPLETED — Kraken2 k35"
echo "======================================================================="
echo "  Resultats Kraken2 : ${SORTIE}"
echo "  Krona merged      : ${SORTIE_KRONA}/krona_merged.html"
echo "  Krona unmerged     : ${SORTIE_KRONA}/krona_unmerged.html"
echo "======================================================================="
