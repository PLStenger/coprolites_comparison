#!/bin/bash

#SBATCH --job-name=04_clumpify_cop
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/04_clumpify.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/04_clumpify.out"

# ==============================================================================
# Script : 04_clumpify.sh
# Author : Pierre-Louis Stenger
# Purpose :
# 1. Deduplication optique et tri des reads par k-mer avec Clumpify
#    sur les fichiers dédupliqués par FastUniq (.fastq non compressés)
# 2. FastQC sur les fichiers produits par Clumpify (.fastq.gz)
# 3. MultiQC sur les résultats FastQC post-clumpify
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

# Répertoire d'entrée : fichiers dédupliqués par FastUniq (.fastq non compressés)
ENTREE="${WORKDIR}/04_fastuniq"

# Répertoire de sortie des fichiers Clumpify (.fastq.gz)
SORTIE="${WORKDIR}/05_clumpify"

# Répertoire pour le contrôle qualité post-clumpify
QUALITE="${WORKDIR}/05_clumpify/controle_qualite"

# Threads
THREADS="${SLURM_CPUS_PER_TASK:-8}"

CLUMPIFY="clumpify.sh"

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
# STEP 1 : Clumpify — deduplication optique + tri par k-mer
#
# Pattern d'entrée  : clean_cop408_grouped_dedup_R1.fastq
# Pattern de sortie : clean_cop408_grouped_dedup_clumpify_R1.fastq.gz
#
# Paramètres (identiques au script de référence) :
#   dedupe=t  → activer la déduplication des reads identiques
#
# Note : Clumpify accepte directement les .fastq non compressés en entrée
#        et produit des .fastq.gz en sortie.
# ==============================================================================

echo "======================================================================="
echo " STEP 1: Clumpify — deduplication optique et tri par k-mer"
echo "======================================================================="

for R1 in "${ENTREE}"/*_dedup_R1.fastq; do

    # Dériver le fichier R2 correspondant
    R2="${R1/_dedup_R1.fastq/_dedup_R2.fastq}"

    if [[ ! -f "${R2}" ]]; then
        echo "  ATTENTION: fichier R2 manquant pour ${R1} — ignoré." >&2
        continue
    fi

    # Nom de base : ex. clean_cop408_grouped_dedup
    base=$(basename "${R1}" _R1.fastq)

    echo ""
    echo "  Traitement de la paire : ${base}"
    echo "    R1 : ${R1}"
    echo "    R2 : ${R2}"

    "${CLUMPIFY}" \
        in="${R1}" \
        in2="${R2}" \
        out="${SORTIE}/${base}_clumpify_R1.fastq.gz" \
        out2="${SORTIE}/${base}_clumpify_R2.fastq.gz" \
        dedupe=t \
        threads="${THREADS}"

    echo "  >>> Clumpify terminé pour l'échantillon : ${base} <<<"

done

echo ""
echo " STEP 1 done."
echo "Analyse Clumpify terminée."

# ==============================================================================
# STEP 2 : FastQC sur les fichiers Clumpify (.fastq.gz), fichier par fichier
# ==============================================================================

echo "======================================================================="
echo " STEP 2: FastQC sur les fichiers post-Clumpify"
echo "======================================================================="

echo "Lancement de FastQC (fichier par fichier)..."

for fq in "${SORTIE}"/*.fastq.gz; do
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
# STEP 3 : MultiQC sur les résultats FastQC post-Clumpify
# ==============================================================================

echo "======================================================================="
echo " STEP 3: MultiQC sur les résultats post-Clumpify"
echo "======================================================================="

echo "Lancement de MultiQC..."

multiqc \
    "${QUALITE}" \
    --outdir "${QUALITE}" \
    --filename "multiqc_report_post_clumpify" \
    --title "MultiQC — Post-Clumpify deduplication (coprolites)" \
    --force

echo " STEP 3 done. Rapport MultiQC dans : ${QUALITE}/multiqc_report_post_clumpify.html"

# ==============================================================================
# RÉSUMÉ FINAL
# ==============================================================================

echo ""
echo "======================================================================="
echo " ALL STEPS COMPLETED SUCCESSFULLY"
echo "======================================================================="
echo ""
echo "  Fichiers Clumpify          : ${SORTIE}"
echo "  FastQC post-clumpify       : ${QUALITE}"
echo "  MultiQC post-clumpify      : ${QUALITE}/multiqc_report_post_clumpify.html"
echo "======================================================================="
