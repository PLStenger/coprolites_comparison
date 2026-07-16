#!/bin/bash

#SBATCH --job-name=02_bbduk_cop
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=8
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/02_bbduk.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/02_bbduk.out"

# ==============================================================================
# Script : 02_bbduk.sh
# Author : Pierre-Louis Stenger
# Purpose :
# 1. Trimming/décontamination des reads avec BBDuk (phiX + adaptateurs)
#    sur les fichiers groupés (concaténés) de chaque échantillon
# 2. FastQC sur les fichiers nettoyés
# 3. MultiQC sur les résultats FastQC post-trim
# ==============================================================================

# set -euo pipefail

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

module load conda/4.12.0
source ~/.bashrc
#conda activate bioinformatic

# ==============================================================================
# PATHS
# ==============================================================================

WORKDIR="/home/plstenge/coprolites_comparison"

# Répertoire d'entrée : fichiers groupés issus du script 01
ENTREE="${WORKDIR}/01_raw_data/grouped"

# Répertoire de sortie des fichiers nettoyés par BBDuk
SORTIE="${WORKDIR}/03_bbduk"

# Répertoire pour le contrôle qualité post-trim (FastQC + MultiQC)
QUALITE="${WORKDIR}/03_bbduk/controle_qualite"

# Threads
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ==============================================================================
# BBDUK SETTINGS
# ==============================================================================

#PHIX="phix"
PHIX="/home/plstenge/bbmap/resources/phix174_ill.ref.fa.gz"

#BBDUK="bbduk.sh"
BBDUK="/home/plstenge/bbmap/bbduk.sh"

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
# STEP 1 : BBDuk sur chaque fichier groupé
#
# Pattern d'entrée : cop408_grouped_R1.fastq.gz / cop408_grouped_R2.fastq.gz
# Pattern de sortie : clean_cop408_grouped_R1.fastq.gz / clean_cop408_grouped_R2.fastq.gz
#
# Paramètres (identiques au script de référence) :
#   ref=phix    → décontamination phiX
#   ktrim=rl    → trim des deux côtés (right + left)
#   k=23        → taille de k-mer
#   mink=11     → k-mer minimal en bout de read
#   hdist=1     → 1 mismatch autorisé
#   stats=...   → fichier de statistiques par échantillon
# ==============================================================================

echo "======================================================================="
echo " STEP 1: BBDuk trimming sur les fichiers groupés"
echo "======================================================================="

echo "DÉBUT DU TRAITEMENT BBDUK"

for r1_file in "${ENTREE}"/*_R1.fastq.gz; do

    # Dériver le chemin R2 correspondant
    r2_file="${r1_file/_R1.fastq.gz/_R2.fastq.gz}"

    if [[ ! -f "${r2_file}" ]]; then
        echo "  ERREUR: Fichier R2 manquant pour ${r1_file}" >&2
        continue
    fi

    # Extraire le nom de base de l'échantillon, ex : cop408_grouped
    file_name=$(basename "${r1_file}" _R1.fastq.gz)

    echo ""
    echo "  Traitement de l'échantillon : ${file_name}"
    echo "    R1 : ${r1_file}"
    echo "    R2 : ${r2_file}"

    "${BBDUK}" \
        in1="${r1_file}" \
        in2="${r2_file}" \
        out1="${SORTIE}/clean_${file_name}_R1.fastq.gz" \
        out2="${SORTIE}/clean_${file_name}_R2.fastq.gz" \
        ref="${PHIX}" \
        ktrim=rl \
        k=23 \
        mink=11 \
        hdist=1 \
        threads="${THREADS}" \
        stats="${SORTIE}/${file_name}_bbduk_stats.txt"

    echo "  >>> BBDuk terminé pour l'échantillon : ${file_name} <<<"

done

echo ""
echo " STEP 1 done."
echo "Analyse BBDuk terminée"

# ==============================================================================
# STEP 2 : FastQC sur les fichiers nettoyés
# ==============================================================================

echo "======================================================================="
echo " STEP 2: FastQC sur les fichiers nettoyés (post-BBDuk)"
echo "======================================================================="

echo "Lancement de FastQC"

for f in "${SORTIE}"/clean_*.fastq.gz; do
    echo "  [FastQC] ${f}"
    fastqc \
        --threads "${THREADS}" \
        --outdir "${QUALITE}" \
        "${f}"
done

echo ""
echo " STEP 2 done. Rapports FastQC dans : ${QUALITE}"

# ==============================================================================
# STEP 3 : MultiQC sur les résultats FastQC post-trim + stats BBDuk
# ==============================================================================

echo "======================================================================="
echo " STEP 3: MultiQC sur les résultats post-BBDuk"
echo "======================================================================="

echo "Lancement de MultiQC"

multiqc \
    "${QUALITE}" \
    "${SORTIE}" \
    --outdir "${QUALITE}" \
    --filename "multiqc_report_post_bbduk" \
    --title "MultiQC — Post-BBDuk trimming (coprolites)" \
    --force

echo " STEP 3 done. Rapport MultiQC dans : ${QUALITE}/multiqc_report_post_bbduk.html"

# ==============================================================================
# RÉSUMÉ FINAL
# ==============================================================================

echo ""
echo "======================================================================="
echo " ALL STEPS COMPLETED SUCCESSFULLY"
echo "======================================================================="
echo ""
echo "  Fichiers nettoyés (BBDuk)  : ${SORTIE}"
echo "  Stats BBDuk par échantillon: ${SORTIE}/*_bbduk_stats.txt"
echo "  FastQC post-trim           : ${QUALITE}"
echo "  MultiQC post-trim          : ${QUALITE}/multiqc_report_post_bbduk.html"
echo "======================================================================="
