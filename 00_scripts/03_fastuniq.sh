#!/bin/bash

#SBATCH --job-name=03_fastuniq_cop
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/03_fastuniq.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/03_fastuniq.out"

# ==============================================================================
# Script : 03_fastuniq.sh
# Author : Pierre-Louis Stenger
# Purpose :
# 1. Déduplication des reads avec FastUniq (nécessite décompression temporaire)
# 2. FastQC sur les fichiers dédupliqués (fichier par fichier)
# 3. MultiQC sur les résultats FastQC post-dedup
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

# Répertoire d'entrée : fichiers nettoyés par BBDuk
ENTREE="${WORKDIR}/03_bbduk"

# Répertoire de sortie des fichiers dédupliqués
SORTIE="${WORKDIR}/04_fastuniq"

# Répertoire pour le contrôle qualité post-dedup (FastQC + MultiQC)
QUALITE="${WORKDIR}/04_fastuniq/controle_qualite"

# Répertoire temporaire pour décompression (FastUniq n'accepte pas le .gz)
TMP="${WORKDIR}/04_fastuniq/tmp"

# Threads pour FastQC
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ==============================================================================
# STEP 0 : Créer les répertoires de sortie
# ==============================================================================

echo "======================================================================="
echo " STEP 0: Creating output directories"
echo "======================================================================="

mkdir -p "${SORTIE}"
mkdir -p "${QUALITE}"
mkdir -p "${TMP}"

echo " Done."

# ==============================================================================
# STEP 1 : FastUniq — déduplication par paire d'échantillon
#
# FastUniq ne lit pas les .fastq.gz → décompression temporaire dans TMP,
# puis suppression immédiate après traitement pour économiser l'espace disque.
#
# Paramètres (identiques au script de référence) :
#   -i listfile  → fichier texte listant R1 et R2 décompressés (un par ligne)
#   -t q         → format FASTQ en entrée
#   -o           → fichier de sortie R1 dédupliqué
#   -p           → fichier de sortie R2 dédupliqué
# ==============================================================================

echo "======================================================================="
echo " STEP 1: FastUniq — déduplication des paires de reads"
echo "======================================================================="

cd "${ENTREE}" || exit 1

for R1_gz in clean_*_R1.fastq.gz; do

    # Dériver le nom de base et le fichier R2 correspondant
    base=$(echo "${R1_gz}" | sed 's/_R1\.fastq\.gz//')
    R2_gz="${base}_R2.fastq.gz"

    if [[ ! -f "${R2_gz}" ]]; then
        echo "  ATTENTION: fichier R2 manquant pour ${base} — ignoré." >&2
        continue
    fi

    echo ""
    echo "  Traitement de la paire : ${base}"
    echo "    R1 : ${ENTREE}/${R1_gz}"
    echo "    R2 : ${ENTREE}/${R2_gz}"

    # Chemins des fichiers temporaires décompressés
    R1_tmp="${TMP}/${base}_R1.fastq"
    R2_tmp="${TMP}/${base}_R2.fastq"
    listfile="${TMP}/${base}.list"

    # Décompression temporaire
    echo "    Décompression temporaire..."
    gzip -dc "${ENTREE}/${R1_gz}" > "${R1_tmp}"
    gzip -dc "${ENTREE}/${R2_gz}" > "${R2_tmp}"

    # Création du fichier liste pour FastUniq
    echo -e "${R1_tmp}\n${R2_tmp}" > "${listfile}"

    # Lancement de FastUniq
    echo "    Lancement FastUniq..."
    fastuniq \
        -i "${listfile}" \
        -t q \
        -o "${SORTIE}/${base}_dedup_R1.fastq" \
        -p "${SORTIE}/${base}_dedup_R2.fastq"

    # Nettoyage des fichiers temporaires immédiatement après traitement
    rm -f "${R1_tmp}" "${R2_tmp}" "${listfile}"

    echo "  >>> FastUniq terminé pour l'échantillon : ${base} <<<"

done

# Nettoyage du répertoire tmp (normalement vide à ce stade)
rm -rf "${TMP}"

echo ""
echo " STEP 1 done."
echo "FastUniq terminé — tous les échantillons traités."

# ==============================================================================
# STEP 2 : FastQC sur les fichiers dédupliqués (fichier par fichier)
#
# Note : FastUniq produit des .fastq (non compressés)
# On traite fichier par fichier pour éviter les pics mémoire (cf. script réf.)
# ==============================================================================

echo "======================================================================="
echo " STEP 2: FastQC sur les fichiers dédupliqués (post-FastUniq)"
echo "======================================================================="

echo "Lancement de FastQC (fichier par fichier)..."

for fq in "${SORTIE}"/*.fastq; do
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
# STEP 3 : MultiQC sur les résultats FastQC post-dedup
#
# On se place dans $QUALITE pour éviter que MultiQC ne parcoure les gros
# fichiers .fastq de $SORTIE (cf. comportement du script de référence).
# ==============================================================================

echo "======================================================================="
echo " STEP 3: MultiQC sur les résultats post-FastUniq"
echo "======================================================================="

echo "Lancement de MultiQC..."

cd "${QUALITE}" || exit 1

multiqc \
    . \
    --outdir . \
    --filename "multiqc_report_post_fastuniq" \
    --title "MultiQC — Post-FastUniq deduplication (coprolites)" \
    --force

echo " STEP 3 done. Rapport MultiQC dans : ${QUALITE}/multiqc_report_post_fastuniq.html"

# ==============================================================================
# RÉSUMÉ FINAL
# ==============================================================================

echo ""
echo "======================================================================="
echo " ALL STEPS COMPLETED SUCCESSFULLY"
echo "======================================================================="
echo ""
echo "  Fichiers dédupliqués (FastUniq) : ${SORTIE}"
echo "  FastQC post-dedup               : ${QUALITE}"
echo "  MultiQC post-dedup              : ${QUALITE}/multiqc_report_post_fastuniq.html"
echo "======================================================================="
