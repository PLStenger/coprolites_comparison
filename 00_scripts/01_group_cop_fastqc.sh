#!/bin/bash

#SBATCH --job-name=01_group_cop_files
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=8
#SBATCH --mem=700G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/01_group_cop_files.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/01_group_cop_files.out"

# ==============================================================================
# Script : 01_group_cop_fastqc.sh
# Author : Pierre-Louis Stenger
# Purpose :
# 1. Concatenate (cat) R1 and R2 reads per sample across all sequencing runs
#    (RUN1 inclus — structure flat : fichiers directement dans le dossier)
# 2. Run FastQC sur chaque fichier FASTQ individuel, par run et par échantillon
# 3. Run MultiQC par échantillon et par run
# 4. Run FastQC sur les fichiers groupés (concaténés)
# 5. Run MultiQC sur les fichiers groupés
# ==============================================================================

# set -euo pipefail

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

module load conda/4.12.0
source ~/.bashrc
conda activate bioinformatic

# ==============================================================================
# PATHS — Input run directories
# ==============================================================================

# RUN1: structure plate — fichiers directement dans le dossier
# Pattern: cop408_R1.fastq.gz  (pas de sous-dossier numéroté)
RUN1="/storage/groups/gdec/shared_paleo/Illumina/01_raw_data"

# RUN2–RUN5: structure sous-dossier — fichiers dans 474_cop408/ etc.
# Pattern: 474_cop408/474_cop408_R1.fastq.gz
RUN2="/storage/groups/gdec/shared_paleo/E1531_final/run1_20250320_AV241601_E1531_Ps5Lane1_Ps6Lane2"
RUN3="/storage/groups/gdec/shared_paleo/E1531_final/run2_20250414_AV241601_E1531_Ps5_Ps6_11022026_CORRECTED"
RUN4="/storage/groups/gdec/shared_paleo/E1531_final/run3_20251008_AV241601_E1531_Ps7_Ps8"
RUN5="/storage/groups/gdec/shared_paleo/E1531_final/run4_20251104_AV241601_E1531_Ps7_Ps8_04112025"

# Labels humains pour chaque run (utilisés dans les noms de dossiers MultiQC)
RUN_LABELS=("RUN1" "RUN2" "RUN3" "RUN4" "RUN5")

# Tableau des chemins dans le même ordre que RUN_LABELS
ALL_RUNS=("${RUN1}" "${RUN2}" "${RUN3}" "${RUN4}" "${RUN5}")

# ==============================================================================
# OUTPUT PATHS
# ==============================================================================

WORKDIR="/home/plstenge/coprolites_comparison"

# Répertoire des fichiers FASTQ concaténés
GROUPED_DIR="${WORKDIR}/01_raw_data/grouped"

# FastQC par séquence individuelle (un rapport par fichier FASTQ, par run)
FASTQC_RAW_DIR="${WORKDIR}/02_quality_check_raw/fastqc_per_run"

# FastQC sur les fichiers groupés
FASTQC_GROUPED_DIR="${WORKDIR}/02_quality_check_raw/fastqc_grouped"

# MultiQC par échantillon × run (sous-dossier par combinaison)
MULTIQC_RAW_DIR="${WORKDIR}/02_quality_check_raw/multiqc_per_run"

# MultiQC sur les fichiers groupés
MULTIQC_GROUPED_DIR="${WORKDIR}/02_quality_check_raw/multiqc_grouped"

# Threads FastQC
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ==============================================================================
# SAMPLES
# Nom court de l'échantillon tel qu'il apparaît dans le nom de fichier
# ==============================================================================

SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")

# ==============================================================================
# HELPER FUNCTION : find_fastq
# Cherche un fichier FASTQ pour un échantillon et un sens de lecture donnés
# dans un répertoire run (structure plate ou 1 niveau de sous-dossier).
# Retourne le chemin complet, ou chaîne vide si non trouvé.
# ==============================================================================

find_fastq() {
    local run_dir="$1"
    local sample="$2"   # e.g. "cop408"
    local read="$3"     # "R1" ou "R2"

    find "${run_dir}" -maxdepth 2 -type f \
        -name "*${sample}*${read}*.fastq.gz" \
        2>/dev/null | head -n 1
}

# ==============================================================================
# STEP 0 : Créer les répertoires de sortie
# ==============================================================================

echo "======================================================================="
echo " STEP 0: Creating output directories"
echo "======================================================================="

mkdir -p "${GROUPED_DIR}"
mkdir -p "${FASTQC_RAW_DIR}"
mkdir -p "${FASTQC_GROUPED_DIR}"
mkdir -p "${MULTIQC_RAW_DIR}"
mkdir -p "${MULTIQC_GROUPED_DIR}"

echo " Done."

# ==============================================================================
# STEP 1 : Concaténation des FASTQ par échantillon à travers tous les runs
#
# Principe :
#   - Pour chaque échantillon, on cherche ses fichiers R1/R2 dans chacun des
#     5 runs (RUN1 compris, avec sa structure plate).
#   - On cat tous les fichiers trouvés dans un seul fichier groupé.
#   - Un echo final confirme la fin de la concaténation pour cet échantillon.
# ==============================================================================

echo "======================================================================="
echo " STEP 1: Concatenating FASTQ files per sample across all runs (RUN1–RUN5)"
echo "======================================================================="

for SAMPLE in "${SAMPLES[@]}"; do

    echo ""
    echo "  Sample: ${SAMPLE}"

    R1_FILES=()
    R2_FILES=()

    for i in "${!ALL_RUNS[@]}"; do
        RUN_DIR="${ALL_RUNS[$i]}"
        RUN_LABEL="${RUN_LABELS[$i]}"

        r1=$(find_fastq "${RUN_DIR}" "${SAMPLE}" "R1")
        r2=$(find_fastq "${RUN_DIR}" "${SAMPLE}" "R2")

        if [[ -n "${r1}" && -n "${r2}" ]]; then
            echo "    [FOUND] ${RUN_LABEL} : ${RUN_DIR}"
            echo "      R1 : ${r1}"
            echo "      R2 : ${r2}"
            R1_FILES+=("${r1}")
            R2_FILES+=("${r2}")
        else
            echo "    [SKIP]  ${RUN_LABEL} : ${SAMPLE} non présent"
        fi
    done

    if [[ ${#R1_FILES[@]} -eq 0 ]]; then
        echo "    WARNING: ${SAMPLE} — aucun fichier trouvé dans aucun run. Ignoré."
        continue
    fi

    OUT_R1="${GROUPED_DIR}/${SAMPLE}_grouped_R1.fastq.gz"
    OUT_R2="${GROUPED_DIR}/${SAMPLE}_grouped_R2.fastq.gz"

    echo "    Concaténation de ${#R1_FILES[@]} fichier(s) R1 -> ${OUT_R1}"
    cat "${R1_FILES[@]}" > "${OUT_R1}"

    echo "    Concaténation de ${#R2_FILES[@]} fichier(s) R2 -> ${OUT_R2}"
    cat "${R2_FILES[@]}" > "${OUT_R2}"

    echo "    >>> CONCATENATION TERMINEE pour ${SAMPLE} <<<"

done

echo ""
echo " STEP 1 done. Fichiers groupés dans : ${GROUPED_DIR}"
echo " >>> TOUTES LES CONCATENATIONS SONT TERMINEES — vous pouvez lancer des analyses en parallèle sur ${GROUPED_DIR} <<<"

# ==============================================================================
# STEP 2 : FastQC sur chaque fichier FASTQ individuel, par run et par échantillon
#
# Principe :
#   - Pour chaque run et chaque échantillon, on cherche les fichiers R1 et R2.
#   - On lance un FastQC INDIVIDUEL sur chaque fichier (pas en batch).
#     Cela donne un rapport par fichier, ex :
#       fastqc_per_run/RUN1/cop408/cop408_R1_fastqc.html
#       fastqc_per_run/RUN1/cop408/cop408_R2_fastqc.html
#   - On range les résultats dans un sous-dossier par run/échantillon.
# ==============================================================================

echo "======================================================================="
echo " STEP 2: FastQC par fichier individuel (par run × échantillon)"
echo "======================================================================="

for i in "${!ALL_RUNS[@]}"; do
    RUN_DIR="${ALL_RUNS[$i]}"
    RUN_LABEL="${RUN_LABELS[$i]}"

    echo ""
    echo "  Run: ${RUN_LABEL} (${RUN_DIR})"

    for SAMPLE in "${SAMPLES[@]}"; do

        r1=$(find_fastq "${RUN_DIR}" "${SAMPLE}" "R1")
        r2=$(find_fastq "${RUN_DIR}" "${SAMPLE}" "R2")

        if [[ -z "${r1}" || -z "${r2}" ]]; then
            echo "    [SKIP] ${SAMPLE} non présent dans ${RUN_LABEL}"
            continue
        fi

        # Répertoire de sortie FastQC pour ce run × échantillon
        OUTDIR_FASTQC="${FASTQC_RAW_DIR}/${RUN_LABEL}/${SAMPLE}"
        mkdir -p "${OUTDIR_FASTQC}"

        echo "    [FastQC] ${SAMPLE} / ${RUN_LABEL}"
        echo "      R1 : ${r1}"
        echo "      R2 : ${r2}"

        # Un appel FastQC par fichier pour avoir un rapport par séquence
        fastqc --threads "${THREADS}" --outdir "${OUTDIR_FASTQC}" "${r1}"
        fastqc --threads "${THREADS}" --outdir "${OUTDIR_FASTQC}" "${r2}"

        echo "      Done: FastQC terminé pour ${SAMPLE} / ${RUN_LABEL}"

    done
done

echo ""
echo " STEP 2 done. Rapports FastQC individuels dans : ${FASTQC_RAW_DIR}"

# ==============================================================================
# STEP 3 : MultiQC par échantillon × run
#
# Pour chaque combinaison échantillon/run qui a des rapports FastQC,
# on lance un MultiQC dans le sous-dossier correspondant.
# Exemple de structure produite :
#   multiqc_per_run/cop408_RUN1/multiqc_report.html
#   multiqc_per_run/cop408_RUN2/multiqc_report.html
#   ...
# ==============================================================================

echo "======================================================================="
echo " STEP 3: MultiQC par échantillon × run"
echo "======================================================================="

for i in "${!ALL_RUNS[@]}"; do
    RUN_LABEL="${RUN_LABELS[$i]}"

    for SAMPLE in "${SAMPLES[@]}"; do

        FASTQC_IN="${FASTQC_RAW_DIR}/${RUN_LABEL}/${SAMPLE}"

        # On vérifie qu'il y a des rapports FastQC pour ce couple
        if [[ ! -d "${FASTQC_IN}" ]] || \
           [[ $(find "${FASTQC_IN}" -maxdepth 1 -name "*_fastqc.zip" | wc -l) -eq 0 ]]; then
            echo "    [SKIP] Pas de rapports FastQC pour ${SAMPLE} / ${RUN_LABEL}"
            continue
        fi

        MULTIQC_OUT="${MULTIQC_RAW_DIR}/${SAMPLE}_${RUN_LABEL}"
        mkdir -p "${MULTIQC_OUT}"

        echo "    [MultiQC] ${SAMPLE} / ${RUN_LABEL} -> ${MULTIQC_OUT}"

        multiqc \
            "${FASTQC_IN}" \
            --outdir "${MULTIQC_OUT}" \
            --filename "multiqc_report_${SAMPLE}_${RUN_LABEL}" \
            --title "MultiQC — ${SAMPLE} ${RUN_LABEL}" \
            --force

        echo "      Done: MultiQC terminé pour ${SAMPLE} / ${RUN_LABEL}"

    done
done

echo ""
echo " STEP 3 done. Rapports MultiQC par run/échantillon dans : ${MULTIQC_RAW_DIR}"

# ==============================================================================
# STEP 4 : FastQC sur les fichiers groupés (concaténés)
# Un rapport par fichier groupé (R1 et R2 séparément).
# ==============================================================================

echo "======================================================================="
echo " STEP 4: FastQC sur les fichiers groupés (concaténés)"
echo "======================================================================="

mapfile -t GROUPED_FILES < <(
    find "${GROUPED_DIR}" -maxdepth 1 -type f -name "*.fastq.gz" | sort
)

if [[ ${#GROUPED_FILES[@]} -eq 0 ]]; then
    echo " ERROR: Aucun fichier groupé trouvé dans ${GROUPED_DIR}."
    echo " Vérifier que STEP 1 s'est bien terminé."
    exit 1
fi

echo "  Trouvé ${#GROUPED_FILES[@]} fichier(s) groupé(s). Lancement FastQC..."

for f in "${GROUPED_FILES[@]}"; do
    echo "    [FastQC grouped] ${f}"
    fastqc --threads "${THREADS}" --outdir "${FASTQC_GROUPED_DIR}" "${f}"
done

echo ""
echo " STEP 4 done. FastQC groupés dans : ${FASTQC_GROUPED_DIR}"

# ==============================================================================
# STEP 5 : MultiQC sur les fichiers groupés
# ==============================================================================

echo "======================================================================="
echo " STEP 5: MultiQC sur les fichiers groupés (concaténés)"
echo "======================================================================="

multiqc \
    "${FASTQC_GROUPED_DIR}" \
    --outdir "${MULTIQC_GROUPED_DIR}" \
    --filename "multiqc_report_grouped" \
    --title "MultiQC — Fichiers groupés (toutes données concaténées)" \
    --force

echo " STEP 5 done. MultiQC groupé dans : ${MULTIQC_GROUPED_DIR}"

# ==============================================================================
# RÉSUMÉ FINAL
# ==============================================================================

echo ""
echo "======================================================================="
echo " ALL STEPS COMPLETED SUCCESSFULLY"
echo "======================================================================="
echo ""
echo "  Fichiers groupés         : ${GROUPED_DIR}"
echo "  FastQC par run/séquence  : ${FASTQC_RAW_DIR}"
echo "  FastQC groupés           : ${FASTQC_GROUPED_DIR}"
echo "  MultiQC par run/échant.  : ${MULTIQC_RAW_DIR}"
echo "  MultiQC groupés          : ${MULTIQC_GROUPED_DIR}/multiqc_report_grouped.html"
echo "======================================================================="
