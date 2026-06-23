#!/bin/bash
#SBATCH --job-name=01_group_cop_files
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/01_group_cop_files.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/01_group_cop_files.out"

# ==============================================================================
# Script: 01_group_cop_files.sh
# Author: Pierre-Louis Stenger
# Purpose:
#   1. Collect all *cop*.fastq.gz files from multiple sequencing runs
#   2. Concatenate (cat) R1 and R2 reads per sample across all runs
#   3. Run FastQC on raw (per-run) files
#   4. Run FastQC on the concatenated (grouped) files
#   5. Run MultiQC on both raw and grouped FastQC outputs
# ==============================================================================

# set -euo pipefail

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

module load conda/4.12.0
source ~/.bashrc
conda activate bioinformatic

# ==============================================================================
# PATHS — Run directories
# ==============================================================================

# RUN1: flat structure — files directly in the folder
#   e.g. cop408_R1.fastq.gz
RUN1="/home/plstenge/coprolites_comparison/01_raw_data/Lot1_Illumina_R1"

# RUN2: subdirectory structure — files inside 474_cop408/, 476_cop412/, etc.
#   e.g. 474_cop408/474_cop408_R1.fastq.gz
RUN2="/storage/groups/gdec/shared_paleo/E1531_final/run1_20250320_AV241601_E1531_Ps5Lane1_Ps6Lane2"

# RUN3: subdirectory structure (same samples as RUN2)
RUN3="/storage/groups/gdec/shared_paleo/E1531_final/run2_20250414_AV241601_E1531_Ps5_Ps6_11022026_CORRECTED"

# RUN4: subdirectory structure (all 5 samples)
RUN4="/storage/groups/gdec/shared_paleo/E1531_final/run3_20251008_AV241601_E1531_Ps7_Ps8"

# RUN5: subdirectory structure (all 5 samples)
RUN5="/storage/groups/gdec/shared_paleo/E1531_final/run4_20251104_AV241601_E1531_Ps7_Ps8_04112025"

# ==============================================================================
# OUTPUT PATHS
# ==============================================================================

WORKDIR="/home/plstenge/coprolites_comparison"

# Directory for grouped (concatenated) FASTQ files
GROUPED_DIR="${WORKDIR}/01_raw_data/grouped"

# FastQC output on raw per-run FASTQ files (before grouping)
FASTQC_RAW_DIR="${WORKDIR}/02_quality_check_raw/fastqc_per_run"

# FastQC output on grouped/concatenated FASTQ files
FASTQC_GROUPED_DIR="${WORKDIR}/02_quality_check_raw/fastqc_grouped"

# MultiQC output on raw FastQC results
MULTIQC_RAW_DIR="${WORKDIR}/02_quality_check_raw/multiqc_per_run"

# MultiQC output on grouped FastQC results
MULTIQC_GROUPED_DIR="${WORKDIR}/02_quality_check_raw/multiqc_grouped"

# Number of threads for FastQC
THREADS=8

# ==============================================================================
# CREATE OUTPUT DIRECTORIES
# ==============================================================================

mkdir -p "${GROUPED_DIR}"
mkdir -p "${FASTQC_RAW_DIR}"
mkdir -p "${FASTQC_GROUPED_DIR}"
mkdir -p "${MULTIQC_RAW_DIR}"
mkdir -p "${MULTIQC_GROUPED_DIR}"

echo "======================================================================="
echo "  STEP 0: Directory structure created"
echo "======================================================================="

# ==============================================================================
# HELPER FUNCTION: find a FASTQ file for a given sample and read (R1 or R2)
#   in a given run directory, handling both flat and subdirectory structures.
#
# Usage: find_fastq <run_dir> <sample_id> <R1|R2>
# Returns: full path to the matching file, or empty string if not found
# ==============================================================================

find_fastq() {
    local run_dir="$1"
    local sample="$2"   # e.g. "cop408"
    local read="$3"     # "R1" or "R2"

    # Search for files matching *<sample>*<read>*.fastq.gz in the run directory,
    # at depth 1 (flat) or depth 2 (one subdirectory level deep).
    # Using find -maxdepth 2 covers both structures.
    find "${run_dir}" -maxdepth 2 -type f \
        -name "*${sample}*${read}*.fastq.gz" \
        ! -name "*.json" \
        2>/dev/null | head -n 1
}

# ==============================================================================
# STEP 1: RUN FastQC ON ALL RAW PER-RUN FILES
# ==============================================================================

echo "======================================================================="
echo "  STEP 1: FastQC on all raw per-run FASTQ files"
echo "======================================================================="

# Collect all cop FASTQ files from all runs and run FastQC on them.
# We use find with -maxdepth 2 to handle both flat and subdirectory structures.

for RUN_DIR in "${RUN1}" "${RUN2}" "${RUN3}" "${RUN4}" "${RUN5}"; do
    echo ""
    echo "  Processing run: ${RUN_DIR}"

    # Find all cop FASTQ files in this run directory (up to 2 levels deep)
    while IFS= read -r -d '' fastq_file; do
        echo "    FastQC (raw): ${fastq_file}"
        fastqc \
            --threads "${THREADS}" \
            --outdir "${FASTQC_RAW_DIR}" \
            "${fastq_file}"
    done < <(find "${RUN_DIR}" -maxdepth 2 -type f \
                -name "*cop*.fastq.gz" \
                ! -name "*.json" \
                -print0 2>/dev/null)
done

echo ""
echo "  STEP 1 done: FastQC raw outputs in ${FASTQC_RAW_DIR}"

# ==============================================================================
# STEP 2: CONCATENATE (cat) PER SAMPLE ACROSS ALL RUNS
# ==============================================================================

echo "======================================================================="
echo "  STEP 2: Concatenating FASTQ files per sample across all runs"
echo "======================================================================="

# List of the 5 samples to process
SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")

for SAMPLE in "${SAMPLES[@]}"; do

    echo ""
    echo "  Sample: ${SAMPLE}"

    # Collect all R1 files across all runs for this sample, in order
    R1_FILES=()
    R2_FILES=()

    for RUN_DIR in "${RUN1}" "${RUN2}" "${RUN3}" "${RUN4}" "${RUN5}"; do
        r1=$(find_fastq "${RUN_DIR}" "${SAMPLE}" "R1")
        r2=$(find_fastq "${RUN_DIR}" "${SAMPLE}" "R2")

        if [[ -n "${r1}" && -n "${r2}" ]]; then
            echo "    Found in: ${RUN_DIR}"
            echo "      R1: ${r1}"
            echo "      R2: ${r2}"
            R1_FILES+=("${r1}")
            R2_FILES+=("${r2}")
        else
            echo "    Not found in: ${RUN_DIR} (skipping)"
        fi
    done

    # Verify that we found at least one run for this sample
    if [[ ${#R1_FILES[@]} -eq 0 ]]; then
        echo "  WARNING: No files found for sample ${SAMPLE}. Skipping."
        continue
    fi

    # Define output paths for concatenated files
    OUT_R1="${GROUPED_DIR}/${SAMPLE}_grouped_R1.fastq.gz"
    OUT_R2="${GROUPED_DIR}/${SAMPLE}_grouped_R2.fastq.gz"

    # Concatenate all R1 files for this sample
    echo "    Concatenating ${#R1_FILES[@]} R1 file(s) -> ${OUT_R1}"
    cat "${R1_FILES[@]}" > "${OUT_R1}"

    # Concatenate all R2 files for this sample
    echo "    Concatenating ${#R2_FILES[@]} R2 file(s) -> ${OUT_R2}"
    cat "${R2_FILES[@]}" > "${OUT_R2}"

    echo "    Done: ${SAMPLE} grouped files written."

done

echo ""
echo "  STEP 2 done: Grouped FASTQ files in ${GROUPED_DIR}"

# ==============================================================================
# STEP 3: RUN FastQC ON GROUPED (CONCATENATED) FILES
# ==============================================================================

echo "======================================================================="
echo "  STEP 3: FastQC on grouped (concatenated) FASTQ files"
echo "======================================================================="

for GROUPED_FILE in "${GROUPED_DIR}"/*.fastq.gz; do
    echo "  FastQC (grouped): ${GROUPED_FILE}"
    fastqc \
        --threads "${THREADS}" \
        --outdir "${FASTQC_GROUPED_DIR}" \
        "${GROUPED_FILE}"
done

echo ""
echo "  STEP 3 done: FastQC grouped outputs in ${FASTQC_GROUPED_DIR}"

# ==============================================================================
# STEP 4: RUN MultiQC ON RAW FastQC OUTPUTS
# ==============================================================================

echo "======================================================================="
echo "  STEP 4: MultiQC on raw (per-run) FastQC outputs"
echo "======================================================================="

multiqc \
    "${FASTQC_RAW_DIR}" \
    --outdir "${MULTIQC_RAW_DIR}" \
    --filename "multiqc_report_raw" \
    --title "MultiQC — Raw per-run FASTQ quality" \
    --force

echo "  STEP 4 done: MultiQC raw report in ${MULTIQC_RAW_DIR}"

# ==============================================================================
# STEP 5: RUN MultiQC ON GROUPED FastQC OUTPUTS
# ==============================================================================

echo "======================================================================="
echo "  STEP 5: MultiQC on grouped (concatenated) FastQC outputs"
echo "======================================================================="

multiqc \
    "${FASTQC_GROUPED_DIR}" \
    --outdir "${MULTIQC_GROUPED_DIR}" \
    --filename "multiqc_report_grouped" \
    --title "MultiQC — Grouped concatenated FASTQ quality" \
    --force

echo "  STEP 5 done: MultiQC grouped report in ${MULTIQC_GROUPED_DIR}"

# ==============================================================================
# SUMMARY
# ==============================================================================

echo ""
echo "======================================================================="
echo "  ALL STEPS COMPLETED SUCCESSFULLY"
echo "======================================================================="
echo ""
echo "  Results summary:"
echo "  -----------------------------------------------------------------------"
echo "  Grouped FASTQ files   : ${GROUPED_DIR}"
echo "  FastQC raw outputs    : ${FASTQC_RAW_DIR}"
echo "  FastQC grouped outputs: ${FASTQC_GROUPED_DIR}"
echo "  MultiQC raw report    : ${MULTIQC_RAW_DIR}/multiqc_report_raw.html"
echo "  MultiQC grouped report: ${MULTIQC_GROUPED_DIR}/multiqc_report_grouped.html"
echo "======================================================================="
