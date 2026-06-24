#!/bin/bash

#SBATCH --job-name=01_group_cop_files
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/01_group_cop_files.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/01_group_cop_files.out"

# ==============================================================================
# Script  : 01_group_cop_fastqc.sh
# Author  : Pierre-Louis Stenger
# Purpose :
#   1. Concatenate (cat) R1 and R2 reads per sample across all sequencing runs
#   2. Run FastQC on the original RUN1 FASTQ files (before grouping)
#   3. Run FastQC on each per-run file (RUN2–RUN5)
#   4. Run FastQC on the grouped (concatenated) files
#   5. Run MultiQC on raw FastQC outputs
#   6. Run MultiQC on grouped FastQC outputs
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

# RUN1: flat structure — files directly in the folder
#   Pattern: cop408_R1.fastq.gz
RUN1="/home/plstenge/coprolites_comparison/01_raw_data/Lot1_Illumina_R1"

# RUN2–RUN5: subdirectory structure — files inside subfolders like 474_cop408/
#   Pattern: 474_cop408/474_cop408_R1.fastq.gz
RUN2="/storage/groups/gdec/shared_paleo/E1531_final/run1_20250320_AV241601_E1531_Ps5Lane1_Ps6Lane2"
RUN3="/storage/groups/gdec/shared_paleo/E1531_final/run2_20250414_AV241601_E1531_Ps5_Ps6_11022026_CORRECTED"
RUN4="/storage/groups/gdec/shared_paleo/E1531_final/run3_20251008_AV241601_E1531_Ps7_Ps8"
RUN5="/storage/groups/gdec/shared_paleo/E1531_final/run4_20251104_AV241601_E1531_Ps7_Ps8_04112025"

# All run directories as an array (used for looping)
ALL_RUNS=("${RUN1}" "${RUN2}" "${RUN3}" "${RUN4}" "${RUN5}")

# ==============================================================================
# OUTPUT PATHS
# ==============================================================================

WORKDIR="/home/plstenge/coprolites_comparison"

# Directory where concatenated grouped FASTQ files will be written
GROUPED_DIR="${WORKDIR}/01_raw_data/grouped"

# FastQC output directory for raw per-run files
FASTQC_RAW_DIR="${WORKDIR}/02_quality_check_raw/fastqc_per_run"

# FastQC output directory for grouped/concatenated files
FASTQC_GROUPED_DIR="${WORKDIR}/02_quality_check_raw/fastqc_grouped"

# MultiQC output directory for raw FastQC results
MULTIQC_RAW_DIR="${WORKDIR}/02_quality_check_raw/multiqc_per_run"

# MultiQC output directory for grouped FastQC results
MULTIQC_GROUPED_DIR="${WORKDIR}/02_quality_check_raw/multiqc_grouped"

# Number of CPU threads available to FastQC
# Must match --cpus-per-task in the SLURM header
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ==============================================================================
# SAMPLES
# ==============================================================================

SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")

# ==============================================================================
# HELPER FUNCTION
# Finds the FASTQ file for a given sample and read direction (R1 or R2)
# within a run directory, handling both flat and one-level subdirectory
# structures (maxdepth 2).
#
# Usage  : find_fastq <run_dir> <sample_id> <R1|R2>
# Returns: full path to the matching file, or empty string if not found
# ==============================================================================

find_fastq() {
    local run_dir="$1"
    local sample="$2"   # e.g. "cop408"
    local read="$3"     # "R1" or "R2"

    find "${run_dir}" -maxdepth 2 -type f \
        -name "*${sample}*${read}*.fastq.gz" \
        2>/dev/null | head -n 1
}

# NOTE: The closing brace above is mandatory.
# Without it, bash treats everything below as part of the function body
# and none of the steps below would ever execute.

# ==============================================================================
# STEP 0: Create output directories
# ==============================================================================

echo "======================================================================="
echo "  STEP 0: Creating output directories"
echo "======================================================================="

mkdir -p "${GROUPED_DIR}"
mkdir -p "${FASTQC_RAW_DIR}"
mkdir -p "${FASTQC_GROUPED_DIR}"
mkdir -p "${MULTIQC_RAW_DIR}"
mkdir -p "${MULTIQC_GROUPED_DIR}"

echo "  Done."

# ==============================================================================
# STEP 1: Concatenate FASTQ files per sample across all runs
# ==============================================================================

echo "======================================================================="
echo "  STEP 1: Concatenating FASTQ files per sample across all runs"
echo "======================================================================="

for SAMPLE in "${SAMPLES[@]}"; do

    echo ""
    echo "  Sample: ${SAMPLE}"

    # Collect R1 and R2 file paths found in each run
    R1_FILES=()
    R2_FILES=()

    for RUN_DIR in "${ALL_RUNS[@]}"; do

        r1=$(find_fastq "${RUN_DIR}" "${SAMPLE}" "R1")
        r2=$(find_fastq "${RUN_DIR}" "${SAMPLE}" "R2")

        if [[ -n "${r1}" && -n "${r2}" ]]; then
            echo "    [FOUND] ${RUN_DIR}"
            echo "      R1 : ${r1}"
            echo "      R2 : ${r2}"
            R1_FILES+=("${r1}")
            R2_FILES+=("${r2}")
        else
            echo "    [SKIP]  ${RUN_DIR} — ${SAMPLE} not present"
        fi

    done

    # Skip the sample if no files were found in any run
    if [[ ${#R1_FILES[@]} -eq 0 ]]; then
        echo "  WARNING: ${SAMPLE} — no files found in any run. Skipping."
        continue
    fi

    # Define output file paths for the concatenated grouped files
    OUT_R1="${GROUPED_DIR}/${SAMPLE}_grouped_R1.fastq.gz"
    OUT_R2="${GROUPED_DIR}/${SAMPLE}_grouped_R2.fastq.gz"

    # Concatenate all R1 files for this sample into a single grouped file
    echo "    Concatenating ${#R1_FILES[@]} R1 file(s) -> ${OUT_R1}"
    cat "${R1_FILES[@]}" > "${OUT_R1}"

    # Concatenate all R2 files for this sample into a single grouped file
    echo "    Concatenating ${#R2_FILES[@]} R2 file(s) -> ${OUT_R2}"
    cat "${R2_FILES[@]}" > "${OUT_R2}"

    echo "    Done: ${SAMPLE} grouped files written."

done

echo ""
echo "  STEP 1 done. Grouped FASTQ files are in: ${GROUPED_DIR}"

# ==============================================================================
# STEP 2: FastQC on all raw per-run FASTQ files
# FastQC is run in batch mode: all matching files are passed at once to a
# single FastQC call per run, which is faster than one call per file.
# ==============================================================================

echo "======================================================================="
echo "  STEP 2: FastQC on all raw per-run FASTQ files"
echo "======================================================================="

for RUN_DIR in "${ALL_RUNS[@]}"; do

    echo ""
    echo "  Run: ${RUN_DIR}"

    # Build the list of cop FASTQ files for this run (flat or subdirectory)
    mapfile -t RAW_FILES < <(
        find "${RUN_DIR}" -maxdepth 2 -type f \
            -name "*cop*.fastq.gz" \
            2>/dev/null | sort
    )

    if [[ ${#RAW_FILES[@]} -eq 0 ]]; then
        echo "    No cop FASTQ files found — skipping."
        continue
    fi

    echo "    Found ${#RAW_FILES[@]} file(s). Launching FastQC..."

    # Pass all files at once; -t controls threads across files
    fastqc \
        --threads "${THREADS}" \
        --outdir "${FASTQC_RAW_DIR}" \
        "${RAW_FILES[@]}"

done

echo ""
echo "  STEP 2 done. FastQC raw outputs in: ${FASTQC_RAW_DIR}"

# ==============================================================================
# STEP 3: FastQC on grouped (concatenated) files
# All 10 grouped files (5 samples × R1/R2) are passed at once.
# ==============================================================================

echo "======================================================================="
echo "  STEP 3: FastQC on grouped (concatenated) FASTQ files"
echo "======================================================================="

mapfile -t GROUPED_FILES < <(
    find "${GROUPED_DIR}" -maxdepth 1 -type f -name "*.fastq.gz" | sort
)

if [[ ${#GROUPED_FILES[@]} -eq 0 ]]; then
    echo "  ERROR: No grouped FASTQ files found in ${GROUPED_DIR}."
    echo "  Check that STEP 1 completed correctly."
    exit 1
fi

echo "  Found ${#GROUPED_FILES[@]} grouped file(s). Launching FastQC..."

fastqc \
    --threads "${THREADS}" \
    --outdir "${FASTQC_GROUPED_DIR}" \
    "${GROUPED_FILES[@]}"

echo ""
echo "  STEP 3 done. FastQC grouped outputs in: ${FASTQC_GROUPED_DIR}"

# ==============================================================================
# STEP 4: MultiQC on raw FastQC outputs
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

echo "  STEP 4 done. MultiQC raw report in: ${MULTIQC_RAW_DIR}"

# ==============================================================================
# STEP 5: MultiQC on grouped FastQC outputs
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

echo "  STEP 5 done. MultiQC grouped report in: ${MULTIQC_GROUPED_DIR}"

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

echo ""
echo "======================================================================="
echo "  ALL STEPS COMPLETED SUCCESSFULLY"
echo "======================================================================="
echo ""
echo "  Grouped FASTQ files    : ${GROUPED_DIR}"
echo "  FastQC raw outputs     : ${FASTQC_RAW_DIR}"
echo "  FastQC grouped outputs : ${FASTQC_GROUPED_DIR}"
echo "  MultiQC raw report     : ${MULTIQC_RAW_DIR}/multiqc_report_raw.html"
echo "  MultiQC grouped report : ${MULTIQC_GROUPED_DIR}/multiqc_report_grouped.html"
echo "======================================================================="
