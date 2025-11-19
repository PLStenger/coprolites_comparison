#!/bin/bash

#SBATCH --job-name=coprolites_5lots
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites/00_scripts/coprolites_5lots.err"
#SBATCH --output="/home/plstenge/coprolites/00_scripts/coprolites_5lots.out"

################################################################################
# Pipeline d'analyse aDNA - Projet Coprolites - 5 LOTS
# Comparaison complète: 2 recettes × 2 technologies × démultiplexage × 2 longueurs
# Author: Pierre-Louis Stenger
# Date: November 2025
#
# Description:
# Ce pipeline analyse 5 lots de séquençage de coprolites:
#
# LOT 1: Recette1 - Illumina (5 échantillons) -- /home/plstenge/coprolites/01_raw_data/Illumina_*
# LOT 2: Recette1 - Aviti sans démultiplexage (3 échantillons) -- /home/plstenge/coprolites/01_raw_data/Aviti_*
# LOT 3: Recette1 - Aviti avec démultiplexage - 2 runs concaténés (3 échantillons)
#     /storage/groups/gdec/shared_paleo/E1531_final/run1_20250320_AV241601_E1531_Ps5Lane1_Ps6Lane2/[dossier]
#     /storage/groups/gdec/shared_paleo/E1531_final/run2_20250414_AV241601_E1531_Ps5_Ps6_14042025/[dossier]
# LOT 4: Recette2 - Aviti fragments courts 2x150bp (5 échantillons) -- /storage/groups/gdec/shared_paleo/E1531_final/run3_20251008_AV241601_E1531_Ps7_Ps8/[dossier]
# LOT 5: Recette2 - Aviti fragments courts 2x75bp  (5 échantillons) -- /storage/groups/gdec/shared_paleo/E1531_final/run4_20251104_AV241601_E1531_Ps7_Ps8_04112025/[dossier]
#
# Échantillons: cop408, cop410, cop412, cop414, cop417
################################################################################

set -eo pipefail

BASE_DIR="/home/plstenge/coprolites"
BBDUK="/home/plstenge/bbmap/bbduk.sh"
CLUMPIFY="/home/plstenge/bbmap/clumpify.sh"
PHIX="/home/plstenge/bbmap/resources/phix174_ill.ref.fa.gz"
KRAKEN2_DB="/home/plstenge/k2_core_nt_20250609"
KRAKENTOOLS_DIR="${BASE_DIR}/07_kraken2/KrakenTools"
THREADS=36

# Définition des lots et mapping
LOTS=("Lot1_Illumina_R1" "Lot2_Aviti_no_demux_R1" "Lot3_Aviti_demux_R1" "Lot4_Aviti_R2_150bp" "Lot5_Aviti_R2_75bp")
declare -A LOT_SAMPLES
LOT_SAMPLES["Lot1_Illumina_R1"]="cop408 cop410 cop412 cop414 cop417"
LOT_SAMPLES["Lot2_Aviti_no_demux_R1"]="cop408 cop412 cop414"
LOT_SAMPLES["Lot3_Aviti_demux_R1"]="cop408 cop412 cop414"
LOT_SAMPLES["Lot4_Aviti_R2_150bp"]="cop408 cop410 cop412 cop414 cop417"
LOT_SAMPLES["Lot5_Aviti_R2_75bp"]="cop408 cop410 cop412 cop414 cop417"

# Mapping dossiers pour demux + SmallFragment
declare -A RUN3_FOLDERS=( ["cop408"]=474_cop408 ["cop410"]=475_cop410 ["cop412"]=476_cop412 ["cop414"]=477_cop414 ["cop417"]=478_cop417 )
declare -A RUN4_FOLDERS=( ["cop408"]=474_cop408 ["cop410"]=475_cop410 ["cop412"]=476_cop412 ["cop414"]=477_cop414 ["cop417"]=478_cop417 )
declare -A DEMUX_FOLDERS=( ["cop408"]=474_cop408 ["cop412"]=476_cop412 ["cop414"]=477_cop414 )

# Chemins sources
RUN1="/storage/groups/gdec/shared_paleo/E1531_final/run1_20250320_AV241601_E1531_Ps5Lane1_Ps6Lane2"
RUN2="/storage/groups/gdec/shared_paleo/E1531_final/run2_20250414_AV241601_E1531_Ps5_Ps6_14042025"
RUN3="/storage/groups/gdec/shared_paleo/E1531_final/run3_20251008_AV241601_E1531_Ps7_Ps8"
RUN4="/storage/groups/gdec/shared_paleo/E1531_final/run4_20251104_AV241601_E1531_Ps7_Ps8_04112025"
RAW_HOME="/home/plstenge/coprolites/01_raw_data"

echo "=========================================="
echo "Pipeline aDNA - Coprolites - 5 LOTS"
echo "Date de début: $(date)"
echo "=========================================="

################################################################################
# ACTIVATION ENVIRONNEMENT
################################################################################
module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics

# Arborescence
for lot in "${LOTS[@]}"; do
  mkdir -p "${BASE_DIR}/01_raw_data/${lot}"
  mkdir -p "${BASE_DIR}/02_quality_check_raw/${lot}"
  mkdir -p "${BASE_DIR}/03_bbduk/${lot}"
  mkdir -p "${BASE_DIR}/04_fastuniq/${lot}"
  mkdir -p "${BASE_DIR}/05_clumpify/${lot}"
  mkdir -p "${BASE_DIR}/06_fastp/${lot}"
  mkdir -p "${BASE_DIR}/07_quality_check_clean/${lot}"
  mkdir -p "${BASE_DIR}/08_kraken2/${lot}"
  mkdir -p "${BASE_DIR}/09_krona/${lot}"
  mkdir -p "${BASE_DIR}/10_mpa_tables/${lot}"
done
mkdir -p "${BASE_DIR}/00_scripts"
mkdir -p "${BASE_DIR}/11_summary_tables"
################################################################################
# ORGANISATION DES DONNÉES
################################################################################
# -------- Lot 1 --------
echo ""; echo "--- LOT 1: Illumina Recette1 ---"
LOT1_DEST="${BASE_DIR}/01_raw_data/Lot1_Illumina_R1"
for sample in ${LOT_SAMPLES["Lot1_Illumina_R1"]}; do
  R1="${RAW_HOME}/Illumina_${sample}_R1.fastq.gz"
  R2="${RAW_HOME}/Illumina_${sample}_R2.fastq.gz"
  ln -sf "$R1" "${LOT1_DEST}/${sample}_R1.fastq.gz"
  ln -sf "$R2" "${LOT1_DEST}/${sample}_R2.fastq.gz"
done
# -------- Lot 2 --------
echo ""; echo "--- LOT 2: Aviti no demux Recette1 ---"
LOT2_DEST="${BASE_DIR}/01_raw_data/Lot2_Aviti_no_demux_R1"
for sample in ${LOT_SAMPLES["Lot2_Aviti_no_demux_R1"]}; do
  R1="${RAW_HOME}/Aviti_${sample}_R1.fastq.gz"
  R2="${RAW_HOME}/Aviti_${sample}_R2.fastq.gz"
  ln -sf "$R1" "${LOT2_DEST}/${sample}_R1.fastq.gz"
  ln -sf "$R2" "${LOT2_DEST}/${sample}_R2.fastq.gz"
done
# -------- Lot 3 --------
echo ""; echo "--- LOT 3: Aviti demux Recette1 ---"
LOT3_DEST="${BASE_DIR}/01_raw_data/Lot3_Aviti_demux_R1"
for sample in ${LOT_SAMPLES["Lot3_Aviti_demux_R1"]}; do
  fldr="${DEMUX_FOLDERS[$sample]}"
  R1a="$RUN1/${fldr}/${fldr}_R1.fastq.gz"
  R2a="$RUN1/${fldr}/${fldr}_R2.fastq.gz"
  R1b="$RUN2/${fldr}/${fldr}_R1.fastq.gz"
  R2b="$RUN2/${fldr}/${fldr}_R2.fastq.gz"
  outR1="${LOT3_DEST}/${sample}_R1.fastq.gz"
  outR2="${LOT3_DEST}/${sample}_R2.fastq.gz"
  if [[ -f "$R1a" ]] && [[ -f "$R1b" ]]; then cat "$R1a" "$R1b" > "$outR1"; fi
  if [[ -f "$R2a" ]] && [[ -f "$R2b" ]]; then cat "$R2a" "$R2b" > "$outR2"; fi
  if [[ ! -f "$outR1" ]] && [[ -f "$R1a" ]]; then ln -sf "$R1a" "$outR1"; elif [[ ! -f "$outR1" ]] && [[ -f "$R1b" ]]; then ln -sf "$R1b" "$outR1"; fi
  if [[ ! -f "$outR2" ]] && [[ -f "$R2a" ]]; then ln -sf "$R2a" "$outR2"; elif [[ ! -f "$outR2" ]] && [[ -f "$R2b" ]]; then ln -sf "$R2b" "$outR2"; fi
done
# -------- Lot 4 --------
echo ""; echo "--- LOT 4: Aviti Recette2 2x150bp ---"
LOT4_DEST="${BASE_DIR}/01_raw_data/Lot4_Aviti_R2_150bp"
for sample in ${LOT_SAMPLES["Lot4_Aviti_R2_150bp"]}; do
  fldr="${RUN3_FOLDERS[$sample]}"
  R1="$RUN3/${fldr}/${fldr}_R1.fastq.gz"
  R2="$RUN3/${fldr}/${fldr}_R2.fastq.gz"
  ln -sf "$R1" "${LOT4_DEST}/${sample}_R1.fastq.gz"
  ln -sf "$R2" "${LOT4_DEST}/${sample}_R2.fastq.gz"
done
# -------- Lot 5 --------
echo ""; echo "--- LOT 5: Aviti Recette2 2x75bp ---"
LOT5_DEST="${BASE_DIR}/01_raw_data/Lot5_Aviti_R2_75bp"
for sample in ${LOT_SAMPLES["Lot5_Aviti_R2_75bp"]}; do
  fldr="${RUN4_FOLDERS[$sample]}"
  R1="$RUN4/${fldr}/${fldr}_R1.fastq.gz"
  R2="$RUN4/${fldr}/${fldr}_R2.fastq.gz"
  ln -sf "$R1" "${LOT5_DEST}/${sample}_R1.fastq.gz"
  ln -sf "$R2" "${LOT5_DEST}/${sample}_R2.fastq.gz"
done
################################################################################
# PIPELINE COMMUN TOUTES LOTS
################################################################################

################################################################################
# ORGANISATION DES DONNÉES - LOT 1: Illumina Recette1
################################################################################

echo ""
echo "=== LOT 1: Organisation Illumina Recette1 ==="
echo ""

LOT1_SOURCE="/home/plstenge/coprolites/01_raw_data"
LOT1_DEST="${BASE_DIR}/01_raw_data/Lot1_Illumina_R1"

for sample in cop408 cop410 cop412 cop414 cop417; do
    R1_source="${LOT1_SOURCE}/Illumina_${sample}_R1.fastq.gz"
    R2_source="${LOT1_SOURCE}/Illumina_${sample}_R2.fastq.gz"
    
    R1_dest="${LOT1_DEST}/${sample}_R1.fastq.gz"
    R2_dest="${LOT1_DEST}/${sample}_R2.fastq.gz"
    
    if [[ -f "$R1_source" ]]; then
        ln -sf "$R1_source" "$R1_dest"
        echo " ✓ ${sample}_R1 lié (Lot1 Illumina)"
    else
        echo " ✗ ATTENTION: ${sample}_R1 non trouvé"
    fi
    
    if [[ -f "$R2_source" ]]; then
        ln -sf "$R2_source" "$R2_dest"
        echo " ✓ ${sample}_R2 lié (Lot1 Illumina)"
    else
        echo " ✗ ATTENTION: ${sample}_R2 non trouvé"
    fi
done

################################################################################
# ORGANISATION DES DONNÉES - LOT 2: Aviti no demux Recette1
################################################################################

echo ""
echo "=== LOT 2: Organisation Aviti no demux Recette1 ==="
echo ""

LOT2_SOURCE="/home/plstenge/coprolites/01_raw_data"
LOT2_DEST="${BASE_DIR}/01_raw_data/Lot2_Aviti_no_demux_R1"

for sample in cop408 cop412 cop414; do
    R1_source="${LOT2_SOURCE}/Aviti_${sample}_R1.fastq.gz"
    R2_source="${LOT2_SOURCE}/Aviti_${sample}_R2.fastq.gz"
    
    R1_dest="${LOT2_DEST}/${sample}_R1.fastq.gz"
    R2_dest="${LOT2_DEST}/${sample}_R2.fastq.gz"
    
    if [[ -f "$R1_source" ]]; then
        ln -sf "$R1_source" "$R1_dest"
        echo " ✓ ${sample}_R1 lié (Lot2 Aviti no demux)"
    else
        echo " ✗ ATTENTION: ${sample}_R1 non trouvé"
    fi
    
    if [[ -f "$R2_source" ]]; then
        ln -sf "$R2_source" "$R2_dest"
        echo " ✓ ${sample}_R2 lié (Lot2 Aviti no demux)"
    else
        echo " ✗ ATTENTION: ${sample}_R2 non trouvé"
    fi
done

################################################################################
# ORGANISATION DES DONNÉES - LOT 3: Aviti demux Recette1 (2 RUNS À CONCATÉNER)
################################################################################

echo ""
echo "=== LOT 3: Organisation Aviti demux Recette1 (concaténation 2 runs) ==="
echo ""

LOT3_DEST="${BASE_DIR}/01_raw_data/Lot3_Aviti_demux_R1"

# Chemins des 2 runs
RUN1_BASE="/storage/groups/gdec/shared_paleo/E1531_final/run1_20250320_AV241601_E1531_Ps5Lane1_Ps6Lane2"
RUN2_BASE="/storage/groups/gdec/shared_paleo/E1531_final/run2_20250414_AV241601_E1531_Ps5_Ps6_14042025"

# Mapping échantillons -> dossiers
declare -A LOT3_FOLDERS
LOT3_FOLDERS["cop408"]="474_cop408"
LOT3_FOLDERS["cop412"]="476_cop412"
LOT3_FOLDERS["cop414"]="477_cop414"

for sample in cop408 cop412 cop414; do
    folder="${LOT3_FOLDERS[$sample]}"
    
    # Fichiers run1
    run1_r1="${RUN1_BASE}/${folder}/${folder}_R1.fastq.gz"
    run1_r2="${RUN1_BASE}/${folder}/${folder}_R2.fastq.gz"
    
    # Fichiers run2
    run2_r1="${RUN2_BASE}/${folder}/${folder}_R1.fastq.gz"
    run2_r2="${RUN2_BASE}/${folder}/${folder}_R2.fastq.gz"
    
    # Destination finale
    dest_r1="${LOT3_DEST}/${sample}_R1.fastq.gz"
    dest_r2="${LOT3_DEST}/${sample}_R2.fastq.gz"
    
    echo " → Concaténation ${sample}..."
    
    # Concaténation R1
    if [[ -f "$run1_r1" ]] && [[ -f "$run2_r1" ]]; then
        cat "$run1_r1" "$run2_r1" > "$dest_r1"
        echo "   ✓ ${sample}_R1 concaténé (run1 + run2)"
    elif [[ -f "$run1_r1" ]]; then
        ln -sf "$run1_r1" "$dest_r1"
        echo "   → ${sample}_R1 run1 seul (run2 absent)"
    elif [[ -f "$run2_r1" ]]; then
        ln -sf "$run2_r1" "$dest_r1"
        echo "   → ${sample}_R1 run2 seul (run1 absent)"
    else
        echo "   ✗ ERREUR: Aucun fichier R1 trouvé pour ${sample}"
    fi
    
    # Concaténation R2
    if [[ -f "$run1_r2" ]] && [[ -f "$run2_r2" ]]; then
        cat "$run1_r2" "$run2_r2" > "$dest_r2"
        echo "   ✓ ${sample}_R2 concaténé (run1 + run2)"
    elif [[ -f "$run1_r2" ]]; then
        ln -sf "$run1_r2" "$dest_r2"
        echo "   → ${sample}_R2 run1 seul (run2 absent)"
    elif [[ -f "$run2_r2" ]]; then
        ln -sf "$run2_r2" "$dest_r2"
        echo "   → ${sample}_R2 run2 seul (run1 absent)"
    else
        echo "   ✗ ERREUR: Aucun fichier R2 trouvé pour ${sample}"
    fi
done

################################################################################
# ORGANISATION DES DONNÉES - LOT 4: Aviti Recette2
################################################################################

echo ""
echo "=== LOT 4: Organisation Aviti Recette2 ==="
echo ""

LOT4_SOURCE="/storage/groups/gdec/shared_paleo/E1531_final/run3_20251008_AV241601_E1531_Ps7_Ps8"
LOT4_DEST="${BASE_DIR}/01_raw_data/Lot4_Aviti_R2"

declare -A LOT4_FOLDERS
LOT4_FOLDERS["cop408"]="474_cop408"
LOT4_FOLDERS["cop410"]="475_cop410"
LOT4_FOLDERS["cop412"]="476_cop412"
LOT4_FOLDERS["cop414"]="477_cop414"
LOT4_FOLDERS["cop417"]="478_cop417"

for sample in cop408 cop410 cop412 cop414 cop417; do
    folder="${LOT4_FOLDERS[$sample]}"
    
    R1_source="${LOT4_SOURCE}/${folder}/${folder}_R1.fastq.gz"
    R2_source="${LOT4_SOURCE}/${folder}/${folder}_R2.fastq.gz"
    
    R1_dest="${LOT4_DEST}/${sample}_R1.fastq.gz"
    R2_dest="${LOT4_DEST}/${sample}_R2.fastq.gz"
    
    if [[ -f "$R1_source" ]]; then
        ln -sf "$R1_source" "$R1_dest"
        echo " ✓ ${sample}_R1 lié (Lot4 Aviti R2)"
    else
        echo " ✗ ATTENTION: ${sample}_R1 non trouvé dans ${folder}"
    fi
    
    if [[ -f "$R2_source" ]]; then
        ln -sf "$R2_source" "$R2_dest"
        echo " ✓ ${sample}_R2 lié (Lot4 Aviti R2)"
    else
        echo " ✗ ATTENTION: ${sample}_R2 non trouvé dans ${folder}"
    fi
done

echo "Organisation des données terminée pour les 4 lots."

################################################################################
# ÉTAPE 1: Contrôle qualité RAW avec FastQC et MultiQC
################################################################################

echo ""
echo "=== ÉTAPE 1: Contrôle qualité RAW (FastQC/MultiQC) ==="
echo ""

for lot in "${LOTS[@]}"; do
    echo "FastQC RAW pour ${lot}..."
    
    INPUT_DIR="${BASE_DIR}/01_raw_data/${lot}"
    OUTPUT_DIR="${BASE_DIR}/02_quality_check_raw/${lot}"
    
    for FILE in "${INPUT_DIR}"/*.fastq.gz; do
        if [[ -f "$FILE" ]]; then
            fastqc "$FILE" -o "$OUTPUT_DIR" -t 4
        fi
    done
    
    # MultiQC
    echo "MultiQC RAW pour ${lot}..."
    cd "$OUTPUT_DIR"
    multiqc . -n "multiqc_raw_${lot}.html" --force
done

echo "Contrôle qualité RAW terminé."

################################################################################
# ÉTAPE 2: Filtrage et trimming avec BBDuk
################################################################################

echo ""
echo "=== ÉTAPE 2: Filtrage et trimming (BBDuk) ==="
echo ""

for lot in "${LOTS[@]}"; do
    echo "BBDuk pour ${lot}..."
    
    INPUT_DIR="${BASE_DIR}/01_raw_data/${lot}"
    OUTPUT_DIR="${BASE_DIR}/03_bbduk/${lot}"
    
    cd "$INPUT_DIR"
    
    for r1_file in *_R1.fastq.gz; do
        r2_file="${r1_file/_R1/_R2}"
        
        if [[ ! -f "$r2_file" ]]; then
            echo " ✗ ERREUR: Fichier R2 manquant pour $r1_file" >&2
            continue
        fi
        
        base_name="${r1_file%%_R1.fastq.gz}"
        echo " → Traitement de ${base_name} (${lot})..."
        
        $BBDUK -Xmx4g \
            in1="$r1_file" \
            in2="$r2_file" \
            out1="${OUTPUT_DIR}/clean_${r1_file}" \
            out2="${OUTPUT_DIR}/clean_${r2_file}" \
            ref=$PHIX \
            ktrim=rl \
            k=23 \
            mink=11 \
            hdist=1 \
            tpe \
            tbo \
            minlen=25 \
            qtrim=r \
            trimq=20 \
            stats="${OUTPUT_DIR}/${base_name}_bbduk_stats.txt"
    done
done

echo "Filtrage BBDuk terminé."

################################################################################
# ÉTAPE 3: Déduplication avec FastUniq
################################################################################

echo ""
echo "=== ÉTAPE 3: Déduplication (FastUniq) ==="
echo ""

TMP="/tmp/fastuniq_coprolites_4lots_$$"
mkdir -p "$TMP"

for lot in "${LOTS[@]}"; do
    echo "FastUniq pour ${lot}..."
    
    INPUT_DIR="${BASE_DIR}/03_bbduk/${lot}"
    OUTPUT_DIR="${BASE_DIR}/04_fastuniq/${lot}"
    
    cd "$INPUT_DIR" || continue
    
    for R1_gz in clean_*_R1.fastq.gz; do
        base=$(echo "$R1_gz" | sed 's/_R1\.fastq\.gz//')
        R2_gz="${base}_R2.fastq.gz"
        
        if [[ -f "$R2_gz" ]]; then
            echo " → Traitement de ${base} (${lot})..."
            
            R1_tmp="${TMP}/${lot}_${base}_R1.fastq"
            R2_tmp="${TMP}/${lot}_${base}_R2.fastq"
            listfile="${TMP}/${lot}_${base}.list"
            
            # Décompression
            zcat "$INPUT_DIR/$R1_gz" > "$R1_tmp"
            zcat "$INPUT_DIR/$R2_gz" > "$R2_tmp"
            
            # Création du fichier liste
            echo -e "${R1_tmp}\n${R2_tmp}" > "$listfile"
            
            # FastUniq
            fastuniq -i "$listfile" -t q \
                -o "${OUTPUT_DIR}/${base}_dedup_R1.fastq" \
                -p "${OUTPUT_DIR}/${base}_dedup_R2.fastq"
            
            # Nettoyage
            rm -f "$R1_tmp" "$R2_tmp" "$listfile"
        else
            echo " ✗ ATTENTION: fichier R2 manquant pour $base"
        fi
    done
done

rm -rf "$TMP"

echo "Déduplication FastUniq terminée."

################################################################################
# ÉTAPE 4: Clumpify - Déduplication optique
################################################################################

echo ""
echo "=== ÉTAPE 4: Clumpify (déduplication optique) ==="
echo ""

for lot in "${LOTS[@]}"; do
    echo "Clumpify pour ${lot}..."
    
    INPUT_DIR="${BASE_DIR}/04_fastuniq/${lot}"
    OUTPUT_DIR="${BASE_DIR}/05_clumpify/${lot}"
    
    for R1 in "${INPUT_DIR}"/*_R1.fastq; do
        R2="${R1/_R1.fastq/_R2.fastq}"
        
        if [[ -f "$R2" ]]; then
            base=$(basename "$R1" _R1.fastq)
            echo " → Traitement de ${base} (${lot})..."
            
            $CLUMPIFY \
                in="$R1" in2="$R2" \
                out="${OUTPUT_DIR}/${base}_clumpify_R1.fastq.gz" \
                out2="${OUTPUT_DIR}/${base}_clumpify_R2.fastq.gz" \
                dedupe=t
        else
            echo " ✗ Fichier R2 manquant pour $R1"
        fi
    done
done

echo "Clumpify terminé."

################################################################################
# ÉTAPE 5: Fastp - Merging et contrôle qualité final
################################################################################

echo ""
echo "=== ÉTAPE 5: Fastp (merging et QC final) ==="
echo ""

for lot in "${LOTS[@]}"; do
    echo "Fastp pour ${lot}..."
    
    INPUT_DIR="${BASE_DIR}/05_clumpify/${lot}"
    OUTPUT_DIR="${BASE_DIR}/06_fastp/${lot}"
    
    for R1 in "${INPUT_DIR}"/*_R1.fastq.gz; do
        R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
        
        if [[ -f "$R2" ]]; then
            base=$(basename "$R1" _R1.fastq.gz)
            echo " → Traitement de ${base} (${lot})..."
            
            fastp \
                -i "$R1" \
                -I "$R2" \
                --merged_out "${OUTPUT_DIR}/${base}_fastp_merged.fastq.gz" \
                --out1 "${OUTPUT_DIR}/${base}_fastp_R1.fastq.gz" \
                --out2 "${OUTPUT_DIR}/${base}_fastp_R2.fastq.gz" \
                --json "${OUTPUT_DIR}/${base}_fastp.json" \
                --html "${OUTPUT_DIR}/${base}_fastp.html" \
                --thread 4 \
                --length_required 30 \
                --qualified_quality_phred 20
        else
            echo " ✗ Fichier R2 manquant pour $R1"
        fi
    done
done

echo "Fastp terminé."

################################################################################
# ÉTAPE 6: Contrôle qualité CLEAN avec FastQC et MultiQC
################################################################################

echo ""
echo "=== ÉTAPE 6: Contrôle qualité CLEAN (FastQC/MultiQC) ==="
echo ""

for lot in "${LOTS[@]}"; do
    echo "FastQC CLEAN pour ${lot}..."
    
    INPUT_DIR="${BASE_DIR}/06_fastp/${lot}"
    OUTPUT_DIR="${BASE_DIR}/07_quality_check_clean/${lot}"
    
    # FastQC sur tous les fichiers fastp (merged, R1, R2)
    for FILE in "${INPUT_DIR}"/*.fastq.gz; do
        if [[ -f "$FILE" ]]; then
            fastqc "$FILE" -o "$OUTPUT_DIR" -t 4
        fi
    done
    
    # MultiQC
    echo "MultiQC CLEAN pour ${lot}..."
    cd "$OUTPUT_DIR"
    multiqc . -n "multiqc_clean_${lot}.html" --force
done

echo "Contrôle qualité CLEAN terminé."

################################################################################
# ÉTAPE 7: Classification taxonomique avec Kraken2
################################################################################

echo ""
echo "=== ÉTAPE 7: Classification taxonomique (Kraken2) ==="
echo ""

for lot in "${LOTS[@]}"; do
    echo "Kraken2 pour ${lot}..."
    
    FASTP_DIR="${BASE_DIR}/06_fastp/${lot}"
    OUT_DIR="${BASE_DIR}/08_kraken2/${lot}"
    
    # Analyse des reads merged
    echo " → Analyse des reads merged..."
    for MERGED in "${FASTP_DIR}"/*_fastp_merged.fastq.gz; do
        if [[ -f "$MERGED" ]]; then
            SAMPLE=$(basename "$MERGED" _fastp_merged.fastq.gz)
            OUT_KRAKEN="${OUT_DIR}/${SAMPLE}_merged.kraken"
            OUT_REPORT="${OUT_DIR}/${SAMPLE}_merged.report"
            
            echo "   • ${SAMPLE} (merged) - ${lot}"
            
            kraken2 --confidence 0.2 --db "$KRAKEN2_DB" --threads $THREADS \
                --output "$OUT_KRAKEN" --report "$OUT_REPORT" "$MERGED"
        fi
    done
    
    # Analyse des reads unmerged
    echo " → Analyse des reads unmerged..."
    for R1 in "${FASTP_DIR}"/*_fastp_R1.fastq.gz; do
        if [[ -f "$R1" ]]; then
            SAMPLE=$(basename "$R1" _fastp_R1.fastq.gz)
            R2="${FASTP_DIR}/${SAMPLE}_fastp_R2.fastq.gz"
            
            if [[ -f "$R2" ]]; then
                OUT_KRAKEN="${OUT_DIR}/${SAMPLE}_unmerged.kraken"
                OUT_REPORT="${OUT_DIR}/${SAMPLE}_unmerged.report"
                
                echo "   • ${SAMPLE} (unmerged) - ${lot}"
                
                kraken2 --confidence 0.2 --paired --db "$KRAKEN2_DB" --threads $THREADS \
                    --output "$OUT_KRAKEN" --report "$OUT_REPORT" "$R1" "$R2"
            fi
        fi
    done
done

echo "Classification Kraken2 terminée."

################################################################################
# ÉTAPE 8: Visualisation avec Krona
################################################################################

echo ""
echo "=== ÉTAPE 8: Visualisation (Krona) ==="
echo ""

for lot in "${LOTS[@]}"; do
    echo "Krona pour ${lot}..."
    
    IN_DIR="${BASE_DIR}/08_kraken2/${lot}"
    OUT_DIR="${BASE_DIR}/09_krona/${lot}"
    
    cd "$IN_DIR"
    
    # Krona pour tous les échantillons combinés
    if ls *.report 1> /dev/null 2>&1; then
        ktImportTaxonomy -t 5 -m 3 -o "${OUT_DIR}/all_samples_krona.html" "${IN_DIR}"/*.report
    fi
    
    # Krona individuel
    for report in "${IN_DIR}"/*.report; do
        if [[ -f "$report" ]]; then
            base=$(basename "$report" .report)
            ktImportTaxonomy -t 5 -m 3 -o "${OUT_DIR}/${base}_krona.html" "$report"
        fi
    done
done

echo "Visualisation Krona terminée."

################################################################################
# ÉTAPE 9: Création des tables MPA
################################################################################

echo ""
echo "=== ÉTAPE 9: Création des tables MPA ==="
echo ""

# Installation KrakenTools si nécessaire
if [[ ! -d "$KRAKENTOOLS_DIR" ]]; then
    echo "Installation de KrakenTools..."
    mkdir -p "${BASE_DIR}/08_kraken2"
    cd "${BASE_DIR}/08_kraken2"
    git clone https://github.com/jenniferlu717/KrakenTools.git
fi

for lot in "${LOTS[@]}"; do
    echo "Conversion MPA pour ${lot}..."
    
    IN_DIR="${BASE_DIR}/08_kraken2/${lot}"
    OUT_DIR="${BASE_DIR}/10_mpa_tables/${lot}"
    
    cd "$IN_DIR"
    
    declare -a mpa_files=()
    
    for report in *.report; do
        if [[ -f "$report" ]]; then
            base=$(basename "$report" .report)
            mpa_file="${OUT_DIR}/${base}.mpa"
            
            echo " → Conversion de ${base}..."
            python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" -r "$report" -o "$mpa_file"
            
            mpa_files+=("$mpa_file")
        fi
    done
    
    # Combinaison
    if [[ ${#mpa_files[@]} -gt 0 ]]; then
        echo " → Combinaison de tous les fichiers MPA..."
        python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" -i "${mpa_files[@]}" -o "${OUT_DIR}/combined_all.tsv"
    fi
done

echo "Création des tables MPA terminée."

################################################################################
# ÉTAPE 10: Génération du tableau récapitulatif des séquences
################################################################################

echo ""
echo "=== ÉTAPE 10: Tableau récapitulatif (avant/après nettoyage) ==="
echo ""

SUMMARY_TABLE="${BASE_DIR}/11_summary_tables/sequences_summary.tsv"

# En-tête du tableau
cat > "$SUMMARY_TABLE" << 'HEADER'
Lot	Sample	Stage	Nb_sequences	Longueur_moyenne	GC_percent
HEADER

# Fonction pour extraire les stats d'un fichier fastq.gz
function extract_stats() {
    local file=$1
    local lot=$2
    local sample=$3
    local stage=$4
    
    if [[ ! -f "$file" ]]; then
        echo "FICHIER_ABSENT" >&2
        return 1
    fi
    
    # Compter le nombre de séquences
    nb_seq=$(zcat "$file" 2>/dev/null | echo $((`wc -l`/4)))
    
    # Calculer longueur moyenne et GC%
    stats=$(zcat "$file" 2>/dev/null | awk 'NR%4==2 {
        total_len += length($0)
        gc_count += gsub(/[GCgc]/, "", $0)
        at_count += gsub(/[ATat]/, "", $0)
        count++
    } END {
        if (count > 0) {
            avg_len = total_len / count
            gc_perc = (gc_count / (gc_count + at_count)) * 100
            printf "%.1f\t%.2f", avg_len, gc_perc
        } else {
            printf "0\t0"
        }
    }')
    
    echo -e "${lot}\t${sample}\t${stage}\t${nb_seq}\t${stats}" >> "$SUMMARY_TABLE"
}

echo "Calcul des statistiques pour chaque échantillon..."

for lot in "${LOTS[@]}"; do
    echo " → Traitement du ${lot}..."
    
    samples="${LOT_SAMPLES[$lot]}"
    
    for sample in $samples; do
        echo "   • ${sample}..."
        
        # RAW
        raw_r1="${BASE_DIR}/01_raw_data/${lot}/${sample}_R1.fastq.gz"
        extract_stats "$raw_r1" "$lot" "$sample" "RAW"
        
        # CLEAN (après Fastp - merged)
        clean_merged="${BASE_DIR}/06_fastp/${lot}/clean_${sample}_dedup_clumpify_fastp_merged.fastq.gz"
        if [[ -f "$clean_merged" ]]; then
            extract_stats "$clean_merged" "$lot" "$sample" "CLEAN_merged"
        fi
        
        # CLEAN (après Fastp - R1 unmerged)
        clean_r1="${BASE_DIR}/06_fastp/${lot}/clean_${sample}_dedup_clumpify_fastp_R1.fastq.gz"
        if [[ -f "$clean_r1" ]]; then
            extract_stats "$clean_r1" "$lot" "$sample" "CLEAN_unmerged"
        fi
    done
done

echo ""
echo "Tableau récapitulatif créé: ${SUMMARY_TABLE}"
echo ""
echo "Aperçu:"
head -20 "$SUMMARY_TABLE" | column -t

################################################################################
# ÉTAPE 11: Génération du rapport global de comparaison
################################################################################

echo ""
echo "=== ÉTAPE 11: Génération du rapport global ==="
echo ""

REPORT_FILE="${BASE_DIR}/00_scripts/pipeline_4lots_report.txt"

cat > "$REPORT_FILE" << 'EOF'
================================================================================
RAPPORT D'ANALYSE - PROJET COPROLITES - 4 LOTS
================================================================================
Date d'analyse: $(date)
Pipeline version: 3.0 (4 lots - environnement conda unique)

================================================================================
DESCRIPTION DU PROJET
================================================================================

Analyse comparative de 4 lots de séquençage de coprolites:

LOT 1: Recette1 - Illumina
  • Technologie: Illumina
  • Protocole: recipes-ShortInsert (2x150bp)
  • Échantillons: cop408, cop410, cop412, cop414, cop417 (5 échantillons)
  • Source: /home/plstenge/coprolites/01_raw_data

LOT 2: Recette1 - Aviti sans démultiplexage
  • Technologie: Aviti (no demultiplexing)
  • Protocole: recipes-ShortInsert (2x150bp)
  • Échantillons: cop408, cop412, cop414 (3 échantillons)
  • Source: /home/plstenge/coprolites/01_raw_data

LOT 3: Recette1 - Aviti avec démultiplexage (2 runs concaténés)
  • Technologie: Aviti (demultiplexing)
  • Protocole: recipes-ShortInsert (2x150bp)
  • Échantillons: cop408, cop412, cop414 (3 échantillons)
  • Sources:
    - Run1: /storage/groups/gdec/shared_paleo/E1531_final/run1_20250320_AV241601_E1531_Ps5Lane1_Ps6Lane2
    - Run2: /storage/groups/gdec/shared_paleo/E1531_final/run2_20250414_AV241601_E1531_Ps5_Ps6_14042025
  • IMPORTANT: Les séquences des 2 runs ont été concaténées

LOT 4: Recette2 - Aviti avec récupération de fragments courts
  • Technologie: Aviti
  • Protocole: recipes-v3.3.x-ShortInsert-SmallFragmentRecovery (2x150bp)
  • Échantillons: cop408, cop410, cop412, cop414, cop417 (5 échantillons)
  • Source: /storage/groups/gdec/shared_paleo/E1531_final/run3_20251008_AV241601_E1531_Ps7_Ps8
  • IMPORTANT: Optimisé pour l'ADN ancien (fragments courts)

================================================================================
ÉTAPES DU PIPELINE
================================================================================

1. Contrôle qualité RAW (FastQC/MultiQC)
   → Analyse de la qualité des données brutes

2. Filtrage et trimming (BBDuk)
   → Élimination des adaptateurs et phiX
   → Trimming qualité (Q≥20)
   → Longueur minimale: 25 bp

3. Déduplication (FastUniq)
   → Suppression des duplicats PCR

4. Clumpify
   → Déduplication optique supplémentaire

5. Fastp
   → Merging des reads paired-end
   → Contrôle qualité final
   → Longueur minimale: 30 bp

6. Contrôle qualité CLEAN (FastQC/MultiQC)
   → Analyse de la qualité après nettoyage complet

7. Classification taxonomique (Kraken2)
   → Base de données: nt (complete)
   → Confidence: 0.2
   → Analyse merged + unmerged

8. Visualisation (Krona)
   → Charts interactifs de composition taxonomique

9. Tables d'assignation (MPA format)
   → Pour analyses statistiques downstream

10. Tableau récapitulatif
    → Nombre de séquences, longueur moyenne, %GC
    → Avant et après nettoyage

IMPORTANT: Tous les outils exécutés depuis un environnement conda unique.

================================================================================
LOCALISATION DES RÉSULTATS
================================================================================

Répertoire principal: /home/plstenge/coprolites

Structure (pour chaque lot):

├── 01_raw_data/{Lot}/                  # Données brutes organisées
├── 02_quality_check_raw/{Lot}/         # QC avant nettoyage
├── 03_bbduk/{Lot}/                     # Reads filtrés
├── 04_fastuniq/{Lot}/                  # Reads dédupliqués
├── 05_clumpify/{Lot}/                  # Déduplication optique
├── 06_fastp/{Lot}/                     # Reads merged/unmerged finaux
├── 07_quality_check_clean/{Lot}/       # QC après nettoyage
├── 08_kraken2/{Lot}/                   # Classification taxonomique
├── 09_krona/{Lot}/                     # Visualisations interactives
├── 10_mpa_tables/{Lot}/                # Tables d'assignation
└── 11_summary_tables/                  # Tableaux récapitulatifs

Les 4 lots:
- Lot1_Illumina_R1
- Lot2_Aviti_no_demux_R1
- Lot3_Aviti_demux_R1
- Lot4_Aviti_R2

================================================================================
TABLEAU RÉCAPITULATIF DES SÉQUENCES
================================================================================

Disponible dans: 11_summary_tables/sequences_summary.tsv

Colonnes:
- Lot: Identifiant du lot
- Sample: Nom de l'échantillon
- Stage: RAW / CLEAN_merged / CLEAN_unmerged
- Nb_sequences: Nombre de séquences
- Longueur_moyenne: Longueur moyenne des reads (bp)
- GC_percent: Pourcentage GC

Ce tableau permet de comparer:
1. La perte de séquences à chaque étape de nettoyage
2. Les différences de longueur entre lots et technologies
3. Les variations de %GC (indicateur de contamination)
4. L'efficacité du merging (aDNA = taux de merging élevé)

================================================================================
ANALYSES DE COMPARAISON RECOMMANDÉES
================================================================================

1. Comparaison des technologies (Illumina vs Aviti):
   - LOT 1 (Illumina) vs LOT 2/3/4 (Aviti)
   - Rendement, qualité, biais taxonomiques?

2. Impact du démultiplexage (Aviti):
   - LOT 2 (no demux) vs LOT 3 (demux)
   - Différences de qualité et rendement?

3. Comparaison des recettes (R1 vs R2):
   - LOT 1/2/3 (Recette1) vs LOT 4 (Recette2)
   - Efficacité de récupération des fragments courts?

4. Effet de la concaténation (LOT 3):
   - Augmentation de profondeur attendue
   - Cohérence entre les 2 runs?

5. Composition taxonomique:
   - Utiliser les fichiers Krona pour visualisation
   - Comparer la diversité entre lots
   - Identifier les biais technologiques

================================================================================
COMMANDES R POUR ANALYSES COMPLÉMENTAIRES
================================================================================

# Chargement du tableau récapitulatif
library(tidyverse)

summary <- read_tsv("/home/plstenge/coprolites/11_summary_tables/sequences_summary.tsv")

# Visualisation du nombre de séquences par lot
summary %>%
  filter(Stage == "RAW") %>%
  ggplot(aes(x = Sample, y = Nb_sequences, fill = Lot)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Nombre de séquences RAW par lot et échantillon",
       y = "Nombre de séquences") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Taux de rétention après nettoyage
retention <- summary %>%
  pivot_wider(names_from = Stage, values_from = Nb_sequences) %>%
  mutate(Retention_pct = (CLEAN_merged + CLEAN_unmerged) / RAW * 100)

ggplot(retention, aes(x = Lot, y = Retention_pct, fill = Lot)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Taux de rétention des séquences après nettoyage",
       y = "Rétention (%)")

# Comparaison des longueurs moyennes
summary %>%
  ggplot(aes(x = Lot, y = Longueur_moyenne, fill = Stage)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribution des longueurs moyennes",
       y = "Longueur moyenne (bp)")

# Comparaison %GC
summary %>%
  ggplot(aes(x = Lot, y = GC_percent, fill = Stage)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribution du %GC",
       y = "%GC")

# Chargement et comparaison des tables MPA
mpa_lot1 <- read_tsv("/home/plstenge/coprolites/10_mpa_tables/Lot1_Illumina_R1/combined_all.tsv")
mpa_lot2 <- read_tsv("/home/plstenge/coprolites/10_mpa_tables/Lot2_Aviti_no_demux_R1/combined_all.tsv")
mpa_lot3 <- read_tsv("/home/plstenge/coprolites/10_mpa_tables/Lot3_Aviti_demux_R1/combined_all.tsv")
mpa_lot4 <- read_tsv("/home/plstenge/coprolites/10_mpa_tables/Lot4_Aviti_R2/combined_all.tsv")

# Richesse taxonomique par lot
richness <- tibble(
  Lot = c("Lot1_Illumina", "Lot2_Aviti_no_demux", "Lot3_Aviti_demux", "Lot4_Aviti_R2"),
  N_taxa = c(nrow(mpa_lot1), nrow(mpa_lot2), nrow(mpa_lot3), nrow(mpa_lot4))
)

ggplot(richness, aes(x = Lot, y = N_taxa, fill = Lot)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Richesse taxonomique par lot",
       y = "Nombre de taxons détectés") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Diagramme de Venn (taxons partagés/uniques)
library(ggvenn)

taxa_lists <- list(
  Lot1_Illumina = mpa_lot1$taxonomy,
  Lot2_Aviti_no_demux = mpa_lot2$taxonomy,
  Lot3_Aviti_demux = mpa_lot3$taxonomy,
  Lot4_Aviti_R2 = mpa_lot4$taxonomy
)

ggvenn(taxa_lists)

================================================================================
INFORMATIONS TECHNIQUES
================================================================================

Base de données Kraken2: /home/plstenge/k2_core_nt_20250609
KrakenTools: /home/plstenge/coprolites/08_kraken2/KrakenTools
Environnement conda: coprolites-pipeline

Threads utilisés: 36
Mémoire allouée: 1000G

Outils et versions:
- FastQC
- MultiQC
- BBMap (BBDuk, Clumpify)
- FastUniq
- Fastp
- Kraken2
- Krona
- KrakenTools

================================================================================
CONTACT
================================================================================

Pour toute question: pierrelouis.stenger@gmail.com

================================================================================
EOF

# Substitution de la date
sed -i "s/\$(date)/$(date)/" "$REPORT_FILE"

cat "$REPORT_FILE"

echo ""
echo "Rapport généré: $REPORT_FILE"

################################################################################
# FIN DU PIPELINE
################################################################################

echo ""
echo "=========================================="
echo "PIPELINE TERMINÉ AVEC SUCCÈS - 5 LOTS"
echo "Date de fin: $(date)"
echo "=========================================="
echo ""
echo "Résultats disponibles dans: ${BASE_DIR}"
echo ""
echo "Fichiers importants:"
echo " • Rapport: ${REPORT_FILE}"
echo " • Tableau récapitulatif: ${SUMMARY_TABLE}"
echo ""
echo "Prochaines étapes:"
echo " 1. Consultez les rapports MultiQC (RAW et CLEAN)"
echo " 2. Analysez le tableau récapitulatif des séquences"
echo " 3. Explorez les visualisations Krona"
echo " 4. Comparez les 5 lots avec R"
echo ""


# Désactivation conda
conda deactivate

# Notification email
mail -s "Pipeline Coprolites 5 lots terminé" pierrelouis.stenger@gmail.com <<< "Pipeline terminé le $(date). Résultats: ${BASE_DIR}"

exit 0
