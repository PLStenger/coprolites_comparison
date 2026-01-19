#!/bin/bash

#SBATCH --job-name=coprolites_merged
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/coprolites_merged.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/coprolites_merged.out"

################################################################################
# Pipeline d'analyse aDNA - Projet Coprolites - MERGED RUNS
# Author: Pierre-Louis Stenger
# Date: Janvier 2026
# 
# Ce script merge les reads de différents runs pour maximiser la profondeur
# Stratégie: merge R1+R2 par run, puis concatenation des merged entre runs
# Analyses séparées pour merged et unmerged
################################################################################

set -eo pipefail

################################################################################
# CONFIGURATION GLOBALE
################################################################################

BASE_DIR="/home/plstenge/coprolites_comparison"
BBDUK="/home/plstenge/bbmap/bbduk.sh"
CLUMPIFY="/home/plstenge/bbmap/clumpify.sh"
PHIX="/home/plstenge/bbmap/resources/phix174_ill.ref.fa.gz"
KRAKEN2_DB="/home/plstenge/k2_core_nt_20250609"
KRAKENTOOLS_DIR="${BASE_DIR}/08_kraken2/KrakenTools"
THREADS=36

################################################################################
# ORGANISATION DES DONNÉES - 5 LOTS À MERGER
################################################################################

echo ""
echo "=========================================="
echo "ÉTAPE 0: Organisation des données - 5 LOTS"
echo "=========================================="

mkdir -p "${BASE_DIR}/01_raw_data_merged"
mkdir -p "${BASE_DIR}/00_scripts"

# Définition des 5 lots sources
declare -a SOURCE_LOTS=(
    "Lot1_illu-4_R1_Ps4_150_Default"
    "Lot2_Run1_R2_Ps6_150_no_filter"
    "Lot4_Run2_R2_Ps6_150_no_filter"
    "Lot6_Run3_R3_Ps8_150_no_filter"
    "Lot9_Run4_R3_Ps8_75_no_filter"
)

# Mapping des lots vers leurs chemins sources
declare -A LOT_SOURCE
declare -A LOT_MODE

LOT_SOURCE["Lot1_illu-4_R1_Ps4_150_Default"]="/home/plstenge/coprolites_comparison/01_raw_data/Lot1_Illumina_R1"
LOT_MODE["Lot1_illu-4_R1_Ps4_150_Default"]="flat"

LOT_SOURCE["Lot2_Run1_R2_Ps6_150_no_filter"]="/storage/groups/gdec/shared_paleo/E1531_final/run1_20250320_AV241601_E1531_Ps5Lane1_Ps6Lane2"
LOT_MODE["Lot2_Run1_R2_Ps6_150_no_filter"]="subdirs"

LOT_SOURCE["Lot4_Run2_R2_Ps6_150_no_filter"]="/storage/groups/gdec/shared_paleo/E1531_final/run2_20250414_AV241601_E1531_Ps5_Ps6_14042025"
LOT_MODE["Lot4_Run2_R2_Ps6_150_no_filter"]="subdirs"

LOT_SOURCE["Lot6_Run3_R3_Ps8_150_no_filter"]="/storage/groups/gdec/shared_paleo/E1531_final/run3_20251008_AV241601_E1531_Ps7_Ps8"
LOT_MODE["Lot6_Run3_R3_Ps8_150_no_filter"]="subdirs"

LOT_SOURCE["Lot9_Run4_R3_Ps8_75_no_filter"]="/storage/groups/gdec/shared_paleo/E1531_final/run4_20251104_AV241601_E1531_Ps7_Ps8_04112025"
LOT_MODE["Lot9_Run4_R3_Ps8_75_no_filter"]="subdirs"

# Organisation des données brutes
shopt -s nullglob

for lot in "${SOURCE_LOTS[@]}"; do
    SRC_DIR="${LOT_SOURCE[$lot]}"
    MODE="${LOT_MODE[$lot]}"
    DEST_DIR="${BASE_DIR}/01_raw_data_merged/${lot}"
    
    echo ""
    echo "------------------------------------------"
    echo "Organisation du ${lot}"
    echo "Source : ${SRC_DIR}"
    echo "Destination : ${DEST_DIR}"
    echo "Mode : ${MODE}"
    echo "------------------------------------------"
    
    # Vérification source
    if [[ ! -d "$SRC_DIR" ]]; then
        echo "⚠ ATTENTION: répertoire source introuvable pour ${lot} : ${SRC_DIR}"
        continue
    fi
    
    mkdir -p "$DEST_DIR"
    
    # Mode flat (Lot1 Illumina)
    if [[ "$MODE" == "flat" ]]; then
        echo "Création de liens symboliques pour ${lot} (mode flat)..."
        for fq in "${SRC_DIR}"/cop*_R[12].fastq.gz; do
            base_fq=$(basename "$fq")
            if [[ -e "${DEST_DIR}/${base_fq}" ]]; then
                echo " ↪ Lien déjà présent pour ${base_fq}, on saute."
            else
                ln -s "$fq" "${DEST_DIR}/${base_fq}"
                echo " ✓ Lien créé: ${base_fq}"
            fi
        done
    
    # Mode subdirs (Lots 2, 4, 6, 9)
    elif [[ "$MODE" == "subdirs" ]]; then
        echo "Recherche des sous-dossiers copXXX dans ${SRC_DIR}..."
        for d in "${SRC_DIR}"/[0-9]*_cop[0-9][0-9][0-9]; do
            [[ -d "$d" ]] || continue
            folder_name=$(basename "$d")
            sample_id="${folder_name#*_}"
            
            R1_SRC="${d}/${folder_name}_R1.fastq.gz"
            R2_SRC="${d}/${folder_name}_R2.fastq.gz"
            
            if [[ -f "$R1_SRC" && -f "$R2_SRC" ]]; then
                R1_DEST="${DEST_DIR}/${sample_id}_R1.fastq.gz"
                R2_DEST="${DEST_DIR}/${sample_id}_R2.fastq.gz"
                
                if [[ -e "$R1_DEST" || -e "$R2_DEST" ]]; then
                    echo " ↪ Liens déjà présents pour ${sample_id} dans ${lot}, on saute."
                else
                    ln -s "$R1_SRC" "$R1_DEST"
                    ln -s "$R2_SRC" "$R2_DEST"
                    echo " ✓ ${sample_id} lié (${folder_name}_R1/R2.fastq.gz → ${sample_id}_R1/R2.fastq.gz)"
                fi
            else
                echo " ⚠ Fichiers R1/R2 manquants pour ${folder_name} dans ${SRC_DIR}"
            fi
        done
    fi
done

shopt -u nullglob

echo ""
echo "=========================================="
echo "Organisation des données terminée"
echo "=========================================="

################################################################################
# ACTIVATION ENVIRONNEMENT CONDA
################################################################################

echo ""
echo "=== Activation environnement conda ==="
module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics
echo "Environnement activé: metagenomics"

echo ""
echo "=== Initialisation de la taxonomie Krona ==="
KRONA_TAX_DIR=$(conda env list | grep metagenomics | awk '{print $NF}')/opt/krona/taxonomy
if [[ ! -d "$KRONA_TAX_DIR" ]] || [[ ! -f "$KRONA_TAX_DIR/taxonomy.tab" ]]; then
    echo "Taxonomie Krona absente. Installation en cours..."
    ktUpdateTaxonomy.sh "$KRONA_TAX_DIR"
    echo "Taxonomie Krona installée avec succès."
else
    echo "Taxonomie Krona déjà installée."
fi

################################################################################
# CRÉATION ARBORESCENCE
################################################################################

echo ""
echo "=== Création de l'arborescence ==="

mkdir -p "${BASE_DIR}/02_quality_check_raw/merged"
mkdir -p "${BASE_DIR}/03_bbduk/merged"
mkdir -p "${BASE_DIR}/04_fastuniq/merged"
mkdir -p "${BASE_DIR}/05_clumpify/merged"
mkdir -p "${BASE_DIR}/06_fastp/merged_reads"
mkdir -p "${BASE_DIR}/06_fastp/unmerged_reads"
mkdir -p "${BASE_DIR}/07_quality_check_clean/merged_reads"
mkdir -p "${BASE_DIR}/07_quality_check_clean/unmerged_reads"
mkdir -p "${BASE_DIR}/08_kraken2/merged_reads"
mkdir -p "${BASE_DIR}/08_kraken2/unmerged_reads"
mkdir -p "${BASE_DIR}/09_krona/merged_reads"
mkdir -p "${BASE_DIR}/09_krona/unmerged_reads"
mkdir -p "${BASE_DIR}/10_mpa_tables/merged_reads"
mkdir -p "${BASE_DIR}/10_mpa_tables/unmerged_reads"
mkdir -p "${BASE_DIR}/11_summary_tables"
mkdir -p "${BASE_DIR}/00_intermediate_merged"
mkdir -p "${BASE_DIR}/12_mapdamage/merged_reads"
mkdir -p "${BASE_DIR}/12_mapdamage/unmerged_reads"

echo "Arborescence créée."

#################################################################################
## ÉTAPE 1: NETTOYAGE PAR LOT (BBDuk + FastUniq + Clumpify)
#################################################################################
#
#echo ""
#echo "=========================================="
#echo "ÉTAPE 1: Nettoyage des reads par lot"
#echo "=========================================="
#
#for lot in "${SOURCE_LOTS[@]}"; do
#    echo ""
#    echo "=========================================="
#    echo "Traitement du ${lot}"
#    echo "=========================================="
#    
#    INPUT_DIR="${BASE_DIR}/01_raw_data_merged/${lot}"
#    
#    if [[ ! -d "$INPUT_DIR" ]]; then
#        echo "⚠ Répertoire introuvable: ${INPUT_DIR}, on skip"
#        continue
#    fi
#    
#    # ------------------------------------------------------------------
#    # BBDuk
#    # ------------------------------------------------------------------
#    echo ""
#    echo "→ BBDuk pour ${lot}..."
#    BBDUK_OUT="${BASE_DIR}/03_bbduk/${lot}"
#    mkdir -p "$BBDUK_OUT"
#    
#    cd "$INPUT_DIR"
#    for r1_file in *_R1.fastq.gz; do
#        r2_file="${r1_file/_R1/_R2}"
#        if [[ ! -f "$r2_file" ]]; then
#            echo " ✗ Fichier R2 manquant pour $r1_file"
#            continue
#        fi
#        
#        base_name="${r1_file%%_R1.fastq.gz}"
#        echo " • Traitement de ${base_name}..."
#        
#        $BBDUK -Xmx4g \
#            in1="$r1_file" \
#            in2="$r2_file" \
#            out1="${BBDUK_OUT}/clean_${r1_file}" \
#            out2="${BBDUK_OUT}/clean_${r2_file}" \
#            ref=$PHIX \
#            ktrim=rl k=23 mink=11 hdist=1 \
#            tpe tbo \
#            minlen=25 \
#            qtrim=r trimq=20 \
#            stats="${BBDUK_OUT}/${base_name}_bbduk_stats.txt"
#    done
#    
#    # ------------------------------------------------------------------
#    # FastUniq
#    # ------------------------------------------------------------------
#    echo ""
#    echo "→ FastUniq pour ${lot}..."
#    FASTUNIQ_OUT="${BASE_DIR}/04_fastuniq/${lot}"
#    mkdir -p "$FASTUNIQ_OUT"
#    
#    TMP="/tmp/fastuniq_${lot}_$$"
#    mkdir -p "$TMP"
#    
#    cd "$BBDUK_OUT" || continue
#    for R1_gz in clean_*_R1.fastq.gz; do
#        base=$(echo "$R1_gz" | sed 's/_R1\.fastq\.gz//')
#        R2_gz="${base}_R2.fastq.gz"
#        
#        if [[ -f "$R2_gz" ]]; then
#            echo " • Déduplication de ${base}..."
#            R1_tmp="${TMP}/${base}_R1.fastq"
#            R2_tmp="${TMP}/${base}_R2.fastq"
#            listfile="${TMP}/${base}.list"
#            
#            zcat "$R1_gz" > "$R1_tmp"
#            zcat "$R2_gz" > "$R2_tmp"
#            echo -e "${R1_tmp}\n${R2_tmp}" > "$listfile"
#            
#            fastuniq -i "$listfile" -t q \
#                -o "${FASTUNIQ_OUT}/${base}_dedup_R1.fastq" \
#                -p "${FASTUNIQ_OUT}/${base}_dedup_R2.fastq"
#            
#            rm -f "$R1_tmp" "$R2_tmp" "$listfile"
#        fi
#    done
#    rm -rf "$TMP"
#    
#    # ------------------------------------------------------------------
#    # Clumpify
#    # ------------------------------------------------------------------
#    echo ""
#    echo "→ Clumpify pour ${lot}..."
#    CLUMPIFY_OUT="${BASE_DIR}/05_clumpify/${lot}"
#    mkdir -p "$CLUMPIFY_OUT"
#    
#    for R1 in "${FASTUNIQ_OUT}"/*_R1.fastq; do
#        R2="${R1/_R1.fastq/_R2.fastq}"
#        if [[ -f "$R2" ]]; then
#            base=$(basename "$R1" _R1.fastq)
#            echo " • Clumpify de ${base}..."
#            
#            $CLUMPIFY \
#                in="$R1" in2="$R2" \
#                out="${CLUMPIFY_OUT}/${base}_clumpify_R1.fastq.gz" \
#                out2="${CLUMPIFY_OUT}/${base}_clumpify_R2.fastq.gz" \
#                dedupe=t
#        fi
#    done
#done
#
#################################################################################
## ÉTAPE 2: MERGE PAR RUN (avec fastp) puis CONCATENATION
#################################################################################
#
#echo ""
#echo "=========================================="
#echo "ÉTAPE 2: Merge R1+R2 par run, puis concatenation"
#echo "=========================================="
#
## Créer répertoire intermédiaire pour les merged par run
#INTERMEDIATE_DIR="${BASE_DIR}/00_intermediate_merged"
#mkdir -p "$INTERMEDIATE_DIR"
#
## Liste des échantillons uniques
#declare -a ALL_SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")
#
## Étape 2a: Merge R1+R2 pour chaque lot/run
#echo ""
#echo "→ Étape 2a: Merge R1+R2 pour chaque run..."
#
#for lot in "${SOURCE_LOTS[@]}"; do
#    echo ""
#    echo "Fastp merge pour ${lot}..."
#    INPUT_DIR="${BASE_DIR}/05_clumpify/${lot}"
#    OUTPUT_DIR="${INTERMEDIATE_DIR}/${lot}"
#    mkdir -p "$OUTPUT_DIR"
#    
#    if [[ ! -d "$INPUT_DIR" ]]; then
#        echo "⚠ Pas de données clumpify pour ${lot}"
#        continue
#    fi
#    
#    for R1 in "${INPUT_DIR}"/*_R1.fastq.gz; do
#        R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
#        if [[ -f "$R2" ]]; then
#            base=$(basename "$R1" _R1.fastq.gz)
#            echo " • Merge de ${base}..."
#            
#            fastp \
#                -i "$R1" -I "$R2" \
#                --merged_out "${OUTPUT_DIR}/${base}_merged.fastq.gz" \
#                --out1 "${OUTPUT_DIR}/${base}_unmerged_R1.fastq.gz" \
#                --out2 "${OUTPUT_DIR}/${base}_unmerged_R2.fastq.gz" \
#                --json "${OUTPUT_DIR}/${base}_fastp.json" \
#                --html "${OUTPUT_DIR}/${base}_fastp.html" \
#                --thread 4 \
#                --length_required 30 \
#                --qualified_quality_phred 20
#        fi
#    done
#done
#
## Étape 2b: Concatenation des fichiers merged entre runs
#echo ""
#echo "→ Étape 2b: Concatenation des fichiers merged et unmerged entre runs..."
#
#FINAL_MERGED_DIR="${BASE_DIR}/06_fastp/merged_reads"
#FINAL_UNMERGED_DIR="${BASE_DIR}/06_fastp/unmerged_reads"
#mkdir -p "$FINAL_MERGED_DIR"
#mkdir -p "$FINAL_UNMERGED_DIR"
#
#for sample in "${ALL_SAMPLES[@]}"; do
#    echo ""
#    echo "=========================================="
#    echo "Concatenation pour ${sample}"
#    echo "=========================================="
#    
#    # Rechercher tous les fichiers merged pour cet échantillon
#    merged_files=()
#    unmerged_r1_files=()
#    unmerged_r2_files=()
#    
#    for lot in "${SOURCE_LOTS[@]}"; do
#        lot_dir="${INTERMEDIATE_DIR}/${lot}"
#        
#        # Chercher les fichiers merged (sans le préfixe "clean_")
#        for mf in "${lot_dir}"/clean_${sample}_dedup_clumpify_merged.fastq.gz; do
#            if [[ -f "$mf" ]]; then
#                merged_files+=("$mf")
#                echo " ✓ Trouvé merged: $(basename $mf) dans ${lot}"
#            fi
#        done
#        
#        # Chercher les fichiers unmerged
#        for uf1 in "${lot_dir}"/clean_${sample}_dedup_clumpify_unmerged_R1.fastq.gz; do
#            if [[ -f "$uf1" ]]; then
#                uf2="${uf1/_R1.fastq.gz/_R2.fastq.gz}"
#                if [[ -f "$uf2" ]]; then
#                    unmerged_r1_files+=("$uf1")
#                    unmerged_r2_files+=("$uf2")
#                    echo " ✓ Trouvé unmerged: $(basename $uf1) dans ${lot}"
#                fi
#            fi
#        done
#    done
#    
#    # Concatenation des merged
#    if [[ ${#merged_files[@]} -gt 0 ]]; then
#        echo " → Concatenation de ${#merged_files[@]} fichiers merged..."
#        cat "${merged_files[@]}" > "${FINAL_MERGED_DIR}/${sample}_final_merged.fastq.gz"
#        echo " ✓ Créé: ${sample}_final_merged.fastq.gz"
#    else
#        echo " ⚠ Aucun fichier merged trouvé pour ${sample}"
#    fi
#    
#    # Concatenation des unmerged
#    if [[ ${#unmerged_r1_files[@]} -gt 0 ]]; then
#        echo " → Concatenation de ${#unmerged_r1_files[@]} paires unmerged..."
#        cat "${unmerged_r1_files[@]}" > "${FINAL_UNMERGED_DIR}/${sample}_final_unmerged_R1.fastq.gz"
#        cat "${unmerged_r2_files[@]}" > "${FINAL_UNMERGED_DIR}/${sample}_final_unmerged_R2.fastq.gz"
#        echo " ✓ Créé: ${sample}_final_unmerged_R1/R2.fastq.gz"
#    else
#        echo " ⚠ Aucun fichier unmerged trouvé pour ${sample}"
#    fi
#done
#
#################################################################################
## ÉTAPE 3: Contrôle qualité MERGED et UNMERGED
#################################################################################
#
#echo ""
#echo "=========================================="
#echo "ÉTAPE 3: Contrôle qualité des données merged et unmerged"
#echo "=========================================="
#
## QC pour les merged
#QC_MERGED_DIR="${BASE_DIR}/07_quality_check_clean/merged_reads"
#mkdir -p "$QC_MERGED_DIR"
#
#echo "FastQC sur les fichiers merged..."
#for FILE in "${FINAL_MERGED_DIR}"/*.fastq.gz; do
#    if [[ -f "$FILE" ]]; then
#        echo " • FastQC sur $(basename $FILE)..."
#        fastqc "$FILE" -o "$QC_MERGED_DIR" -t 4
#    fi
#done
#
#echo "MultiQC merged..."
#cd "$QC_MERGED_DIR"
#multiqc . -n "multiqc_merged_reads.html" --force
#
## QC pour les unmerged
#QC_UNMERGED_DIR="${BASE_DIR}/07_quality_check_clean/unmerged_reads"
#mkdir -p "$QC_UNMERGED_DIR"
#
#echo "FastQC sur les fichiers unmerged..."
#for FILE in "${FINAL_UNMERGED_DIR}"/*.fastq.gz; do
#    if [[ -f "$FILE" ]]; then
#        echo " • FastQC sur $(basename $FILE)..."
#        fastqc "$FILE" -o "$QC_UNMERGED_DIR" -t 4
#    fi
#done
#
#echo "MultiQC unmerged..."
#cd "$QC_UNMERGED_DIR"
#multiqc . -n "multiqc_unmerged_reads.html" --force
#
#################################################################################
## ÉTAPE 4: Classification Kraken2 - MERGED
#################################################################################
#
#echo ""
#echo "=========================================="
#echo "ÉTAPE 4a: Classification taxonomique (Kraken2) - MERGED READS"
#echo "=========================================="
#
#KRAKEN_MERGED="${BASE_DIR}/08_kraken2/merged_reads"
#mkdir -p "$KRAKEN_MERGED"
#
#for sample in "${ALL_SAMPLES[@]}"; do
#    echo ""
#    echo "Kraken2 pour ${sample} (merged)..."
#    
#    # Analyse merged
#    MERGED="${FINAL_MERGED_DIR}/${sample}_final_merged.fastq.gz"
#    if [[ -f "$MERGED" ]]; then
#        echo " → Analyse merged..."
#        kraken2 --confidence 0.2 --db "$KRAKEN2_DB" --threads $THREADS \
#            --output "${KRAKEN_MERGED}/${sample}_merged.kraken" \
#            --report "${KRAKEN_MERGED}/${sample}_merged.report" \
#            "$MERGED"
#    else
#        echo " ⚠ Fichier merged absent pour ${sample}"
#    fi
#done
#
#################################################################################
## ÉTAPE 4b: Classification Kraken2 - UNMERGED
#################################################################################
#
#echo ""
#echo "=========================================="
#echo "ÉTAPE 4b: Classification taxonomique (Kraken2) - UNMERGED READS"
#echo "=========================================="
#
#KRAKEN_UNMERGED="${BASE_DIR}/08_kraken2/unmerged_reads"
#mkdir -p "$KRAKEN_UNMERGED"
#
#for sample in "${ALL_SAMPLES[@]}"; do
#    echo ""
#    echo "Kraken2 pour ${sample} (unmerged)..."
#    
#    # Analyse unmerged
#    R1="${FINAL_UNMERGED_DIR}/${sample}_final_unmerged_R1.fastq.gz"
#    R2="${FINAL_UNMERGED_DIR}/${sample}_final_unmerged_R2.fastq.gz"
#    if [[ -f "$R1" && -f "$R2" ]]; then
#        echo " → Analyse unmerged..."
#        kraken2 --confidence 0.2 --paired --db "$KRAKEN2_DB" --threads $THREADS \
#            --output "${KRAKEN_UNMERGED}/${sample}_unmerged.kraken" \
#            --report "${KRAKEN_UNMERGED}/${sample}_unmerged.report" \
#            "$R1" "$R2"
#    else
#        echo " ⚠ Fichiers unmerged absents pour ${sample}"
#    fi
#done
#
#################################################################################
## ÉTAPE 5: Visualisation Krona - MERGED
#################################################################################
#
#echo ""
#echo "=========================================="
#echo "ÉTAPE 5a: Visualisation (Krona) - MERGED READS"
#echo "=========================================="
#
#KRONA_MERGED="${BASE_DIR}/09_krona/merged_reads"
#mkdir -p "$KRONA_MERGED"
#
#cd "$KRAKEN_MERGED"
#
## Krona combiné pour tous les échantillons merged
#if ls *.report 1> /dev/null 2>&1; then
#    echo " → Génération Krona combiné merged..."
#    ktImportTaxonomy -t 5 -m 3 -o "${KRONA_MERGED}/all_samples_merged_krona.html" *.report
#fi
#
## Krona individuel pour chaque échantillon merged
#for report in *.report; do
#    if [[ -f "$report" ]]; then
#        base=$(basename "$report" .report)
#        echo " → Génération Krona pour ${base}..."
#        ktImportTaxonomy -t 5 -m 3 -o "${KRONA_MERGED}/${base}_krona.html" "$report"
#    fi
#done
#
#################################################################################
## ÉTAPE 5b: Visualisation Krona - UNMERGED
#################################################################################
#
#echo ""
#echo "=========================================="
#echo "ÉTAPE 5b: Visualisation (Krona) - UNMERGED READS"
#echo "=========================================="
#
#KRONA_UNMERGED="${BASE_DIR}/09_krona/unmerged_reads"
#mkdir -p "$KRONA_UNMERGED"
#
#cd "$KRAKEN_UNMERGED"
#
## Krona combiné pour tous les échantillons unmerged
#if ls *.report 1> /dev/null 2>&1; then
#    echo " → Génération Krona combiné unmerged..."
#    ktImportTaxonomy -t 5 -m 3 -o "${KRONA_UNMERGED}/all_samples_unmerged_krona.html" *.report
#fi
#
## Krona individuel pour chaque échantillon unmerged
#for report in *.report; do
#    if [[ -f "$report" ]]; then
#        base=$(basename "$report" .report)
#        echo " → Génération Krona pour ${base}..."
#        ktImportTaxonomy -t 5 -m 3 -o "${KRONA_UNMERGED}/${base}_krona.html" "$report"
#    fi
#done
#
#################################################################################
## ÉTAPE 6: Tables MPA - MERGED
#################################################################################
#
#echo ""
#echo "=========================================="
#echo "ÉTAPE 6a: Création des tables MPA - MERGED READS"
#echo "=========================================="
#
#if [[ ! -d "$KRAKENTOOLS_DIR" ]]; then
#    echo "Installation de KrakenTools..."
#    mkdir -p "${BASE_DIR}/08_kraken2"
#    cd "${BASE_DIR}/08_kraken2"
#    git clone https://github.com/jenniferlu717/KrakenTools.git
#fi
#
#MPA_MERGED="${BASE_DIR}/10_mpa_tables/merged_reads"
#mkdir -p "$MPA_MERGED"
#
#cd "$KRAKEN_MERGED"
#
#declare -a mpa_files_merged=()
#for report in *.report; do
#    if [[ -f "$report" ]]; then
#        base=$(basename "$report" .report)
#        mpa_file="${MPA_MERGED}/${base}.mpa"
#        echo " → Conversion de ${base}..."
#        python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" -r "$report" -o "$mpa_file"
#        mpa_files_merged+=("$mpa_file")
#    fi
#done
#
#if [[ ${#mpa_files_merged[@]} -gt 0 ]]; then
#    echo " → Combinaison de tous les fichiers MPA merged..."
#    python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" -i "${mpa_files_merged[@]}" -o "${MPA_MERGED}/combined_all_merged.tsv"
#fi
#
#################################################################################
## ÉTAPE 6b: Tables MPA - UNMERGED
#################################################################################
#
#echo ""
#echo "=========================================="
#echo "ÉTAPE 6b: Création des tables MPA - UNMERGED READS"
#echo "=========================================="
#
#MPA_UNMERGED="${BASE_DIR}/10_mpa_tables/unmerged_reads"
#mkdir -p "$MPA_UNMERGED"
#
#cd "$KRAKEN_UNMERGED"
#
#declare -a mpa_files_unmerged=()
#for report in *.report; do
#    if [[ -f "$report" ]]; then
#        base=$(basename "$report" .report)
#        mpa_file="${MPA_UNMERGED}/${base}.mpa"
#        echo " → Conversion de ${base}..."
#        python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" -r "$report" -o "$mpa_file"
#        mpa_files_unmerged+=("$mpa_file")
#    fi
#done
#
#if [[ ${#mpa_files_unmerged[@]} -gt 0 ]]; then
#    echo " → Combinaison de tous les fichiers MPA unmerged..."
#    python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" -i "${mpa_files_unmerged[@]}" -o "${MPA_UNMERGED}/combined_all_unmerged.tsv"
#fi
#
#################################################################################
## ÉTAPE 7: Tableau récapitulatif
#################################################################################
#
#echo ""
#echo "=========================================="
#echo "ÉTAPE 7: Tableau récapitulatif"
#echo "=========================================="
#
#SUMMARY_TABLE="${BASE_DIR}/11_summary_tables/sequences_summary_merged.tsv"
#
#cat > "$SUMMARY_TABLE" << 'HEADER'
#Sample	Type	Nb_sequences	Longueur_moyenne	GC_percent
#HEADER
#
#function extract_stats() {
#    local file=$1
#    local sample=$2
#    local type=$3
#    
#    if [[ ! -f "$file" ]]; then
#        return 1
#    fi
#    
#    nb_seq=$(zcat "$file" 2>/dev/null | echo $((`wc -l`/4)))
#    stats=$(zcat "$file" 2>/dev/null | awk 'NR%4==2 {
#        total_len += length($0)
#        gc_count += gsub(/[GCgc]/, "", $0)
#        at_count += gsub(/[ATat]/, "", $0)
#        count++
#    } END {
#        if (count > 0) {
#            avg_len = total_len / count
#            gc_perc = (gc_count / (gc_count + at_count)) * 100
#            printf "%.1f\t%.2f", avg_len, gc_perc
#        } else {
#            printf "0\t0"
#        }
#    }')
#    
#    echo -e "${sample}\t${type}\t${nb_seq}\t${stats}" >> "$SUMMARY_TABLE"
#}
#
#echo "Calcul des statistiques..."
#for sample in "${ALL_SAMPLES[@]}"; do
#    echo " → ${sample}..."
#    
#    merged="${FINAL_MERGED_DIR}/${sample}_final_merged.fastq.gz"
#    if [[ -f "$merged" ]]; then
#        extract_stats "$merged" "$sample" "merged"
#    fi
#    
#    unmerged_r1="${FINAL_UNMERGED_DIR}/${sample}_final_unmerged_R1.fastq.gz"
#    if [[ -f "$unmerged_r1" ]]; then
#        extract_stats "$unmerged_r1" "$sample" "unmerged_R1"
#    fi
#    
#    unmerged_r2="${FINAL_UNMERGED_DIR}/${sample}_final_unmerged_R2.fastq.gz"
#    if [[ -f "$unmerged_r2" ]]; then
#        extract_stats "$unmerged_r2" "$sample" "unmerged_R2"
#    fi
#done
#
#echo ""
#echo "Tableau récapitulatif créé: ${SUMMARY_TABLE}"
#echo ""
#echo "Aperçu:"
#head -20 "$SUMMARY_TABLE" | column -t

################################################################################
# ÉTAPE 8: MapDamage (Ovis aries & Capra hircus & Phragmites)
################################################################################

echo ""
echo "=========================================="
echo "ÉTAPE 8: MapDamage - Analyse des dommages"
echo "=========================================="

# Activer l'environnement mapdamage
conda deactivate
conda activate mapdamage_py39

DAMAGE_MERGED_BASE="${BASE_DIR}/12_mapdamage/merged_reads"
DAMAGE_UNMERGED_BASE="${BASE_DIR}/12_mapdamage/unmerged_reads"
mkdir -p "$DAMAGE_MERGED_BASE"
mkdir -p "$DAMAGE_UNMERGED_BASE"

LOGFILE="${BASE_DIR}/00_scripts/mapdamage_merged_$(date +%Y%m%d_%H%M%S).txt"
touch "$LOGFILE"
MAPPING_INFO="${BASE_DIR}/11_summary_tables/mapping_bwa_info_merged.tsv"

echo "Script MapDamage started at $(date)" | tee -a "$LOGFILE"

# Initialiser le fichier de mapping info
echo -e "Sample\tSpecies\tType\tTotal_Reads\tMapped_Reads\tMapping_Rate(%)" > "$MAPPING_INFO"

# Définition des génomes de référence
declare -A TAXONS=(
    ["Ovis_aries"]="9940:/home/plstenge/genomes/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa"
    ["Capra_hircus"]="9925:/home/plstenge/genomes/Capra_hircus.ARS1.dna.toplevel.fa"
    ["Bos_taurus"]="9903:/home/plstenge/genomes/Bos_taurus/GCF_002263795.3_ARS-UCD2.0_genomic.fna"
     ["Alnus_glutinosa"]="3517:/home/plstenge/genomes/Alnus_glutinosa_genome_assembly_dhAlnGlut1.fa"
     ["Phragmites_australis"]="29695:/home/plstenge/genomes/Phragmites_australis_GCA_040373225.1.fasta"
     ["Corylus_avellana"]="13451:/home/plstenge/genomes/Corylus_avellana_CavTom2PMs_1_0.fasta"
)


################################################################################
# FONCTION: Nettoyer les noms de séquences dupliqués dans un FASTA
################################################################################
fix_duplicate_fasta_headers() {
    local input_fasta="$1"
    local output_fasta="$2"
    
    echo " → Vérification des headers dupliqués dans le génome..."
    
    # Vérifier s'il y a des duplications
    duplicates=$(grep "^>" "$input_fasta" | sort | uniq -d | wc -l)
    
    if [[ $duplicates -gt 0 ]]; then
        echo " ⚠ ${duplicates} headers dupliqués détectés. Correction en cours..."
        
        # Créer un FASTA avec des noms uniques
        awk '/^>/ {
            header=$0
            count[header]++
            if (count[header] > 1) {
                print header "_dup" count[header]
            } else {
                print header
            }
            next
        }
        {print}' "$input_fasta" > "$output_fasta"
        
        echo " ✓ Headers corrigés dans: $output_fasta"
        return 0
    else
        echo " ✓ Aucun header dupliqué détecté"
        # Si pas de duplication, créer un lien symbolique
        ln -sf "$input_fasta" "$output_fasta" 2>/dev/null || cp "$input_fasta" "$output_fasta"
        return 1
    fi
}

# Fonction pour calculer le taux de mapping
calculate_mapping_rate() {
    local bam_file="$1"
    local sample_name="$2"
    local species="$3"
    local type="$4"
    
    if [[ -f "$bam_file" ]]; then
        local total_reads=$(samtools view -c "$bam_file" 2>/dev/null || echo 0)
        local mapped_reads=$(samtools view -c -F 4 "$bam_file" 2>/dev/null || echo 0)
        local mapping_rate=0
        
        if [[ $total_reads -gt 0 ]]; then
            mapping_rate=$(echo "scale=2; $mapped_reads * 100 / $total_reads" | bc)
        fi
        
        echo -e "${sample_name}\t${species}\t${type}\t${total_reads}\t${mapped_reads}\t${mapping_rate}" >> "$MAPPING_INFO"
        echo "✓ Mapping stats for ${sample_name}_${species}_${type}: ${mapped_reads}/${total_reads} (${mapping_rate}%)" | tee -a "$LOGFILE"
    fi
}

# Fonction MapDamage sécurisée
run_mapdamage_safe() {
    local bam_file="$1"
    local ref_fasta="$2"
    local output_dir="$3"
    local sample_name="$4"
    
    echo " → MapDamage pour ${sample_name}..."
    
    if [[ ! -s "$bam_file" ]]; then
        echo " ✗ BAM vide ou absent, on skip MapDamage"
        return 0
    fi
    
    mkdir -p "$output_dir"
    
    mapDamage -i "$bam_file" -r "$ref_fasta" --folder "$output_dir" --no-stats 2>>"$LOGFILE" || {
        echo " ⚠ MapDamage a échoué pour ${sample_name} (erreur ignorée)"
        return 0
    }
    
    echo " ✓ MapDamage terminé pour ${sample_name}"
    return 0
}

################################################################################
# PRÉ-TRAITEMENT: Vérifier et corriger les génomes de référence
################################################################################

echo ""
echo "→ Vérification des génomes de référence..."

declare -A FIXED_GENOMES

for GROUP in "${!TAXONS[@]}"; do
    ORIGINAL_REF="${TAXONS[$GROUP]#*:}"
    FIXED_REF="/home/plstenge/genomes/$(basename ${ORIGINAL_REF%.fasta}).fixed.fasta"
    
    echo ""
    echo "-- Traitement du génome: ${GROUP} --"
    
    # Vérifier que le génome original existe
    if [[ ! -f "$ORIGINAL_REF" ]]; then
        echo " ✗ Génome original non trouvé: $ORIGINAL_REF"
        continue
    fi
    
    # Corriger les duplications si nécessaire
    fix_duplicate_fasta_headers "$ORIGINAL_REF" "$FIXED_REF"
    
    # Stocker le chemin du génome corrigé
    FIXED_GENOMES[$GROUP]="$FIXED_REF"
    
    # Vérifier/créer les index BWA
    if [[ ! -f "${FIXED_REF}.bwt" ]]; then
        echo " → Indexation BWA du génome corrigé..."
        bwa index "$FIXED_REF" 2>>"$LOGFILE" || {
            echo " ✗ Indexation BWA échouée pour ${GROUP}"
            continue
        }
        echo " ✓ Index BWA créé"
    else
        echo " ✓ Index BWA déjà présent"
    fi
    
    # Vérifier/créer l'index samtools
    if [[ ! -f "${FIXED_REF}.fai" ]]; then
        echo " → Indexation samtools du génome corrigé..."
        samtools faidx "$FIXED_REF" 2>>"$LOGFILE" || {
            echo " ✗ Indexation samtools échouée pour ${GROUP}"
            continue
        }
        echo " ✓ Index samtools créé"
    else
        echo " ✓ Index samtools déjà présent"
    fi
done

shopt -s nullglob

################################################################################
# ÉTAPE 8a: MapDamage - MERGED READS
################################################################################

echo ""
echo "→ MapDamage pour MERGED READS..."

KRAKEN_MERGED="${BASE_DIR}/08_kraken2/merged_reads"
FINAL_MERGED_DIR="${BASE_DIR}/06_fastp/merged_reads"

for KRAKEN_FILE in "${KRAKEN_MERGED}"/*.kraken; do
    [[ -f "$KRAKEN_FILE" ]] || continue
    KRAKEN_BASE_NAME=$(basename "$KRAKEN_FILE" .kraken)
    
    echo ""
    echo "=========================================="
    echo "Traitement: ${KRAKEN_BASE_NAME} (merged)"
    echo "=========================================="
    
    SAMPLE="${KRAKEN_BASE_NAME%_merged}"
    FASTQ_FILE="${FINAL_MERGED_DIR}/${SAMPLE}_final_merged.fastq.gz"
    
    if [[ ! -f "$FASTQ_FILE" || ! -s "$FASTQ_FILE" ]]; then
        echo "⚠ Fichier FASTQ merged absent ou vide: ${FASTQ_FILE}"
        continue
    fi
    
    # Boucle sur les espèces
    for GROUP in "${!TAXONS[@]}"; do
        REF_FASTA="${FIXED_GENOMES[$GROUP]}"
        [[ -z "$REF_FASTA" ]] && continue
        
        TAX_ID="${TAXONS[$GROUP]%:*}"
        
        DAMAGE_DIR="${DAMAGE_MERGED_BASE}/${GROUP}"
        mkdir -p "$DAMAGE_DIR"
        
        echo ""
        echo "-- Espèce: ${GROUP} (TaxID: ${TAX_ID}) --"
        
        OUT_FASTQ="${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.fastq"
        
        echo "→ Extraction reads merged..."
        python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
            -k "$KRAKEN_FILE" -s "$FASTQ_FILE" -t "$TAX_ID" \
            -o "$OUT_FASTQ" --fastq-output 2>>"$LOGFILE"
        
        if [[ -f "$OUT_FASTQ" && -s "$OUT_FASTQ" ]]; then
            echo "→ Mapping BWA single-end..."
            bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF_FASTA" "$OUT_FASTQ" > "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.sai" 2>>"$LOGFILE"
            
            bwa samse "$REF_FASTA" "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.sai" "$OUT_FASTQ" 2>>"$LOGFILE" | \
            samtools view -bS - 2>>"$LOGFILE" | \
            samtools sort -o "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.sorted.bam" - 2>>"$LOGFILE"
            
            samtools index "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.sorted.bam" 2>>"$LOGFILE"
            
            calculate_mapping_rate "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.sorted.bam" "$SAMPLE" "$GROUP" "merged"
            
            rm -f "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.sai"
            
            run_mapdamage_safe "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.sorted.bam" "$REF_FASTA" \
                "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_mapDamage" "$KRAKEN_BASE_NAME"
        else
            echo " ⚠ Aucun read extrait pour ${GROUP}"
        fi
    done
done

################################################################################
# ÉTAPE 8b: MapDamage - UNMERGED READS
################################################################################

echo ""
echo "→ MapDamage pour UNMERGED READS..."

KRAKEN_UNMERGED="${BASE_DIR}/08_kraken2/unmerged_reads"
FINAL_UNMERGED_DIR="${BASE_DIR}/06_fastp/unmerged_reads"

for KRAKEN_FILE in "${KRAKEN_UNMERGED}"/*.kraken; do
    [[ -f "$KRAKEN_FILE" ]] || continue
    KRAKEN_BASE_NAME=$(basename "$KRAKEN_FILE" .kraken)
    
    echo ""
    echo "=========================================="
    echo "Traitement: ${KRAKEN_BASE_NAME} (unmerged)"
    echo "=========================================="
    
    SAMPLE="${KRAKEN_BASE_NAME%_unmerged}"
    R1_FILE="${FINAL_UNMERGED_DIR}/${SAMPLE}_final_unmerged_R1.fastq.gz"
    R2_FILE="${FINAL_UNMERGED_DIR}/${SAMPLE}_final_unmerged_R2.fastq.gz"
    
    if [[ ! -f "$R1_FILE" || ! -f "$R2_FILE" || ! -s "$R1_FILE" || ! -s "$R2_FILE" ]]; then
        echo "⚠ Fichiers FASTQ unmerged absents ou vides"
        continue
    fi
    
    # Boucle sur les espèces
    for GROUP in "${!TAXONS[@]}"; do
        REF_FASTA="${FIXED_GENOMES[$GROUP]}"
        [[ -z "$REF_FASTA" ]] && continue
        
        TAX_ID="${TAXONS[$GROUP]%:*}"
        
        DAMAGE_DIR="${DAMAGE_UNMERGED_BASE}/${GROUP}"
        mkdir -p "$DAMAGE_DIR"
        
        echo ""
        echo "-- Espèce: ${GROUP} (TaxID: ${TAX_ID}) --"
        
        OUT_R1="${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_R1.fastq"
        OUT_R2="${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_R2.fastq"
        
        echo "→ Extraction reads paired-end..."
        python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
            -k "$KRAKEN_FILE" -s "$R1_FILE" -s2 "$R2_FILE" -t "$TAX_ID" \
            -o "$OUT_R1" -o2 "$OUT_R2" --fastq-output 2>>"$LOGFILE"
        
        if [[ -f "$OUT_R1" && -f "$OUT_R2" && -s "$OUT_R1" && -s "$OUT_R2" ]]; then
            echo "→ Mapping BWA paired-end..."
            bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF_FASTA" "$OUT_R1" > "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_R1.sai" 2>>"$LOGFILE"
            bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF_FASTA" "$OUT_R2" > "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_R2.sai" 2>>"$LOGFILE"
            
            bwa sampe "$REF_FASTA" \
                "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_R1.sai" \
                "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_R2.sai" \
                "$OUT_R1" "$OUT_R2" 2>>"$LOGFILE" | \
            samtools view -bS - 2>>"$LOGFILE" | \
            samtools sort -o "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.sorted.bam" - 2>>"$LOGFILE"
            
            samtools index "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.sorted.bam" 2>>"$LOGFILE"
            
            calculate_mapping_rate "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.sorted.bam" "$SAMPLE" "$GROUP" "unmerged"
            
            rm -f "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_R1.sai" \
                  "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_R2.sai"
            
            run_mapdamage_safe "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.sorted.bam" "$REF_FASTA" \
                "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_mapDamage" "$KRAKEN_BASE_NAME"
        else
            echo " ⚠ Aucun read extrait pour ${GROUP}"
        fi
    done
done

shopt -u nullglob

echo ""
echo "=========================================="
echo "MapDamage terminé!"
echo "=========================================="
echo ""
echo "Résultats merged: ${DAMAGE_MERGED_BASE}"
echo "Résultats unmerged: ${DAMAGE_UNMERGED_BASE}"
echo "Statistiques: ${MAPPING_INFO}"

################################################################################
# FIN
################################################################################

echo ""
echo "=========================================="
echo "PIPELINE TERMINÉ AVEC SUCCÈS"
echo "Date de fin: $(date)"
echo "=========================================="
echo ""
echo "Résultats: ${BASE_DIR}"
echo "Tableau récapitulatif: ${SUMMARY_TABLE}"
echo "Statistiques mapping: ${MAPPING_INFO}"
echo ""

conda deactivate

echo "Pipeline terminé le $(date). Résultats: ${BASE_DIR}" | \
    mail -s "Pipeline Coprolites MERGED - Terminé" pierrelouis.stenger@gmail.com 2>/dev/null || true

exit 0
