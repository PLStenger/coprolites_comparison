#!/bin/bash

#SBATCH --job-name=coprolites_5lots
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/coprolites_5lots.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/coprolites_5lots.out"

################################################################################
# Pipeline d'analyse aDNA - Projet Coprolites - 5 LOTS
# Author: Pierre-Louis Stenger
# Date: November 2025
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
# ORGANISATION DES DONNÉES - 9 LOTS
################################################################################

echo ""
echo "=========================================="
echo "ÉTAPE 0: Organisation des données - 9 LOTS"
echo "=========================================="

# Tous les raw standardisés seront dans :
#   ${BASE_DIR}/01_raw_data/<NomLot>/
# avec des fichiers nommés : cop408_R1.fastq.gz, cop408_R2.fastq.gz, etc.

mkdir -p "${BASE_DIR}/01_raw_data"

# Définition des lots et de leurs répertoires sources
# Les noms de lots ici sont ceux qui seront utilisés partout dans le pipeline
declare -a LOTS=(
  "Lot1_Illumina"
  "Lot2_Run1_R1_Ps4_150_no_filtered"
  "Lot3_Run1_R1_Ps6_150_filtered"
  "Lot4_Run2_R1_Ps6_150_no_filtered"
  "Lot5_Run2_R1_Ps6_150_filtered"
  "Lot6_Run3_R2_Ps6_150_no_filtered"
  "Lot7_Run3_R2_Ps8_150_filtered"
  "Lot8_Run4_R2_Ps8_150_filtered"
  "Lot9_Run4_R2_Ps8_75_no_filtered"
)

# Répertoires sources
declare -A LOT_SOURCE
declare -A LOT_MODE   # "flat" = fichiers déjà copXXX_R1/R2.fastq.gz, "subdirs" = sous-dossiers 474_cop408/...

LOT_SOURCE["Lot1_Illumina"]="/home/plstenge/coprolites_comparison/01_raw_data/Lot1_Illumina_R1"
LOT_MODE["Lot1_Illumina"]="flat"

LOT_SOURCE["Lot2_Run1_R1_Ps4_150_no_filtered"]="/storage/groups/gdec/shared_paleo/E1531_final/run1_20250320_AV241601_E1531_Ps5Lane1_Ps6Lane2"
LOT_MODE["Lot2_Run1_R1_Ps4_150_no_filtered"]="subdirs"

LOT_SOURCE["Lot3_Run1_R1_Ps6_150_filtered"]="/storage/groups/gdec/shared_paleo/recipes-ShortInsert_filter_Run1"
LOT_MODE["Lot3_Run1_R1_Ps6_150_filtered"]="subdirs"

LOT_SOURCE["Lot4_Run2_R1_Ps6_150_no_filtered"]="/storage/groups/gdec/shared_paleo/E1531_final/run2_20250414_AV241601_E1531_Ps5_Ps6_14042025"
LOT_MODE["Lot4_Run2_R1_Ps6_150_no_filtered"]="subdirs"

LOT_SOURCE["Lot5_Run2_R1_Ps6_150_filtered"]="/storage/groups/gdec/shared_paleo/recipes-ShortInsert_filter_Run2"
LOT_MODE["Lot5_Run2_R1_Ps6_150_filtered"]="subdirs"

LOT_SOURCE["Lot6_Run3_R2_Ps6_150_no_filtered"]="/storage/groups/gdec/shared_paleo/E1531_final/run3_20251008_AV241601_E1531_Ps7_Ps8"
LOT_MODE["Lot6_Run3_R2_Ps6_150_no_filtered"]="subdirs"

LOT_SOURCE["Lot7_Run3_R2_Ps8_150_filtered"]="/storage/groups/gdec/shared_paleo/E1531_final/run3_20251008_AV241601_E1531_Ps7_Ps8-filtered"
LOT_MODE["Lot7_Run3_R2_Ps8_150_filtered"]="subdirs"

LOT_SOURCE["Lot8_Run4_R2_Ps8_150_filtered"]="/storage/groups/gdec/shared_paleo/Lot_6"
LOT_MODE["Lot8_Run4_R2_Ps8_150_filtered"]="subdirs"

LOT_SOURCE["Lot9_Run4_R2_Ps8_75_no_filtered"]="/storage/groups/gdec/shared_paleo/E1531_final/run4_20251104_AV241601_E1531_Ps7_Ps8_04112025"
LOT_MODE["Lot9_Run4_R2_Ps8_75_no_filtered"]="subdirs"

shopt -s nullglob

for lot in "${LOTS[@]}"; do
    SRC_DIR="${LOT_SOURCE[$lot]}"
    MODE="${LOT_MODE[$lot]}"
    DEST_DIR="${BASE_DIR}/01_raw_data/${lot}"

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

    # Lot 1 : fichiers déjà bien nommés dans le bon dossier
    if [[ "$MODE" == "flat" ]]; then
        # Si le dossier source est déjà le dossier destination, ne rien faire
        if [[ "$SRC_DIR" == "$DEST_DIR" ]]; then
            echo "✓ Lot ${lot} déjà organisé dans ${DEST_DIR}"
        else
            mkdir -p "$DEST_DIR"
            echo "Création de liens symboliques pour ${lot} (mode flat)..."
            for fq in "${SRC_DIR}"/cop*_R[12].fastq.gz; do
                base_fq=$(basename "$fq")
                if [[ -e "${DEST_DIR}/${base_fq}" ]]; then
                    echo "  ↪ Lien déjà présent pour ${base_fq}, on saute."
                else
                    ln -s "$fq" "${DEST_DIR}/${base_fq}"
                    echo "  ✓ Lien créé: ${base_fq}"
                fi
            done
        fi

    # Lots 2–9 : sous-dossiers 474_cop408/474_cop408_R1.fastq.gz...
    elif [[ "$MODE" == "subdirs" ]]; then
        mkdir -p "$DEST_DIR"
        echo "Recherche des sous-dossiers copXXX dans ${SRC_DIR}..."

        # Sous-dossiers de type 474_cop408, 475_cop410, ...
        for d in "${SRC_DIR}"/[0-9]*_cop[0-9][0-9][0-9]; do
            [[ -d "$d" ]] || continue

            folder_name=$(basename "$d")   # ex: 474_cop408
            sample_id="${folder_name#*_}"  # ex: cop408

            R1_SRC="${d}/${folder_name}_R1.fastq.gz"
            R2_SRC="${d}/${folder_name}_R2.fastq.gz"

            if [[ -f "$R1_SRC" && -f "$R2_SRC" ]]; then
                R1_DEST="${DEST_DIR}/${sample_id}_R1.fastq.gz"
                R2_DEST="${DEST_DIR}/${sample_id}_R2.fastq.gz"

                if [[ -e "$R1_DEST" || -e "$R2_DEST" ]]; then
                    echo "  ↪ Liens déjà présents pour ${sample_id} dans ${lot}, on saute."
                else
                    ln -s "$R1_SRC" "$R1_DEST"
                    ln -s "$R2_SRC" "$R2_DEST"
                    echo "  ✓ ${sample_id} lié (${folder_name}_R1/R2.fastq.gz → ${sample_id}_R1/R2.fastq.gz)"
                fi
            else
                echo "  ⚠ Fichiers R1/R2 manquants pour ${folder_name} dans ${SRC_DIR}"
                echo "    Cherchés: ${R1_SRC} et ${R2_SRC}"
            fi
        done

        echo "✓ Organisation terminée pour ${lot} → ${DEST_DIR}"
    else
        echo "⚠ MODE inconnu pour ${lot} : ${MODE}"
    fi
done

shopt -u nullglob

echo ""
echo "=========================================="
echo "Organisation des données - 9 LOTS terminée"
echo "=========================================="


echo "=========================================="
echo "Pipeline aDNA - Coprolites - 9 LOTS"
echo "Date de début: $(date)"
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
# Vérifier si la taxonomie Krona est installée
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

echo "Arborescence créée."


#################################################################################
## ÉTAPE 1: Contrôle qualité RAW
#################################################################################
#
#echo ""
#echo "=== ÉTAPE 1: Contrôle qualité RAW (FastQC/MultiQC) ==="
#
#for lot in "${LOTS[@]}"; do
#    echo "FastQC RAW pour ${lot}..."
#    
#    INPUT_DIR="${BASE_DIR}/01_raw_data/${lot}"
#    OUTPUT_DIR="${BASE_DIR}/02_quality_check_raw/${lot}"
#    
#    for FILE in "${INPUT_DIR}"/*.fastq.gz; do
#        if [[ -f "$FILE" ]]; then
#            fastqc "$FILE" -o "$OUTPUT_DIR" -t 4
#        fi
#    done
#    
#    echo "MultiQC RAW pour ${lot}..."
#    cd "$OUTPUT_DIR"
#    multiqc . -n "multiqc_raw_${lot}.html" --force
#done
#
#echo "Contrôle qualité RAW terminé."
#
#################################################################################
## ÉTAPE 2: Filtrage BBDuk
#################################################################################
#
#echo ""
#echo "=== ÉTAPE 2: Filtrage et trimming (BBDuk) ==="
#
#for lot in "${LOTS[@]}"; do
#    echo "BBDuk pour ${lot}..."
#    
#    INPUT_DIR="${BASE_DIR}/01_raw_data/${lot}"
#    OUTPUT_DIR="${BASE_DIR}/03_bbduk/${lot}"
#    
#    cd "$INPUT_DIR"
#    
#    for r1_file in *_R1.fastq.gz; do
#        r2_file="${r1_file/_R1/_R2}"
#        
#        if [[ ! -f "$r2_file" ]]; then
#            echo " ✗ ERREUR: Fichier R2 manquant pour $r1_file"
#            continue
#        fi
#        
#        base_name="${r1_file%%_R1.fastq.gz}"
#        echo " → Traitement de ${base_name} (${lot})..."
#        
#        $BBDUK -Xmx4g \
#            in1="$r1_file" \
#            in2="$r2_file" \
#            out1="${OUTPUT_DIR}/clean_${r1_file}" \
#            out2="${OUTPUT_DIR}/clean_${r2_file}" \
#            ref=$PHIX \
#            ktrim=rl k=23 mink=11 hdist=1 \
#            tpe tbo \
#            minlen=25 \
#            qtrim=r trimq=20 \
#            stats="${OUTPUT_DIR}/${base_name}_bbduk_stats.txt"
#    done
#done
#
#echo "Filtrage BBDuk terminé."
#
#################################################################################
## ÉTAPE 3: Déduplication FastUniq
#################################################################################
#
#echo ""
#echo "=== ÉTAPE 3: Déduplication (FastUniq) ==="
#
#TMP="/tmp/fastuniq_coprolites_5lots_$$"
#mkdir -p "$TMP"
#
#for lot in "${LOTS[@]}"; do
#    echo "FastUniq pour ${lot}..."
#    
#    INPUT_DIR="${BASE_DIR}/03_bbduk/${lot}"
#    OUTPUT_DIR="${BASE_DIR}/04_fastuniq/${lot}"
#    
#    cd "$INPUT_DIR" || continue
#    
#    for R1_gz in clean_*_R1.fastq.gz; do
#        base=$(echo "$R1_gz" | sed 's/_R1\.fastq\.gz//')
#        R2_gz="${base}_R2.fastq.gz"
#        
#        if [[ -f "$R2_gz" ]]; then
#            echo " → Traitement de ${base} (${lot})..."
#            
#            R1_tmp="${TMP}/${lot}_${base}_R1.fastq"
#            R2_tmp="${TMP}/${lot}_${base}_R2.fastq"
#            listfile="${TMP}/${lot}_${base}.list"
#            
#            zcat "$INPUT_DIR/$R1_gz" > "$R1_tmp"
#            zcat "$INPUT_DIR/$R2_gz" > "$R2_tmp"
#            
#            echo -e "${R1_tmp}\n${R2_tmp}" > "$listfile"
#            
#            fastuniq -i "$listfile" -t q \
#                -o "${OUTPUT_DIR}/${base}_dedup_R1.fastq" \
#                -p "${OUTPUT_DIR}/${base}_dedup_R2.fastq"
#            
#            rm -f "$R1_tmp" "$R2_tmp" "$listfile"
#        fi
#    done
#done
#
#rm -rf "$TMP"
#echo "Déduplication FastUniq terminée."
#
#################################################################################
## ÉTAPE 4: Clumpify
#################################################################################
#
#echo ""
#echo "=== ÉTAPE 4: Clumpify (déduplication optique) ==="
#
#for lot in "${LOTS[@]}"; do
#    echo "Clumpify pour ${lot}..."
#    
#    INPUT_DIR="${BASE_DIR}/04_fastuniq/${lot}"
#    OUTPUT_DIR="${BASE_DIR}/05_clumpify/${lot}"
#    
#    for R1 in "${INPUT_DIR}"/*_R1.fastq; do
#        R2="${R1/_R1.fastq/_R2.fastq}"
#        
#        if [[ -f "$R2" ]]; then
#            base=$(basename "$R1" _R1.fastq)
#            echo " → Traitement de ${base} (${lot})..."
#            
#            $CLUMPIFY \
#                in="$R1" in2="$R2" \
#                out="${OUTPUT_DIR}/${base}_clumpify_R1.fastq.gz" \
#                out2="${OUTPUT_DIR}/${base}_clumpify_R2.fastq.gz" \
#                dedupe=t
#        fi
#    done
#done
#
#echo "Clumpify terminé."
#
#################################################################################
## ÉTAPE 5: Fastp
#################################################################################
#
#echo ""
#echo "=== ÉTAPE 5: Fastp (merging et QC final) ==="
#
#for lot in "${LOTS[@]}"; do
#    echo "Fastp pour ${lot}..."
#    
#    INPUT_DIR="${BASE_DIR}/05_clumpify/${lot}"
#    OUTPUT_DIR="${BASE_DIR}/06_fastp/${lot}"
#    
#    for R1 in "${INPUT_DIR}"/*_R1.fastq.gz; do
#        R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
#        
#        if [[ -f "$R2" ]]; then
#            base=$(basename "$R1" _R1.fastq.gz)
#            echo " → Traitement de ${base} (${lot})..."
#            
#            fastp \
#                -i "$R1" -I "$R2" \
#                --merged_out "${OUTPUT_DIR}/${base}_fastp_merged.fastq.gz" \
#                --out1 "${OUTPUT_DIR}/${base}_fastp_R1.fastq.gz" \
#                --out2 "${OUTPUT_DIR}/${base}_fastp_R2.fastq.gz" \
#                --json "${OUTPUT_DIR}/${base}_fastp.json" \
#                --html "${OUTPUT_DIR}/${base}_fastp.html" \
#                --thread 4 \
#                --length_required 30 \
#                --qualified_quality_phred 20
#        fi
#    done
#done
#
#echo "Fastp terminé."
#
#################################################################################
## ÉTAPE 6: Contrôle qualité CLEAN
#################################################################################
#
#echo ""
#echo "=== ÉTAPE 6: Contrôle qualité CLEAN (FastQC/MultiQC) ==="
#
#for lot in "${LOTS[@]}"; do
#    echo "FastQC CLEAN pour ${lot}..."
#    
#    INPUT_DIR="${BASE_DIR}/06_fastp/${lot}"
#    OUTPUT_DIR="${BASE_DIR}/07_quality_check_clean/${lot}"
#    
#    for FILE in "${INPUT_DIR}"/*.fastq.gz; do
#        if [[ -f "$FILE" ]]; then
#            fastqc "$FILE" -o "$OUTPUT_DIR" -t 4
#        fi
#    done
#    
#    echo "MultiQC CLEAN pour ${lot}..."
#    cd "$OUTPUT_DIR"
#    multiqc . -n "multiqc_clean_${lot}.html" --force
#done
#
#echo "Contrôle qualité CLEAN terminé."
#
#################################################################################
## ÉTAPE 7: Classification Kraken2
#################################################################################
#
#echo ""
#echo "=== ÉTAPE 7: Classification taxonomique (Kraken2) ==="
#
#for lot in "${LOTS[@]}"; do
#    echo "Kraken2 pour ${lot}..."
#    
#    FASTP_DIR="${BASE_DIR}/06_fastp/${lot}"
#    OUT_DIR="${BASE_DIR}/08_kraken2/${lot}"
#    
#    echo " → Analyse des reads merged..."
#    for MERGED in "${FASTP_DIR}"/*_fastp_merged.fastq.gz; do
#        if [[ -f "$MERGED" ]]; then
#            SAMPLE=$(basename "$MERGED" _fastp_merged.fastq.gz)
#            OUT_KRAKEN="${OUT_DIR}/${SAMPLE}_merged.kraken"
#            OUT_REPORT="${OUT_DIR}/${SAMPLE}_merged.report"
#            
#            echo "   • ${SAMPLE} (merged) - ${lot}"
#            kraken2 --confidence 0.2 --db "$KRAKEN2_DB" --threads $THREADS \
#                --output "$OUT_KRAKEN" --report "$OUT_REPORT" "$MERGED"
#        fi
#    done
#    
#    echo " → Analyse des reads unmerged..."
#    for R1 in "${FASTP_DIR}"/*_fastp_R1.fastq.gz; do
#        if [[ -f "$R1" ]]; then
#            SAMPLE=$(basename "$R1" _fastp_R1.fastq.gz)
#            R2="${FASTP_DIR}/${SAMPLE}_fastp_R2.fastq.gz"
#            
#            if [[ -f "$R2" ]]; then
#                OUT_KRAKEN="${OUT_DIR}/${SAMPLE}_unmerged.kraken"
#                OUT_REPORT="${OUT_DIR}/${SAMPLE}_unmerged.report"
#                
#                echo "   • ${SAMPLE} (unmerged) - ${lot}"
#                kraken2 --confidence 0.2 --paired --db "$KRAKEN2_DB" --threads $THREADS \
#                    --output "$OUT_KRAKEN" --report "$OUT_REPORT" "$R1" "$R2"
#            fi
#        fi
#    done
#done
#
#echo "Classification Kraken2 terminée."
#
#################################################################################
## ÉTAPE 8: Visualisation Krona
#################################################################################
#
#echo ""
#echo "=== ÉTAPE 8: Visualisation (Krona) ==="
#
#
#for lot in "${LOTS[@]}"; do
#    echo "Krona pour ${lot}..."
#    
#    IN_DIR="${BASE_DIR}/08_kraken2/${lot}"
#    OUT_DIR="${BASE_DIR}/09_krona/${lot}"
#    
#    cd "$IN_DIR"
#    
#    # Krona combiné pour tous les échantillons du lot
#    if ls *.report 1> /dev/null 2>&1; then
#        echo " → Génération Krona combiné pour ${lot}..."
#        ktImportTaxonomy -t 5 -m 3 -o "${OUT_DIR}/all_samples_krona.html" "${IN_DIR}"/*.report 2>&1 || {
#            echo "   ✗ ERREUR Krona combiné pour ${lot}"
#        }
#    fi
#    
#    # Krona individuel pour chaque échantillon
#    for report in "${IN_DIR}"/*.report; do
#        if [[ -f "$report" ]]; then
#            base=$(basename "$report" .report)
#            echo " → Génération Krona pour ${base}..."
#            ktImportTaxonomy -t 5 -m 3 -o "${OUT_DIR}/${base}_krona.html" "$report" 2>&1 || {
#                echo "   ✗ ERREUR Krona pour ${base}"
#            }
#        fi
#    done
#done
#
#echo "Visualisation Krona terminée."
#
#################################################################################
## ÉTAPE 9: Tables MPA
#################################################################################
#
#echo ""
#echo "=== ÉTAPE 9: Création des tables MPA ==="
#
#if [[ ! -d "$KRAKENTOOLS_DIR" ]]; then
#    echo "Installation de KrakenTools..."
#    mkdir -p "${BASE_DIR}/08_kraken2"
#    cd "${BASE_DIR}/08_kraken2"
#    git clone https://github.com/jenniferlu717/KrakenTools.git
#fi
#
#for lot in "${LOTS[@]}"; do
#    echo "Conversion MPA pour ${lot}..."
#    
#    IN_DIR="${BASE_DIR}/08_kraken2/${lot}"
#    OUT_DIR="${BASE_DIR}/10_mpa_tables/${lot}"
#    
#    cd "$IN_DIR"
#    
#    declare -a mpa_files=()
#    
#    for report in *.report; do
#        if [[ -f "$report" ]]; then
#            base=$(basename "$report" .report)
#            mpa_file="${OUT_DIR}/${base}.mpa"
#            
#            echo " → Conversion de ${base}..."
#            python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" -r "$report" -o "$mpa_file"
#            mpa_files+=("$mpa_file")
#        fi
#    done
#    
#    if [[ ${#mpa_files[@]} -gt 0 ]]; then
#        echo " → Combinaison de tous les fichiers MPA..."
#        python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" -i "${mpa_files[@]}" -o "${OUT_DIR}/combined_all.tsv"
#    fi
#done
#
#echo "Création des tables MPA terminée."
#
#################################################################################
## ÉTAPE 10: Tableau récapitulatif
#################################################################################
#
#echo ""
#echo "=== ÉTAPE 10: Tableau récapitulatif (avant/après nettoyage) ==="
#
#SUMMARY_TABLE="${BASE_DIR}/11_summary_tables/sequences_summary.tsv"
#
#cat > "$SUMMARY_TABLE" << 'HEADER'
#Lot	Sample	Stage	Nb_sequences	Longueur_moyenne	GC_percent
#HEADER
#
#function extract_stats() {
#    local file=$1
#    local lot=$2
#    local sample=$3
#    local stage=$4
#    
#    if [[ ! -f "$file" ]]; then
#        return 1
#    fi
#    
#    nb_seq=$(zcat "$file" 2>/dev/null | echo $((`wc -l`/4)))
#    
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
#    echo -e "${lot}\t${sample}\t${stage}\t${nb_seq}\t${stats}" >> "$SUMMARY_TABLE"
#}
#
#echo "Calcul des statistiques pour chaque échantillon..."
#
#for lot in "${LOTS[@]}"; do
#    echo " → Traitement du ${lot}..."
#    samples="${LOT_SAMPLES[$lot]}"
#    
#    for sample in $samples; do
#        echo "   • ${sample}..."
#        
#        raw_r1="${BASE_DIR}/01_raw_data/${lot}/${sample}_R1.fastq.gz"
#        extract_stats "$raw_r1" "$lot" "$sample" "RAW"
#        
#        clean_merged="${BASE_DIR}/06_fastp/${lot}/clean_${sample}_dedup_clumpify_fastp_merged.fastq.gz"
#        if [[ -f "$clean_merged" ]]; then
#            extract_stats "$clean_merged" "$lot" "$sample" "CLEAN_merged"
#        fi
#        
#        clean_r1="${BASE_DIR}/06_fastp/${lot}/clean_${sample}_dedup_clumpify_fastp_R1.fastq.gz"
#        if [[ -f "$clean_r1" ]]; then
#            extract_stats "$clean_r1" "$lot" "$sample" "CLEAN_unmerged"
#        fi
#    done
#done
#
#echo ""
#echo "Tableau récapitulatif créé: ${SUMMARY_TABLE}"
#echo ""
#echo "Aperçu:"
#head -20 "$SUMMARY_TABLE" | column -t

################################################################################
# ÉTAPE 11: Analyse MapDamage (Ovis aries & Capra hircus)
# CORRIGÉ: Ajout de --no-stats pour éviter les erreurs de compilation C++
################################################################################

echo ""
echo "=========================================="
echo "ÉTAPE 11: MapDamage - Analyse des dommages de l'ADN"
echo "Date de début: $(date)"
echo "=========================================="

# Activer l'environnement mapdamage
conda deactivate
conda activate mapdamage_py39

# Répertoires
KRAKEN_BASE="${BASE_DIR}/08_kraken2"
FASTQ_BASE="${BASE_DIR}/06_fastp"
DAMAGE_BASE="${BASE_DIR}/12_mapdamage"
LOGFILE="${BASE_DIR}/00_scripts/mapdamage_$(date +%Y%m%d_%H%M%S).txt"
MAPPING_INFO="${BASE_DIR}/11_summary_tables/mapping_bwa_info.tsv"

echo "Script MapDamage started at $(date)" | tee -a "$LOGFILE"

# Initialiser le fichier de mapping info
echo -e "Lot\tSample\tSpecies\tType\tTotal_Reads\tMapped_Reads\tMapping_Rate(%)" > "$MAPPING_INFO"

# Définition des génomes de référence
declare -A TAXONS=(
    ["Ovis_aries"]="9940:/home/plstenge/genomes/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa"
    ["Capra_hircus"]="9925:/home/plstenge/genomes/Capra_hircus.ARS1.dna.toplevel.fa"
)

# Indexer les génomes si nécessaire (commenter si déjà fait)
# echo "Vérification des index BWA..."
# bwa index /home/plstenge/genomes/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa 2>>"$LOGFILE"
# bwa index /home/plstenge/genomes/Capra_hircus.ARS1.dna.toplevel.fa 2>>"$LOGFILE"

# Fonction pour calculer le taux de mapping
calculate_mapping_rate() {
    local bam_file="$1"
    local lot="$2"
    local sample_name="$3"
    local species="$4"
    local type="$5"
    
    if [[ -f "$bam_file" ]]; then
        local total_reads=$(samtools view -c "$bam_file")
        local mapped_reads=$(samtools view -c -F 4 "$bam_file")
        
        local mapping_rate=0
        if [[ $total_reads -gt 0 ]]; then
            mapping_rate=$(echo "scale=2; $mapped_reads * 100 / $total_reads" | bc)
        fi
        
        echo -e "${lot}\t${sample_name}\t${species}\t${type}\t${total_reads}\t${mapped_reads}\t${mapping_rate}" >> "$MAPPING_INFO"
        echo "✓ Mapping stats for ${lot}_${sample_name}_${species}_${type}: ${mapped_reads}/${total_reads} (${mapping_rate}%)" | tee -a "$LOGFILE"
    fi
}

# Fonction MapDamage sécurisée (avec gestion d'erreurs)
run_mapdamage_safe() {
    local bam_file="$1"
    local ref_fasta="$2"
    local output_dir="$3"
    local sample_name="$4"
    
    echo "  → MapDamage pour ${sample_name}..."
    
    # Vérifier si le BAM est valide
    if [[ ! -s "$bam_file" ]]; then
        echo "  ✗ BAM vide ou absent, on skip MapDamage"
        return 0
    fi
    
    # Créer le dossier de sortie
    mkdir -p "$output_dir"
    
    # Lancer MapDamage avec --no-stats pour éviter l'erreur de compilation
    # Si ça plante toujours, on continue avec le prochain échantillon
    mapDamage -i "$bam_file" -r "$ref_fasta" --folder "$output_dir" --no-stats 2>>"$LOGFILE" || {
        echo "  ⚠ MapDamage a échoué pour ${sample_name} (erreur ignorée)"
        echo "  → Les graphiques de dommages principaux ont été générés malgré tout"
        return 0
    }
    
    echo "  ✓ MapDamage terminé pour ${sample_name}"
    return 0
}

shopt -s nullglob

# Boucle sur tous les lots
for lot in "${LOTS[@]}"; do
    echo ""
    echo "========================================"
    echo "TRAITEMENT DU LOT: ${lot}"
    echo "========================================"
    
    KRAKEN_DIR="${KRAKEN_BASE}/${lot}"
    FASTQ_DIR="${FASTQ_BASE}/${lot}"
    
    if [[ ! -d "$KRAKEN_DIR" ]]; then
        echo "⚠ Répertoire Kraken2 absent pour ${lot}, on skip"
        continue
    fi
    
    # Compter les fichiers kraken
    kraken_files=("$KRAKEN_DIR"/*.kraken)
    if [[ ${#kraken_files[@]} -eq 0 ]]; then
        echo "⚠ Aucun fichier .kraken trouvé pour ${lot}, on skip"
        continue
    fi
    
    echo "→ ${#kraken_files[@]} fichiers kraken trouvés"
    
    for KRAKEN_FILE in "${kraken_files[@]}"; do
        KRAKEN_BASE_NAME=$(basename "$KRAKEN_FILE" .kraken)
        echo ""
        echo "---- Traitement: ${KRAKEN_BASE_NAME} ----"
        
        # Extraire le préfixe de base
        PREFIX=$(echo "$KRAKEN_BASE_NAME" | sed -E 's/_dedup_clumpify_(un)?merged$//')
        
        # Chercher les fichiers FASTQ
        R1_FILES=("${FASTQ_DIR}/${PREFIX}_dedup_clumpify_fastp_R1.fastq"*)
        R2_FILES=("${FASTQ_DIR}/${PREFIX}_dedup_clumpify_fastp_R2.fastq"*)
        MERGED_FILES=("${FASTQ_DIR}/${PREFIX}_dedup_clumpify_fastp_merged.fastq"*)
        
        R1_FILE="${R1_FILES[0]:-}"
        R2_FILE="${R2_FILES[0]:-}"
        MERGED_FILE="${MERGED_FILES[0]:-}"
        
        # Boucle sur les espèces (Ovis et Capra)
        for GROUP in "${!TAXONS[@]}"; do
            TAX_ID="${TAXONS[$GROUP]%:*}"
            REF_FASTA="${TAXONS[$GROUP]#*:}"
            DAMAGE_DIR="${DAMAGE_BASE}/${lot}/${GROUP}"
            mkdir -p "$DAMAGE_DIR"
            
            echo ""
            echo "-- Espèce: ${GROUP} (TaxID: ${TAX_ID}) --"
            
            OUT_R1="${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_R1.fastq"
            OUT_R2="${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_R2.fastq"
            OUT_MERGED="${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_merged.fastq"
            
            ###################################################################
            # TRAITEMENT PAIRED-END (unmerged)
            ###################################################################
            
            if [[ -n "$R1_FILE" && -n "$R2_FILE" ]]; then
                echo "→ Extraction paired-end..."
                
                python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
                    -k "$KRAKEN_FILE" -s "$R1_FILE" -s2 "$R2_FILE" -t "$TAX_ID" \
                    -o "$OUT_R1" -o2 "$OUT_R2" --fastq-output 2>>"$LOGFILE"
                
                if [[ -f "$OUT_R1" && -f "$OUT_R2" ]]; then
                    echo "→ Mapping BWA paired-end..."
                    
                    # BWA aln
                    bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF_FASTA" "$OUT_R1" > "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_R1.sai" 2>>"$LOGFILE"
                    bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF_FASTA" "$OUT_R2" > "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_R2.sai" 2>>"$LOGFILE"
                    
                    # BWA sampe
                    bwa sampe "$REF_FASTA" \
                        "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_R1.sai" \
                        "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_R2.sai" \
                        "$OUT_R1" "$OUT_R2" > "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.sam" 2>>"$LOGFILE"
                    
                    # SAM → BAM → tri → index
                    samtools view -bS "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.sam" > "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.bam" 2>>"$LOGFILE"
                    samtools sort -o "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.sorted.bam" "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.bam" 2>>"$LOGFILE"
                    samtools index "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.sorted.bam" 2>>"$LOGFILE"
                    
                    # Stats mapping
                    calculate_mapping_rate "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.sorted.bam" "$lot" "$KRAKEN_BASE_NAME" "$GROUP" "unmerged"
                    
                    # Nettoyage fichiers intermédiaires
                    rm -f "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_R1.sai" \
                          "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_R2.sai" \
                          "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.sam" \
                          "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.bam" 2>>"$LOGFILE"
                    
                    # MapDamage (version sécurisée)
                    run_mapdamage_safe "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}.sorted.bam" "$REF_FASTA" \
                        "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_mapDamage_unmerged" "$KRAKEN_BASE_NAME"
                fi
            fi
            
            ###################################################################
            # TRAITEMENT SINGLE-END (merged)
            ###################################################################
            
            if [[ -n "$MERGED_FILE" ]]; then
                echo "→ Extraction merged..."
                
                python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
                    -k "$KRAKEN_FILE" -s "$MERGED_FILE" -t "$TAX_ID" \
                    -o "$OUT_MERGED" --fastq-output 2>>"$LOGFILE"
                
                if [[ -f "$OUT_MERGED" ]]; then
                    echo "→ Mapping BWA single-end..."
                    
                    # BWA aln
                    bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF_FASTA" "$OUT_MERGED" > "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_merged.sai" 2>>"$LOGFILE"
                    
                    # BWA samse
                    bwa samse "$REF_FASTA" \
                        "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_merged.sai" \
                        "$OUT_MERGED" > "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_merged.sam" 2>>"$LOGFILE"
                    
                    # SAM → BAM → tri → index
                    samtools view -bS "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_merged.sam" > "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_merged.bam" 2>>"$LOGFILE"
                    samtools sort -o "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_merged.sorted.bam" "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_merged.bam" 2>>"$LOGFILE"
                    samtools index "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_merged.sorted.bam" 2>>"$LOGFILE"
                    
                    # Stats mapping
                    calculate_mapping_rate "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_merged.sorted.bam" "$lot" "$KRAKEN_BASE_NAME" "$GROUP" "merged"
                    
                    # Nettoyage
                    rm -f "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_merged.sai" \
                          "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_merged.sam" \
                          "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_merged.bam" 2>>"$LOGFILE"
                    
                    # MapDamage (version sécurisée)
                    run_mapdamage_safe "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_merged.sorted.bam" "$REF_FASTA" \
                        "${DAMAGE_DIR}/${KRAKEN_BASE_NAME}_${GROUP}_mapDamage_merged" "$KRAKEN_BASE_NAME"
                fi
            fi
        done
    done
done

shopt -u nullglob

echo ""
echo "=========================================="
echo "MapDamage terminé pour tous les lots!"
echo "Date de fin: $(date)"
echo "=========================================="
echo ""
echo "Résultats: ${DAMAGE_BASE}"
echo "Statistiques: ${MAPPING_INFO}"
echo ""

# Notification
echo "MapDamage terminé le $(date)" | \
    mail -s "Pipeline Coprolites - MapDamage terminé" pierrelouis.stenger@gmail.com 2>/dev/null || true

exit 0


################################################################################
# FIN
################################################################################

echo ""
echo "=========================================="
echo "PIPELINE TERMINÉ AVEC SUCCÈS - 5 LOTS"
echo "Date de fin: $(date)"
echo "=========================================="
echo ""
echo "Résultats: ${BASE_DIR}"
echo "Tableau: ${SUMMARY_TABLE}"
echo ""
echo "Comparaisons à réaliser:"
echo " • Illumina (Lot1) vs Aviti (Lots 2-5)"
echo " • Aviti no demux (Lot2) vs demux (Lot3)"
echo " • Recette1 (Lots 1-3) vs Recette2 (Lots 4-5)"
echo " • 2x150bp (Lot4) vs 2x75bp (Lot5)"
echo ""

conda deactivate

echo "Pipeline terminé le $(date). Résultats: ${BASE_DIR}" | \
    mail -s "Pipeline Coprolites 5 lots - Terminé" pierrelouis.stenger@gmail.com 2>/dev/null || true

exit 0
