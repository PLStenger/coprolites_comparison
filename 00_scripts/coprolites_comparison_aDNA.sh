#!/bin/bash

#SBATCH --job-name=coprolites_aDNA
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/coprolites_comparison_aDNA.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/coprolites_comparison_aDNA.out"

################################################################################
# Pipeline d'analyse aDNA - Projet Coprolites
# Comparaison de recettes de séquençage
# Author: Pierre-Louis Stenger
# Date: November 2025
#
# Description:
# Ce pipeline analyse des échantillons de coprolites séquencés avec 2 recettes:
# - Recette1 (R1): recipes-ShortInsert (2x150bp) - Aviti uniquement
# - Recette2 (R2): recipes-v3.3.x-ShortInsert-SmallFragmentRecovery (2x150bp)
#   → Aviti + Illumina
#
# Échantillons analysés:
# - cop408, cop410, cop412, cop414, cop417
#
# Technologies:
# - Recette1: Aviti uniquement
# - Recette2: Aviti + Illumina (analysés séparément et combinés)
################################################################################

set -eo pipefail

################################################################################
# CONFIGURATION GLOBALE
################################################################################

# Chemins de base
BASE_DIR="/home/plstenge/coprolites_comparison"
RECIPE1_BASE="/storage/groups/gdec/shared_paleo/E1531_final/run3_20251008_AV241601_E1531_Ps7_Ps8"
RECIPE2_BASE="/home/plstenge/coprolites_comparison/01_raw_data"

# Outils BBMap (chemins absolus)
BBDUK="/home/plstenge/bbmap/bbduk.sh"
CLUMPIFY="/home/plstenge/bbmap/clumpify.sh"
PHIX="/home/plstenge/bbmap/resources/phix174_ill.ref.fa.gz"

# Base de données Kraken2
KRAKEN2_DB="/home/plstenge/k2_core_nt_20250609"

# KrakenTools
KRAKENTOOLS_DIR="${BASE_DIR}/07_kraken2/KrakenTools"

# Paramètres
THREADS=36

# Définition des échantillons
declare -a SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")

# Définition des IDs de dossiers pour Recette1 (Aviti)
declare -A RECIPE1_FOLDERS=(
    ["cop408"]="474_cop408"
    ["cop410"]="475_cop410"
    ["cop412"]="476_cop412"
    ["cop414"]="477_cop414"
    ["cop417"]="478_cop417"
)

echo "=========================================="
echo "Pipeline aDNA - Coprolites"
echo "Date de début: $(date)"
echo "=========================================="

################################################################################
# ACTIVATION DE L'ENVIRONNEMENT CONDA UNIQUE
################################################################################

echo ""
echo "=== Activation de l'environnement conda unique ==="
echo ""

# Chargement du module conda
module load conda/4.12.0
source ~/.bashrc

# Activation de l'environnement unique contenant tous les outils
# IMPORTANT: Cet environnement doit être créé au préalable avec:
# conda create -n coprolites-pipeline -c bioconda -c conda-forge \
#   fastqc multiqc bbmap fastuniq fastp kraken2 krona python=3.9
# pip install krakentools (ou installation manuelle des scripts python)

conda activate coprolites-pipeline

echo "Environnement conda activé: coprolites-pipeline"
echo "Tous les outils seront utilisés depuis cet environnement unique."

################################################################################
# ÉTAPE 0: Création de l'arborescence et organisation des données
################################################################################

echo ""
echo "=== ÉTAPE 0: Préparation de l'arborescence ==="
echo ""

# Création de la structure de dossiers
mkdir -p "${BASE_DIR}/00_scripts"
mkdir -p "${BASE_DIR}/01_raw_data/recette1_aviti"
mkdir -p "${BASE_DIR}/01_raw_data/recette2_aviti"
mkdir -p "${BASE_DIR}/01_raw_data/recette2_illumina"
mkdir -p "${BASE_DIR}/01_raw_data/recette2_aviti_illumina_combined"
mkdir -p "${BASE_DIR}/02_quality_check"
mkdir -p "${BASE_DIR}/03_bbduk"
mkdir -p "${BASE_DIR}/04_fastuniq"
mkdir -p "${BASE_DIR}/05_clumpify"
mkdir -p "${BASE_DIR}/06_fastp"
mkdir -p "${BASE_DIR}/07_kraken2"
mkdir -p "${BASE_DIR}/08_krona"
mkdir -p "${BASE_DIR}/09_mpa_tables"

echo "Arborescence créée."

################################################################################
# Copie/liens symboliques des données Recette 1 (Aviti uniquement)
################################################################################

echo ""
echo "=== Copie des données Recette1 (Aviti) ==="
echo ""

for sample in "${SAMPLES[@]}"; do
    folder="${RECIPE1_FOLDERS[$sample]}"
    source_dir="${RECIPE1_BASE}/${folder}"
    
    R1_source="${source_dir}/${folder}_R1.fastq.gz"
    R2_source="${source_dir}/${folder}_R2.fastq.gz"
    
    R1_dest="${BASE_DIR}/01_raw_data/recette1_aviti/${sample}_R1.fastq.gz"
    R2_dest="${BASE_DIR}/01_raw_data/recette1_aviti/${sample}_R2.fastq.gz"
    
    if [[ -f "$R1_source" ]]; then
        ln -sf "$R1_source" "$R1_dest"
        echo " ✓ ${sample}_R1 lié (Recette1 Aviti)"
    else
        echo " ✗ ATTENTION: ${sample}_R1 non trouvé dans ${source_dir}"
    fi
    
    if [[ -f "$R2_source" ]]; then
        ln -sf "$R2_source" "$R2_dest"
        echo " ✓ ${sample}_R2 lié (Recette1 Aviti)"
    else
        echo " ✗ ATTENTION: ${sample}_R2 non trouvé dans ${source_dir}"
    fi
done

################################################################################
# Organisation des données Recette 2 (Aviti et Illumina)
################################################################################

echo ""
echo "=== Organisation des données Recette2 (Aviti + Illumina) ==="
echo ""

# Les fichiers sont déjà dans /home/plstenge/coprolites_comparison/01_raw_data
# On crée des liens symboliques vers les sous-dossiers appropriés

for sample in "${SAMPLES[@]}"; do
    # Aviti Recette2
    aviti_r1_source="${RECIPE2_BASE}/Aviti_${sample}_R1.fastq.gz"
    aviti_r2_source="${RECIPE2_BASE}/Aviti_${sample}_R2.fastq.gz"
    
    if [[ -f "$aviti_r1_source" ]]; then
        ln -sf "$aviti_r1_source" "${BASE_DIR}/01_raw_data/recette2_aviti/${sample}_R1.fastq.gz"
        echo " ✓ ${sample} Aviti_R1 lié (Recette2)"
    fi
    
    if [[ -f "$aviti_r2_source" ]]; then
        ln -sf "$aviti_r2_source" "${BASE_DIR}/01_raw_data/recette2_aviti/${sample}_R2.fastq.gz"
        echo " ✓ ${sample} Aviti_R2 lié (Recette2)"
    fi
    
    # Illumina Recette2
    illumina_r1_source="${RECIPE2_BASE}/Illumina_${sample}_R1.fastq.gz"
    illumina_r2_source="${RECIPE2_BASE}/Illumina_${sample}_R2.fastq.gz"
    
    if [[ -f "$illumina_r1_source" ]]; then
        ln -sf "$illumina_r1_source" "${BASE_DIR}/01_raw_data/recette2_illumina/${sample}_R1.fastq.gz"
        echo " ✓ ${sample} Illumina_R1 lié (Recette2)"
    fi
    
    if [[ -f "$illumina_r2_source" ]]; then
        ln -sf "$illumina_r2_source" "${BASE_DIR}/01_raw_data/recette2_illumina/${sample}_R2.fastq.gz"
        echo " ✓ ${sample} Illumina_R2 lié (Recette2)"
    fi
done

################################################################################
# Création des fichiers combinés Aviti+Illumina pour Recette2
################################################################################

echo ""
echo "=== Fusion Aviti + Illumina (Recette2) ==="
echo ""

for sample in "${SAMPLES[@]}"; do
    aviti_r1="${BASE_DIR}/01_raw_data/recette2_aviti/${sample}_R1.fastq.gz"
    aviti_r2="${BASE_DIR}/01_raw_data/recette2_aviti/${sample}_R2.fastq.gz"
    illumina_r1="${BASE_DIR}/01_raw_data/recette2_illumina/${sample}_R1.fastq.gz"
    illumina_r2="${BASE_DIR}/01_raw_data/recette2_illumina/${sample}_R2.fastq.gz"
    
    combined_r1="${BASE_DIR}/01_raw_data/recette2_aviti_illumina_combined/${sample}_R1.fastq.gz"
    combined_r2="${BASE_DIR}/01_raw_data/recette2_aviti_illumina_combined/${sample}_R2.fastq.gz"
    
    # Fusion R1
    if [[ -f "$aviti_r1" ]] && [[ -f "$illumina_r1" ]]; then
        cat "$aviti_r1" "$illumina_r1" > "$combined_r1"
        echo " ✓ ${sample}_R1 fusionné (Aviti + Illumina)"
    elif [[ -f "$aviti_r1" ]]; then
        ln -sf "$aviti_r1" "$combined_r1"
        echo " → ${sample}_R1 Aviti seul (pas d'Illumina)"
    elif [[ -f "$illumina_r1" ]]; then
        ln -sf "$illumina_r1" "$combined_r1"
        echo " → ${sample}_R1 Illumina seul (pas d'Aviti)"
    fi
    
    # Fusion R2
    if [[ -f "$aviti_r2" ]] && [[ -f "$illumina_r2" ]]; then
        cat "$aviti_r2" "$illumina_r2" > "$combined_r2"
        echo " ✓ ${sample}_R2 fusionné (Aviti + Illumina)"
    elif [[ -f "$aviti_r2" ]]; then
        ln -sf "$aviti_r2" "$combined_r2"
        echo " → ${sample}_R2 Aviti seul (pas d'Illumina)"
    elif [[ -f "$illumina_r2" ]]; then
        ln -sf "$illumina_r2" "$combined_r2"
        echo " → ${sample}_R2 Illumina seul (pas d'Aviti)"
    fi
done

echo "Organisation des données terminée."

################################################################################
# ÉTAPE 1: Contrôle qualité avec FastQC et MultiQC
################################################################################

echo ""
echo "=== ÉTAPE 1: Contrôle qualité (FastQC/MultiQC) ==="
echo ""

# Définition des stratégies à analyser
declare -a STRATEGIES=("recette1_aviti" "recette2_aviti" "recette2_illumina" "recette2_aviti_illumina_combined")

for strategy in "${STRATEGIES[@]}"; do
    echo "FastQC pour ${strategy}..."
    
    INPUT_DIR="${BASE_DIR}/01_raw_data/${strategy}"
    OUTPUT_DIR="${BASE_DIR}/02_quality_check/${strategy}"
    mkdir -p "$OUTPUT_DIR"
    
    for FILE in "${INPUT_DIR}"/*.fastq.gz; do
        if [[ -f "$FILE" ]]; then
            fastqc "$FILE" -o "$OUTPUT_DIR" -t 4
        fi
    done
    
    # MultiQC
    echo "MultiQC pour ${strategy}..."
    cd "$OUTPUT_DIR"
    multiqc . -n "multiqc_${strategy}.html" --force
done

echo "Contrôle qualité terminé."

################################################################################
# ÉTAPE 2: Filtrage et trimming avec BBDuk
################################################################################

echo ""
echo "=== ÉTAPE 2: Filtrage et trimming (BBDuk) ==="
echo ""

for strategy in "${STRATEGIES[@]}"; do
    echo "BBDuk pour ${strategy}..."
    
    INPUT_DIR="${BASE_DIR}/01_raw_data/${strategy}"
    OUTPUT_DIR="${BASE_DIR}/03_bbduk/${strategy}"
    mkdir -p "$OUTPUT_DIR"
    
    cd "$INPUT_DIR"
    
    for r1_file in *_R1.fastq.gz; do
        r2_file="${r1_file/_R1/_R2}"
        
        if [[ ! -f "$r2_file" ]]; then
            echo " ✗ ERREUR: Fichier R2 manquant pour $r1_file" >&2
            continue
        fi
        
        base_name="${r1_file%%_R1.fastq.gz}"
        echo " → Traitement de ${base_name}..."
        
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

TMP="/tmp/fastuniq_coprolites_tmp_$$"
mkdir -p "$TMP"

for strategy in "${STRATEGIES[@]}"; do
    echo "FastUniq pour ${strategy}..."
    
    INPUT_DIR="${BASE_DIR}/03_bbduk/${strategy}"
    OUTPUT_DIR="${BASE_DIR}/04_fastuniq/${strategy}"
    mkdir -p "$OUTPUT_DIR"
    
    cd "$INPUT_DIR" || continue
    
    for R1_gz in clean_*_R1.fastq.gz; do
        base=$(echo "$R1_gz" | sed 's/_R1\.fastq\.gz//')
        R2_gz="${base}_R2.fastq.gz"
        
        if [[ -f "$R2_gz" ]]; then
            echo " → Traitement de ${base}..."
            
            R1_tmp="${TMP}/${base}_R1.fastq"
            R2_tmp="${TMP}/${base}_R2.fastq"
            listfile="${TMP}/${base}.list"
            
            # Décompression
            zcat "$INPUT_DIR/$R1_gz" > "$R1_tmp"
            zcat "$INPUT_DIR/$R2_gz" > "$R2_tmp"
            
            # Création du fichier liste pour FastUniq
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
# ÉTAPE 4: Clumpify - Déduplication optique supplémentaire
################################################################################

echo ""
echo "=== ÉTAPE 4: Clumpify (déduplication optique) ==="
echo ""

for strategy in "${STRATEGIES[@]}"; do
    echo "Clumpify pour ${strategy}..."
    
    INPUT_DIR="${BASE_DIR}/04_fastuniq/${strategy}"
    OUTPUT_DIR="${BASE_DIR}/05_clumpify/${strategy}"
    mkdir -p "$OUTPUT_DIR"
    
    for R1 in "${INPUT_DIR}"/*_R1.fastq; do
        R2="${R1/_R1.fastq/_R2.fastq}"
        
        if [[ -f "$R2" ]]; then
            base=$(basename "$R1" _R1.fastq)
            echo " → Traitement de ${base}..."
            
            $CLUMPIFY \
                in="$R1" in2="$R2" \
                out="${OUTPUT_DIR}/${base}_clumpify_R1.fastq.gz" \
                out2="${OUTPUT_DIR}/${base}_clumpify_R2.fastq.gz" \
                dedupe=t
        else
            echo " ✗ Fichier R2 manquant pour $R1, ignoré."
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

for strategy in "${STRATEGIES[@]}"; do
    echo "Fastp pour ${strategy}..."
    
    INPUT_DIR="${BASE_DIR}/05_clumpify/${strategy}"
    OUTPUT_DIR="${BASE_DIR}/06_fastp/${strategy}"
    mkdir -p "$OUTPUT_DIR"
    
    for R1 in "${INPUT_DIR}"/*_R1.fastq.gz; do
        R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
        
        if [[ -f "$R2" ]]; then
            base=$(basename "$R1" _R1.fastq.gz)
            echo " → Traitement de ${base}..."
            
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
            echo " ✗ Fichier R2 manquant pour $R1, ignoré."
        fi
    done
done

echo "Fastp terminé."

################################################################################
# ÉTAPE 6: Classification taxonomique avec Kraken2
################################################################################

echo ""
echo "=== ÉTAPE 6: Classification taxonomique (Kraken2) ==="
echo ""

for strategy in "${STRATEGIES[@]}"; do
    echo "Kraken2 pour ${strategy}..."
    
    FASTP_DIR="${BASE_DIR}/06_fastp/${strategy}"
    OUT_DIR="${BASE_DIR}/07_kraken2/${strategy}"
    mkdir -p "$OUT_DIR"
    
    # Analyse des reads merged (single-end)
    echo " → Analyse des reads merged..."
    for MERGED in "${FASTP_DIR}"/*_fastp_merged.fastq.gz; do
        if [[ -f "$MERGED" ]]; then
            SAMPLE=$(basename "$MERGED" _fastp_merged.fastq.gz)
            OUT_KRAKEN="${OUT_DIR}/${SAMPLE}_merged.kraken"
            OUT_REPORT="${OUT_DIR}/${SAMPLE}_merged.report"
            
            echo "   • ${SAMPLE} (merged)"
            
            kraken2 --confidence 0.2 --db "$KRAKEN2_DB" --threads $THREADS \
                --output "$OUT_KRAKEN" --report "$OUT_REPORT" "$MERGED"
        fi
    done
    
    # Analyse des reads unmerged (paired-end)
    echo " → Analyse des reads unmerged..."
    for R1 in "${FASTP_DIR}"/*_fastp_R1.fastq.gz; do
        if [[ -f "$R1" ]]; then
            SAMPLE=$(basename "$R1" _fastp_R1.fastq.gz)
            R2="${FASTP_DIR}/${SAMPLE}_fastp_R2.fastq.gz"
            
            if [[ -f "$R2" ]]; then
                OUT_KRAKEN="${OUT_DIR}/${SAMPLE}_unmerged.kraken"
                OUT_REPORT="${OUT_DIR}/${SAMPLE}_unmerged.report"
                
                echo "   • ${SAMPLE} (unmerged)"
                
                kraken2 --confidence 0.2 --paired --db "$KRAKEN2_DB" --threads $THREADS \
                    --output "$OUT_KRAKEN" --report "$OUT_REPORT" "$R1" "$R2"
            fi
        fi
    done
done

echo "Classification Kraken2 terminée."

################################################################################
# ÉTAPE 7: Visualisation avec Krona
################################################################################

echo ""
echo "=== ÉTAPE 7: Visualisation (Krona) ==="
echo ""

for strategy in "${STRATEGIES[@]}"; do
    echo "Krona pour ${strategy}..."
    
    IN_DIR="${BASE_DIR}/07_kraken2/${strategy}"
    OUT_DIR="${BASE_DIR}/08_krona/${strategy}"
    mkdir -p "$OUT_DIR"
    
    cd "$IN_DIR"
    
    # Krona pour tous les échantillons combinés
    if ls *.report 1> /dev/null 2>&1; then
        ktImportTaxonomy -t 5 -m 3 -o "${OUT_DIR}/all_samples_krona.html" "${IN_DIR}"/*.report
    fi
    
    # Krona individuel pour chaque échantillon
    for report in "${IN_DIR}"/*.report; do
        if [[ -f "$report" ]]; then
            base=$(basename "$report" .report)
            ktImportTaxonomy -t 5 -m 3 -o "${OUT_DIR}/${base}_krona.html" "$report"
        fi
    done
done

echo "Visualisation Krona terminée."

################################################################################
# ÉTAPE 8: Création des tables d'assignation taxonomique (format MPA)
################################################################################

echo ""
echo "=== ÉTAPE 8: Création des tables MPA ==="
echo ""

# Installation ou vérification de KrakenTools si nécessaire
if [[ ! -d "$KRAKENTOOLS_DIR" ]]; then
    echo "Installation de KrakenTools..."
    mkdir -p "${BASE_DIR}/07_kraken2"
    cd "${BASE_DIR}/07_kraken2"
    git clone https://github.com/jenniferlu717/KrakenTools.git
fi

for strategy in "${STRATEGIES[@]}"; do
    echo "Conversion MPA pour ${strategy}..."
    
    IN_DIR="${BASE_DIR}/07_kraken2/${strategy}"
    OUT_DIR="${BASE_DIR}/09_mpa_tables/${strategy}"
    mkdir -p "$OUT_DIR"
    
    cd "$IN_DIR"
    
    # Conversion kreport vers format MPA
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
    
    # Combinaison de tous les fichiers MPA en une seule table
    if [[ ${#mpa_files[@]} -gt 0 ]]; then
        echo " → Combinaison de tous les fichiers MPA..."
        python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" -i "${mpa_files[@]}" -o "${OUT_DIR}/combined_all.tsv"
    fi
done

echo "Création des tables MPA terminée."

################################################################################
# ÉTAPE 9: Génération du rapport de comparaison
################################################################################

echo ""
echo "=== ÉTAPE 9: Génération du rapport de comparaison ==="
echo ""

REPORT_FILE="${BASE_DIR}/00_scripts/pipeline_comparison_report.txt"

cat > "$REPORT_FILE" << 'EOF'
================================================================================
RAPPORT D'ANALYSE - PROJET COPROLITES
================================================================================
Date d'analyse: $(date)
Pipeline version: 2.0 (environnement conda unique)

================================================================================
DESCRIPTION DU PROJET
================================================================================

Analyse comparative de deux recettes de séquençage pour échantillons de coprolites:

Échantillons:
- cop408
- cop410
- cop412
- cop414
- cop417

Stratégies de séquençage comparées:

1. Recette1 - Aviti uniquement (recette1_aviti)
   → recipes-ShortInsert (2x150bp)
   → Technologie: Aviti

2. Recette2 - Aviti (recette2_aviti)
   → recipes-v3.3.x-ShortInsert-SmallFragmentRecovery (2x150bp)
   → Technologie: Aviti
   → Optimisé pour récupération de fragments courts (aDNA)

3. Recette2 - Illumina (recette2_illumina)
   → recipes-v3.3.x-ShortInsert-SmallFragmentRecovery (2x150bp)
   → Technologie: Illumina
   → Optimisé pour récupération de fragments courts (aDNA)

4. Recette2 - Aviti + Illumina combinés (recette2_aviti_illumina_combined)
   → Fusion des datasets Aviti et Illumina de la Recette2
   → Maximise la profondeur de séquençage

================================================================================
ÉTAPES DU PIPELINE
================================================================================

1. Contrôle qualité (FastQC/MultiQC)
2. Filtrage et trimming (BBDuk) - Élimination des adaptateurs et phiX
3. Déduplication (FastUniq) - Suppression des duplicats PCR
4. Clumpify - Déduplication optique supplémentaire
5. Fastp - Merging des reads et contrôle qualité final
6. Classification taxonomique (Kraken2) - Base de données nt
7. Visualisation (Krona) - Charts interactifs de composition taxonomique
8. Tables d'assignation (MPA format) - Pour analyses downstream

IMPORTANT: Tous les outils sont exécutés depuis un environnement conda unique
           (coprolites-pipeline) pour éviter les erreurs de transition.

================================================================================
LOCALISATION DES RÉSULTATS
================================================================================

Répertoire principal: /home/plstenge/coprolites_comparison

Structure des résultats:

├── 01_raw_data/                        # Données brutes organisées
│   ├── recette1_aviti/                 # Recette1 (Aviti)
│   ├── recette2_aviti/                 # Recette2 (Aviti)
│   ├── recette2_illumina/              # Recette2 (Illumina)
│   └── recette2_aviti_illumina_combined/  # Recette2 (Aviti+Illumina)

├── 02_quality_check/                   # Rapports FastQC/MultiQC
├── 03_bbduk/                           # Reads filtrés
├── 04_fastuniq/                        # Reads dédupliqués
├── 05_clumpify/                        # Déduplication optique
├── 06_fastp/                           # Reads merged et unmerged finaux

├── 07_kraken2/                         # Classification taxonomique
│   ├── recette1_aviti/
│   ├── recette2_aviti/
│   ├── recette2_illumina/
│   └── recette2_aviti_illumina_combined/

├── 08_krona/                           # Visualisations interactives
│   ├── recette1_aviti/
│   ├── recette2_aviti/
│   ├── recette2_illumina/
│   └── recette2_aviti_illumina_combined/

└── 09_mpa_tables/                      # Tables d'assignation
    ├── recette1_aviti/
    ├── recette2_aviti/
    ├── recette2_illumina/
    └── recette2_aviti_illumina_combined/

================================================================================
ANALYSES DE COMPARAISON RECOMMANDÉES
================================================================================

Pour comparer les performances des 4 stratégies:

1. Nombre de reads à chaque étape:
   - Comparez le taux de rétention entre les stratégies
   - Vérifiez le taux de merging (indicateur de qualité aDNA)
   - Comparez Aviti vs Illumina dans Recette2

2. Composition taxonomique:
   - Ouvrez les fichiers Krona pour visualisation interactive
   - Comparez la diversité détectée entre recettes et technologies
   - Identifiez les taxons spécifiques à chaque stratégie

3. Analyse quantitative:
   - Utilisez les fichiers combined_all.tsv dans 09_mpa_tables/
   - Importez dans R pour analyses statistiques et visualisations

4. Évaluation comparative:
   - Recette2 vs Recette1: Amélioration de la récupération de fragments courts?
   - Aviti vs Illumina (Recette2): Différences technologiques?
   - Combined vs séparé: Valeur ajoutée de la fusion?

================================================================================
COMMANDES POUR ANALYSES COMPLÉMENTAIRES EN R
================================================================================

# Chargement et comparaison des tables MPA
library(tidyverse)

# Chargement des 4 stratégies
mpa_r1_aviti <- read_tsv("/home/plstenge/coprolites_comparison/09_mpa_tables/recette1_aviti/combined_all.tsv")
mpa_r2_aviti <- read_tsv("/home/plstenge/coprolites_comparison/09_mpa_tables/recette2_aviti/combined_all.tsv")
mpa_r2_illumina <- read_tsv("/home/plstenge/coprolites_comparison/09_mpa_tables/recette2_illumina/combined_all.tsv")
mpa_r2_combined <- read_tsv("/home/plstenge/coprolites_comparison/09_mpa_tables/recette2_aviti_illumina_combined/combined_all.tsv")

# Comparaison de la richesse taxonomique
richness_comparison <- tibble(
  Strategy = c("Recette1_Aviti", "Recette2_Aviti", "Recette2_Illumina", "Recette2_Combined"),
  N_taxa = c(nrow(mpa_r1_aviti), nrow(mpa_r2_aviti), nrow(mpa_r2_illumina), nrow(mpa_r2_combined))
)

# Visualisation
ggplot(richness_comparison, aes(x = Strategy, y = N_taxa, fill = Strategy)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Richesse taxonomique par stratégie",
       y = "Nombre de taxons détectés") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Analyse de Venn diagram (taxons partagés/uniques)
library(ggvenn)

taxa_lists <- list(
  Recette1_Aviti = mpa_r1_aviti$taxonomy,
  Recette2_Aviti = mpa_r2_aviti$taxonomy,
  Recette2_Illumina = mpa_r2_illumina$taxonomy,
  Recette2_Combined = mpa_r2_combined$taxonomy
)

ggvenn(taxa_lists)

================================================================================
INFORMATIONS TECHNIQUES
================================================================================

Base de données Kraken2: /home/plstenge/k2_core_nt_20250609
KrakenTools: /home/plstenge/coprolites_comparison/07_kraken2/KrakenTools
Environnement conda: coprolites-pipeline

Threads utilisés: 36
Mémoire allouée: 1000G

================================================================================
CONTACT
================================================================================

Pour toute question: pierrelouis.stenger@gmail.com

================================================================================
EOF

# Substitution de la date dans le rapport
sed -i "s/\$(date)/$(date)/" "$REPORT_FILE"

cat "$REPORT_FILE"

echo ""
echo "Rapport généré: $REPORT_FILE"

################################################################################
# ÉTAPE 10: Statistiques finales
################################################################################

echo ""
echo "=== ÉTAPE 10: Génération des statistiques ==="
echo ""

STATS_FILE="${BASE_DIR}/00_scripts/pipeline_statistics.txt"

cat > "$STATS_FILE" << EOF
================================================================================
STATISTIQUES D'ANALYSE - PROJET COPROLITES
================================================================================
Date: $(date)

EOF

for strategy in "${STRATEGIES[@]}"; do
    echo "Statistiques pour ${strategy}..." | tee -a "$STATS_FILE"
    echo "----------------------------------------" | tee -a "$STATS_FILE"
    
    # Compter les fichiers à différentes étapes
    RAW_COUNT=$(ls "${BASE_DIR}/01_raw_data/${strategy}"/*.fastq.gz 2>/dev/null | wc -l)
    BBDUK_COUNT=$(ls "${BASE_DIR}/03_bbduk/${strategy}"/clean_*.fastq.gz 2>/dev/null | wc -l)
    FASTP_MERGED_COUNT=$(ls "${BASE_DIR}/06_fastp/${strategy}"/*_merged.fastq.gz 2>/dev/null | wc -l)
    KRAKEN_COUNT=$(ls "${BASE_DIR}/07_kraken2/${strategy}"/*.report 2>/dev/null | wc -l)
    
    echo "  Fichiers bruts: ${RAW_COUNT}" | tee -a "$STATS_FILE"
    echo "  Fichiers après BBDuk: ${BBDUK_COUNT}" | tee -a "$STATS_FILE"
    echo "  Fichiers merged (Fastp): ${FASTP_MERGED_COUNT}" | tee -a "$STATS_FILE"
    echo "  Rapports Kraken2: ${KRAKEN_COUNT}" | tee -a "$STATS_FILE"
    echo "" | tee -a "$STATS_FILE"
done

cat "$STATS_FILE"

################################################################################
# FIN DU PIPELINE
################################################################################

echo ""
echo "=========================================="
echo "PIPELINE TERMINÉ AVEC SUCCÈS"
echo "Date de fin: $(date)"
echo "=========================================="
echo ""
echo "Résultats disponibles dans: ${BASE_DIR}"
echo "Rapport: ${REPORT_FILE}"
echo "Statistiques: ${STATS_FILE}"
echo ""
echo "Prochaines étapes:"
echo " 1. Consultez les rapports MultiQC dans 02_quality_check/"
echo " 2. Explorez les visualisations Krona dans 08_krona/"
echo " 3. Analysez les tables MPA dans 09_mpa_tables/ avec R"
echo " 4. Comparez les 4 stratégies:"
echo "    - Recette1 (Aviti)"
echo "    - Recette2 (Aviti)"
echo "    - Recette2 (Illumina)"
echo "    - Recette2 (Aviti + Illumina combinés)"
echo ""
echo "Environnement conda utilisé: coprolites-pipeline"
echo "  (Un seul environnement, aucune erreur de transition!)"
echo ""

# Désactivation de l'environnement conda
conda deactivate

# Envoi d'un email de notification (si configuré)
echo "Pipeline Coprolites terminé le $(date). Résultats: ${BASE_DIR}" | \
    mail -s "Pipeline Coprolites - Terminé" pierrelouis.stenger@gmail.com 2>/dev/null || true

exit 0
