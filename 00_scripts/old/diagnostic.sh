#!/bin/bash

#SBATCH --job-name=coprolites_mapdamage
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=500G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/mapdamage.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/mapdamage.out"

set -eo pipefail

BASE_DIR="/home/plstenge/coprolites_comparison"
KRAKENTOOLS_DIR="${BASE_DIR}/08_kraken2/KrakenTools"

################################################################################
# ÉTAPE 8: MapDamage - VERSION MODULES (sans conda)
################################################################################

echo ""
echo "=========================================="
echo "ÉTAPE 8: MapDamage - Analyse des dommages"
echo "=========================================="
echo ""

# ========== CHARGER LES MODULES ==========
echo "Chargement des modules..."

# Essayer les modules disponibles (adapter selon votre cluster)
module load bwa/0.7.17 2>/dev/null || \
    module load bwa 2>/dev/null || \
    { echo "⚠ bwa non trouvé via modules"; }

module load samtools/1.15 2>/dev/null || \
    module load samtools 2>/dev/null || \
    { echo "⚠ samtools non trouvé via modules"; }

module load mapdamage 2>/dev/null || \
    { echo "ℹ mapdamage pas dispo via modules (sera installé via conda)"; }

# Si mapdamage pas dispo via modules, utiliser conda
if ! command -v mapDamage &>/dev/null; then
    echo "Activation conda pour mapdamage..."
    source ~/.bashrc 2>/dev/null
    module load conda/4.12.0 2>/dev/null
    eval "$(conda shell.bash hook)" 2>/dev/null
    conda activate mapdamage_py39 2>/dev/null || \
        { echo "✗ ERREUR: Impossible d'activer mapdamage"; exit 1; }
fi

# Vérification finale
echo ""
echo "Vérification des outils:"
for tool in bwa samtools; do
    if command -v "$tool" &>/dev/null; then
        echo "✓ $tool: $(which $tool)"
    else
        echo "✗ $tool: NON TROUVÉ"
    fi
done

if command -v mapDamage &>/dev/null; then
    echo "✓ mapDamage: $(which mapDamage)"
else
    echo "⚠ mapDamage: non essentiel, continuer..."
fi

# ========== SETUP RÉPERTOIRES ==========

DAMAGE_MERGED_BASE="${BASE_DIR}/12_mapdamage/merged_reads"
DAMAGE_UNMERGED_BASE="${BASE_DIR}/12_mapdamage/unmerged_reads"
mkdir -p "$DAMAGE_MERGED_BASE"
mkdir -p "$DAMAGE_UNMERGED_BASE"

LOGFILE="${BASE_DIR}/00_scripts/mapdamage_$(date +%Y%m%d_%H%M%S).txt"
touch "$LOGFILE"
MAPPING_INFO="${BASE_DIR}/11_summary_tables/mapping_bwa_info_merged.tsv"

echo "$(date): Script MapDamage démarré" | tee -a "$LOGFILE"

# Initialiser le fichier de mapping info
echo -e "Sample\tSpecies\tType\tTotal_Reads\tMapped_Reads\tMapping_Rate(%)" > "$MAPPING_INFO"

# Génomes de référence
declare -A TAXONS=(
    ["Ovis_aries"]="9940:/home/plstenge/genomes/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa"
    ["Capra_hircus"]="9925:/home/plstenge/genomes/Capra_hircus.ARS1.dna.toplevel.fa"
    ["Phragmites_australis"]="35393:/home/plstenge/genomes/Phragmites_australis_GCA_040373225.1.fasta"
)

################################################################################
# FONCTIONS
################################################################################

calculate_mapping_rate() {
    local bam_file="$1"
    local sample_name="$2"
    local species="$3"
    local type="$4"
    
    if [[ -f "$bam_file" ]]; then
        local total=$(samtools view -c "$bam_file" 2>/dev/null || echo 0)
        local mapped=$(samtools view -c -F 4 "$bam_file" 2>/dev/null || echo 0)
        local rate=0
        
        if [[ $total -gt 0 ]]; then
            rate=$(echo "scale=1; $mapped * 100 / $total" | bc 2>/dev/null || echo "0")
        fi
        
        echo -e "${sample_name}\t${species}\t${type}\t${total}\t${mapped}\t${rate}" >> "$MAPPING_INFO"
        echo "  ✓ ${species}: ${mapped}/${total} (${rate}%)" | tee -a "$LOGFILE"
    fi
}

run_mapdamage_safe() {
    local bam_file="$1"
    local ref_fasta="$2"
    local output_dir="$3"
    
    if [[ ! -s "$bam_file" ]]; then
        return 0
    fi
    
    if command -v mapDamage &>/dev/null; then
        mkdir -p "$output_dir"
        mapDamage -i "$bam_file" -r "$ref_fasta" --folder "$output_dir" --no-stats 2>>"$LOGFILE" || true
    fi
}

################################################################################
# INDEXATION DES GÉNOMES
################################################################################

echo ""
echo "→ Préparation des génomes..."

declare -A FIXED_GENOMES

for GROUP in "${!TAXONS[@]}"; do
    ORIGINAL_REF="${TAXONS[$GROUP]#*:}"
    FIXED_REF="/home/plstenge/genomes/$(basename "${ORIGINAL_REF%.*}").fixed.fa"
    
    echo ""
    echo "  ${GROUP}..."
    
    if [[ ! -f "$ORIGINAL_REF" ]]; then
        echo "    ✗ Non trouvé: $ORIGINAL_REF"
        continue
    fi
    
    # Créer lien
    ln -sf "$ORIGINAL_REF" "$FIXED_REF" 2>/dev/null
    FIXED_GENOMES[$GROUP]="$FIXED_REF"
    
    # Indexation BWA avec timeout
    if [[ ! -f "${FIXED_REF}.bwt" ]]; then
        echo "    → Indexation BWA (timeout: 30min)..."
        timeout 1800 bwa index "$FIXED_REF" 2>&1 | \
            grep -v "^\[" | tail -3 >> "$LOGFILE" || {
            if [[ $? -eq 124 ]]; then
                echo "    ⚠ TIMEOUT (génome trop gros)"
            fi
        }
        
        if [[ -f "${FIXED_REF}.bwt" ]]; then
            echo "    ✓ Index BWA créé"
        else
            echo "    ℹ Pas de .bwt (BWA peut être lent sur ce cluster)"
        fi
    else
        echo "    ✓ Index BWA existant"
    fi
    
    # Index samtools (toujours)
    samtools faidx "$FIXED_REF" 2>/dev/null || true
done

################################################################################
# MAPPINGS
################################################################################

shopt -s nullglob

KRAKEN_MERGED="${BASE_DIR}/08_kraken2/merged_reads"
FINAL_MERGED_DIR="${BASE_DIR}/06_fastp/merged_reads"

echo ""
echo "=========================================="
echo "ÉTAPE 8a: Mapping MERGED"
echo "=========================================="

if [[ -d "$KRAKEN_MERGED" ]]; then
    for KRAKEN_FILE in "$KRAKEN_MERGED"/*.kraken; do
        [[ -f "$KRAKEN_FILE" ]] || continue
        KRAKEN_BASE=$(basename "$KRAKEN_FILE" .kraken)
        SAMPLE="${KRAKEN_BASE%_merged}"
        FASTQ="${FINAL_MERGED_DIR}/${SAMPLE}_final_merged.fastq.gz"
        
        [[ ! -f "$FASTQ" || ! -s "$FASTQ" ]] && continue
        
        echo ""
        echo "  ${SAMPLE} (merged)..."
        
        for GROUP in "${!TAXONS[@]}"; do
            REF="${FIXED_GENOMES[$GROUP]}"
            TAX_ID="${TAXONS[$GROUP]%:*}"
            [[ -z "$REF" ]] && continue
            
            OUTDIR="${DAMAGE_MERGED_BASE}/${GROUP}"
            mkdir -p "$OUTDIR"
            
            OUTFQ="${OUTDIR}/${KRAKEN_BASE}_${GROUP}.fastq"
            
            # Extraction
            python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
                -k "$KRAKEN_FILE" -s "$FASTQ" -t "$TAX_ID" \
                -o "$OUTFQ" --fastq-output 2>>"$LOGFILE"
            
            [[ ! -f "$OUTFQ" || ! -s "$OUTFQ" ]] && continue
            
            # Mapping
            bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF" "$OUTFQ" \
                > "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sai" 2>>"$LOGFILE"
            
            bwa samse "$REF" "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sai" "$OUTFQ" 2>>"$LOGFILE" | \
                samtools view -bS - 2>>"$LOGFILE" | \
                samtools sort -o "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" - 2>>"$LOGFILE"
            
            samtools index "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" 2>>"$LOGFILE"
            
            calculate_mapping_rate "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" "$SAMPLE" "$GROUP" "merged"
            
            rm -f "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sai" "$OUTFQ"
            
            run_mapdamage_safe "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" "$REF" \
                "${OUTDIR}/${KRAKEN_BASE}_${GROUP}_mapDamage"
        done
    done
fi

KRAKEN_UNMERGED="${BASE_DIR}/08_kraken2/unmerged_reads"
FINAL_UNMERGED_DIR="${BASE_DIR}/06_fastp/unmerged_reads"

echo ""
echo "=========================================="
echo "ÉTAPE 8b: Mapping UNMERGED"
echo "=========================================="

if [[ -d "$KRAKEN_UNMERGED" ]]; then
    for KRAKEN_FILE in "$KRAKEN_UNMERGED"/*.kraken; do
        [[ -f "$KRAKEN_FILE" ]] || continue
        KRAKEN_BASE=$(basename "$KRAKEN_FILE" .kraken)
        SAMPLE="${KRAKEN_BASE%_unmerged}"
        R1="${FINAL_UNMERGED_DIR}/${SAMPLE}_final_unmerged_R1.fastq.gz"
        R2="${FINAL_UNMERGED_DIR}/${SAMPLE}_final_unmerged_R2.fastq.gz"
        
        [[ ! -f "$R1" || ! -f "$R2" || ! -s "$R1" || ! -s "$R2" ]] && continue
        
        echo ""
        echo "  ${SAMPLE} (unmerged)..."
        
        for GROUP in "${!TAXONS[@]}"; do
            REF="${FIXED_GENOMES[$GROUP]}"
            TAX_ID="${TAXONS[$GROUP]%:*}"
            [[ -z "$REF" ]] && continue
            
            OUTDIR="${DAMAGE_UNMERGED_BASE}/${GROUP}"
            mkdir -p "$OUTDIR"
            
            OUT_R1="${OUTDIR}/${KRAKEN_BASE}_${GROUP}_R1.fastq"
            OUT_R2="${OUTDIR}/${KRAKEN_BASE}_${GROUP}_R2.fastq"
            
            # Extraction
            python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
                -k "$KRAKEN_FILE" -s "$R1" -s2 "$R2" -t "$TAX_ID" \
                -o "$OUT_R1" -o2 "$OUT_R2" --fastq-output 2>>"$LOGFILE"
            
            [[ ! -f "$OUT_R1" || ! -f "$OUT_R2" || ! -s "$OUT_R1" || ! -s "$OUT_R2" ]] && continue
            
            # Mapping
            bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF" "$OUT_R1" \
                > "${OUTDIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" 2>>"$LOGFILE"
            bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF" "$OUT_R2" \
                > "${OUTDIR}/${KRAKEN_BASE}_${GROUP}_R2.sai" 2>>"$LOGFILE"
            
            bwa sampe "$REF" \
                "${OUTDIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" \
                "${OUTDIR}/${KRAKEN_BASE}_${GROUP}_R2.sai" \
                "$OUT_R1" "$OUT_R2" 2>>"$LOGFILE" | \
                samtools view -bS - 2>>"$LOGFILE" | \
                samtools sort -o "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" - 2>>"$LOGFILE"
            
            samtools index "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" 2>>"$LOGFILE"
            
            calculate_mapping_rate "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" "$SAMPLE" "$GROUP" "unmerged"
            
            rm -f "${OUTDIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" \
                  "${OUTDIR}/${KRAKEN_BASE}_${GROUP}_R2.sai" \
                  "$OUT_R1" "$OUT_R2"
            
            run_mapdamage_safe "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" "$REF" \
                "${OUTDIR}/${KRAKEN_BASE}_${GROUP}_mapDamage"
        done
    done
fi

shopt -u nullglob

echo ""
echo "=========================================="
echo "FIN - $(date)"
echo "=========================================="
echo ""
echo "Résultats: $MAPPING_INFO"
echo "Log complet: $LOGFILE"
echo ""

exit 0
