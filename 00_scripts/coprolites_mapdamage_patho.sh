#!/bin/bash

#SBATCH --job-name=coprolites_mapdamage
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=500G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/mapdamage.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/mapdamage.out"

echo ""
echo "=========================================="
echo "ÉTAPE 8: MapDamage - Analyse des dommages (UNMERGED READS)"
echo "=========================================="
echo ""

# ========== INITIALISATION CONDA POUR SLURM ==========
echo "Initialisation de conda..."
module load conda/4.12.0
source ~/.bashrc
conda activate mapdamage_py39
echo "✓ Environnement mapdamage_py39 activé"

# ========== CONFIGURATION GLOBALE ==========
BASE_DIR="/home/plstenge/coprolites_comparison"
KRAKENTOOLS_DIR="${BASE_DIR}/08_kraken2/KrakenTools"
THREADS=36

echo ""
echo "Configuration:"
echo "  BASE_DIR: $BASE_DIR"
echo "  KRAKENTOOLS_DIR: $KRAKENTOOLS_DIR"

# ========== VÉRIFICATION DES OUTILS ==========
echo ""
echo "Vérification des outils disponibles:"
for tool in bwa samtools mapDamage python3; do
    if command -v "$tool" &>/dev/null; then
        echo "  ✓ $tool: $(which $tool)"
    else
        echo "  ✗ $tool: NON TROUVÉ"
        [[ "$tool" == "mapDamage" ]] && { echo "ERREUR CRITIQUE: mapDamage requis"; exit 1; }
    fi
done

echo ""
echo "Versions:"
bwa 2>&1 | head -1
samtools --version | head -1
mapDamage --version 2>&1 | head -1 || echo "mapDamage: version inconnue"

# ========== VÉRIFICATION DES RÉPERTOIRES SOURCE ==========
echo ""
echo "Vérification des répertoires source:"
KRAKEN_UNMERGED="${BASE_DIR}/08_kraken2/unmerged_reads"
FINAL_UNMERGED_DIR="${BASE_DIR}/06_fastp/unmerged_reads"

echo "  KRAKEN_UNMERGED: $KRAKEN_UNMERGED"
if [[ -d "$KRAKEN_UNMERGED" ]]; then
    N_KRAKEN=$(ls -1 "$KRAKEN_UNMERGED"/*.kraken 2>/dev/null | wc -l)
    echo "  ✓ Existe ($N_KRAKEN fichiers .kraken)"
else
    echo "  ✗ MANQUANT"
    exit 1
fi

echo "  FINAL_UNMERGED_DIR: $FINAL_UNMERGED_DIR"
if [[ -d "$FINAL_UNMERGED_DIR" ]]; then
    N_FQ=$(ls -1 "$FINAL_UNMERGED_DIR"/*.fastq.gz 2>/dev/null | wc -l)
    echo "  ✓ Existe ($N_FQ fichiers fastq.gz)"
else
    echo "  ✗ MANQUANT"
    exit 1
fi

# ========== SETUP RÉPERTOIRES ==========
DAMAGE_UNMERGED_BASE="${BASE_DIR}/12_mapdamage/unmerged_reads"
mkdir -p "$DAMAGE_UNMERGED_BASE"

LOGFILE="${BASE_DIR}/00_scripts/mapdamage_$(date +%Y%m%d_%H%M%S).txt"
touch "$LOGFILE"

MAPPING_INFO="${BASE_DIR}/11_summary_tables/mapping_bwa_info_unmerged.tsv"
mkdir -p "${BASE_DIR}/11_summary_tables"

echo "$(date): Script MapDamage UNMERGED démarré" | tee -a "$LOGFILE"

# Initialiser le fichier de mapping info
echo -e "Sample\tSpecies\tType\tTotal_Reads\tMapped_Reads\tMapping_Rate(%)" > "$MAPPING_INFO"

# ========== GÉNOMES DE RÉFÉRENCE - CONFIGURATION ==========
# Chaque entrée : TAXONS["nom_espece"]="taxid:/chemin/vers/genome.fasta"
declare -A TAXONS=(
    ["Rhodococcus_globerulus"]="33008:/home/plstenge/genomes/Rhodococcus_globerulus/GCF_019334125.1_ASM1933412v1_genomic.fna"
    ["Corynebacterium_occultum"]="2675219:/home/plstenge/genomes/Corynebacterium_occultum/GCF_009734425.1_ASM973442v1_genomic.fna"
    ["Williamsia_herbipolensis"]="1603258:/home/plstenge/genomes/Williamsia_herbipolensis/GCF_036153515.1_ASM3615351v1_genomic.fna"
    ["Kroppenstedtia_eburnea"]="714067:/home/plstenge/genomes/Kroppenstedtia_eburnea/GCF_013282215.1_ASM1328221v1_genomic.fna"
)
	
	TAXXXIOOOOONN IIIIIIIDDDD

echo ""
echo "Génomes de référence configurés: ${#TAXONS[@]}"
echo ""

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
        echo "    Stats: ${mapped}/${total} mappés (${rate}%)" | tee -a "$LOGFILE"
    fi
}

# FONCTION CORRIGÉE : --no-stats ajouté pour TOUS les modes (paired ET single)
run_mapdamage_with_length_dist() {
    local bam_file="$1"
    local ref_fasta="$2"
    local output_dir="$3"
    local read_type="$4"
    
    if [[ ! -s "$bam_file" ]]; then
        echo "    ⚠ Fichier BAM vide ou inexistant"
        return 0
    fi
    
    mkdir -p "$output_dir"
    echo -n "    → MapDamage en cours..."
    
    # CORRECTION PRINCIPALE : Désactiver les statistiques bayésiennes pour éviter l'erreur Rcpp/GCC
    # Cela contourne le problème de compilation tout en conservant les graphiques de dommages
    if mapDamage -i "$bam_file" -r "$ref_fasta" --folder "$output_dir" --no-stats 2>>"$LOGFILE"; then
        echo " ✓ OK (graphiques de dommages générés, stats bayésiennes désactivées)"
    else
        echo " ⚠ (échoué - non bloquant)"
    fi
}

check_bwa_index_complete() {
    local ref="$1"
    if [[ -f "${ref}.bwt" ]] && \
       [[ -f "${ref}.amb" ]] && \
       [[ -f "${ref}.ann" ]] && \
       [[ -f "${ref}.pac" ]] && \
       [[ -f "${ref}.sa" ]]; then
        return 0  # Complet
    fi
    return 1  # Incomplet
}

################################################################################
# PHASE 1 : INDEXATION PRÉALABLE DES GÉNOMES
################################################################################

echo ""
echo "=========================================="
echo "PHASE 1: Préparation et indexation des génomes de référence"
echo "=========================================="
echo ""

declare -A VALID_GENOMES
INDEXATION_ERRORS=0

for GROUP in "${!TAXONS[@]}"; do
    ORIGINAL_REF="${TAXONS[$GROUP]#*:}"
    echo -n "→ ${GROUP}: "
    
    # Vérifier que le FASTA existe
    if [[ ! -f "$ORIGINAL_REF" ]]; then
        echo "✗ Génome FASTA non trouvé: $ORIGINAL_REF"
        ((INDEXATION_ERRORS++))
        continue
    fi
    
    GENOME_SIZE=$(du -h "$ORIGINAL_REF" | cut -f1)
    echo -n "($GENOME_SIZE) ... "
    
    # ====== INDEX BWA ======
    if check_bwa_index_complete "$ORIGINAL_REF"; then
        echo -n "[BWA:✓] "
    else
        echo -n "[BWA:⏳ indexation] "
        if timeout 3600 bwa index "$ORIGINAL_REF" 2>>"$LOGFILE"; then
            if check_bwa_index_complete "$ORIGINAL_REF"; then
                echo -n "[✓] "
            else
                echo "✗ [ERREUR: Index BWA incomplet après indexation]"
                ((INDEXATION_ERRORS++))
                continue
            fi
        else
            echo "✗ [ERREUR: BWA index timeout ou erreur]"
            ((INDEXATION_ERRORS++))
            continue
        fi
    fi
    
    # ====== INDEX SAMTOOLS ======
    if [[ ! -f "${ORIGINAL_REF}.fai" ]]; then
        echo -n "[SAM:⏳] "
        if samtools faidx "$ORIGINAL_REF" 2>>"$LOGFILE"; then
            echo "[✓]"
        else
            echo "✗ [ERREUR: samtools faidx échoué]"
            ((INDEXATION_ERRORS++))
            continue
        fi
    else
        echo "[SAM:✓]"
    fi
    
    # Stocker le génome comme valide
    VALID_GENOMES[$GROUP]="$ORIGINAL_REF"
done

echo ""
if [[ $INDEXATION_ERRORS -gt 0 ]]; then
    echo "⚠ WARNING: $INDEXATION_ERRORS génomes ont échoué l'indexation"
fi
echo "✓ Génomes prêts pour le mapping: ${#VALID_GENOMES[@]}"

if [[ ${#VALID_GENOMES[@]} -eq 0 ]]; then
    echo "✗ ERREUR CRITIQUE: Aucun génome valide disponible!"
    exit 1
fi

echo ""
echo "Génomes valides pour la suite :"
for sp in "${!VALID_GENOMES[@]}"; do
    echo "  ✓ $sp"
done
echo ""

################################################################################
# PHASE 2 : MAPPINGS - UNMERGED READS (PAIRED-END)
################################################################################

shopt -s nullglob

echo ""
echo "=========================================="
echo "PHASE 2: Analyse UNMERGED READS (Paired-End) - MapDamage"
echo "=========================================="
echo ""

if [[ -d "$KRAKEN_UNMERGED" ]] && ls "$KRAKEN_UNMERGED"/*.kraken >/dev/null 2>&1; then
    for KRAKEN_FILE in "$KRAKEN_UNMERGED"/*.kraken; do
        KRAKEN_BASE=$(basename "$KRAKEN_FILE" .kraken)
        SAMPLE="${KRAKEN_BASE%_unmerged}"
        
        R1="${FINAL_UNMERGED_DIR}/${SAMPLE}_final_unmerged_R1.fastq.gz"
        R2="${FINAL_UNMERGED_DIR}/${SAMPLE}_final_unmerged_R2.fastq.gz"
        
        # Vérifier que les fichiers existent et ne sont pas vides
        if [[ ! -f "$R1" || ! -f "$R2" || ! -s "$R1" || ! -s "$R2" ]]; then
            echo ""
            echo "⚠ Fichiers unmerged manquants ou vides pour $SAMPLE"
            continue
        fi
        
        echo ""
        echo "┌─────────────────────────────────────────"
        echo "│ ${SAMPLE} (unmerged - $(ls -lh $R1 | awk '{print $5}') + $(ls -lh $R2 | awk '{print $5}'))"
        echo "└─────────────────────────────────────────"
        
        for GROUP in "${!VALID_GENOMES[@]}"; do
            REF="${VALID_GENOMES[$GROUP]}"
            TAX_ID="${TAXONS[$GROUP]%:*}"
            
            if [[ -z "$REF" || ! -f "$REF" ]]; then
                echo "    ⚠ Génome non disponible pour ${GROUP}"
                continue
            fi
            
            OUTDIR="${DAMAGE_UNMERGED_BASE}/${GROUP}/paired_end"
            mkdir -p "$OUTDIR"
            
            OUT_R1="${OUTDIR}/${KRAKEN_BASE}_${GROUP}_R1.fastq"
            OUT_R2="${OUTDIR}/${KRAKEN_BASE}_${GROUP}_R2.fastq"
            
            echo ""
            echo "  → ${GROUP} (TaxID: ${TAX_ID})"
            echo "    Extraction reads paired-end..."
            
            # Extraire les reads pour ce taxon
            python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
                -k "$KRAKEN_FILE" \
                -s "$R1" -s2 "$R2" \
                -t "$TAX_ID" \
                -o "$OUT_R1" -o2 "$OUT_R2" \
                --fastq-output 2>>"$LOGFILE"
            
            # Vérifier qu'on a des résultats
            if [[ ! -f "$OUT_R1" || ! -f "$OUT_R2" || ! -s "$OUT_R1" || ! -s "$OUT_R2" ]]; then
                echo "    ⚠ Aucun read extrait pour ${GROUP}"
                rm -f "$OUT_R1" "$OUT_R2" 2>/dev/null
                continue
            fi
            
            READ_COUNT=$(grep -c "^@" "$OUT_R1" 2>/dev/null || echo 0)
            echo "    ✓ ${READ_COUNT} paires extraites"
            
            # Mapping BWA paired-end
            echo "    Mapping BWA paired-end..."
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
            
            if [[ ! -f "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" ]]; then
                echo "    ✗ Erreur lors du mapping"
                continue
            fi
            
            samtools index "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" 2>>"$LOGFILE"
            echo "    ✓ Mapping terminé"
            
            # Calculer taux de mapping
            calculate_mapping_rate "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" "$SAMPLE" "$GROUP" "unmerged_PE"
            
            # MapDamage (paired-end : pas de distribution de longueurs)
            run_mapdamage_with_length_dist "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" "$REF" \
                "${OUTDIR}/${KRAKEN_BASE}_${GROUP}_mapDamage" "paired"
            
            # Nettoyage (garder les BAM indexés pour futures analyses)
            rm -f "${OUTDIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" \
                  "${OUTDIR}/${KRAKEN_BASE}_${GROUP}_R2.sai" \
                  "$OUT_R1" "$OUT_R2"
        done
    done
else
    echo "✗ ERREUR: Aucun fichier .kraken trouvé dans ${KRAKEN_UNMERGED}"
    exit 1
fi

################################################################################
# PHASE 3 : ANALYSE SINGLE-END R1 UNIQUEMENT (POUR DISTRIBUTION DE LONGUEURS)
################################################################################

echo ""
echo "=========================================="
echo "PHASE 3: Analyse UNMERGED R1 (Single-End) - pour distributions de longueurs"
echo "=========================================="
echo ""

if [[ -d "$KRAKEN_UNMERGED" ]] && ls "$KRAKEN_UNMERGED"/*.kraken >/dev/null 2>&1; then
    for KRAKEN_FILE in "$KRAKEN_UNMERGED"/*.kraken; do
        KRAKEN_BASE=$(basename "$KRAKEN_FILE" .kraken)
        SAMPLE="${KRAKEN_BASE%_unmerged}"
        
        R1="${FINAL_UNMERGED_DIR}/${SAMPLE}_final_unmerged_R1.fastq.gz"
        
        # Vérifier que le fichier existe et n'est pas vide
        if [[ ! -f "$R1" || ! -s "$R1" ]]; then
            continue
        fi
        
        echo ""
        echo "┌─────────────────────────────────────────"
        echo "│ ${SAMPLE} (R1 only - $(ls -lh $R1 | awk '{print $5}'))"
        echo "└─────────────────────────────────────────"
        
        for GROUP in "${!VALID_GENOMES[@]}"; do
            REF="${VALID_GENOMES[$GROUP]}"
            TAX_ID="${TAXONS[$GROUP]%:*}"
            
            if [[ -z "$REF" || ! -f "$REF" ]]; then
                continue
            fi
            
            OUTDIR="${DAMAGE_UNMERGED_BASE}/${GROUP}/single_end_R1"
            mkdir -p "$OUTDIR"
            
            OUT_R1="${OUTDIR}/${KRAKEN_BASE}_${GROUP}_R1.fastq"
            
            echo ""
            echo "  → ${GROUP} (Single-end R1)"
            echo "    Extraction R1 uniquement..."
            
            # Extraire R1 uniquement
            python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
                -k "$KRAKEN_FILE" \
                -s "$R1" \
                -t "$TAX_ID" \
                -o "$OUT_R1" \
                --fastq-output 2>>"$LOGFILE"
            
            if [[ ! -f "$OUT_R1" || ! -s "$OUT_R1" ]]; then
                echo "    ⚠ Aucun read extrait pour ${GROUP}"
                rm -f "$OUT_R1" 2>/dev/null
                continue
            fi
            
            READ_COUNT=$(grep -c "^@" "$OUT_R1" 2>/dev/null || echo 0)
            echo "    ✓ ${READ_COUNT} reads extraits"
            
            # Mapping BWA single-end
            echo "    Mapping BWA single-end..."
            bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF" "$OUT_R1" \
                > "${OUTDIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" 2>>"$LOGFILE"
            
            bwa samse "$REF" \
                "${OUTDIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" \
                "$OUT_R1" 2>>"$LOGFILE" | \
            samtools view -bS - 2>>"$LOGFILE" | \
            samtools sort -o "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" - 2>>"$LOGFILE"
            
            if [[ ! -f "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" ]]; then
                echo "    ✗ Erreur lors du mapping"
                continue
            fi
            
            samtools index "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" 2>>"$LOGFILE"
            echo "    ✓ Mapping terminé"
            
            # Calculer taux de mapping
            calculate_mapping_rate "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" "$SAMPLE" "$GROUP" "unmerged_SE_R1"
            
            # MapDamage (single-end : normalement AVEC distribution de longueurs, mais stats désactivées)
            run_mapdamage_with_length_dist "${OUTDIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam}" "$REF" \
                "${OUTDIR}/${KRAKEN_BASE}_${GROUP}_mapDamage" "single"
            
            # Nettoyage
            rm -f "${OUTDIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" "$OUT_R1"
        done
    done
fi

shopt -u nullglob

################################################################################
# RÉSUMÉ ET FIN
################################################################################

echo ""
echo "=========================================="
echo "PIPELINE MAPDAMAGE TERMINÉ"
echo "Date: $(date)"
echo "=========================================="
echo ""

echo "📊 Statistiques de mapping: $MAPPING_INFO"
echo "📝 Log complet: $LOGFILE"
echo "📁 Fichiers BAM + MapDamage:"
echo "   - Paired-end: $DAMAGE_UNMERGED_BASE/*/paired_end/"
echo "   - Single-end R1: $DAMAGE_UNMERGED_BASE/*/single_end_R1/"
echo ""

# Afficher un aperçu des résultats
if [[ -f "$MAPPING_INFO" ]]; then
    echo "Aperçu des résultats:"
    echo ""
    column -t "$MAPPING_INFO"
    echo ""
fi

echo "✓ Script terminé avec succès"
echo ""

conda deactivate

# Envoyer notification
echo "Pipeline MapDamage terminé le $(date +%Y-%m-%d\ %H:%M:%S). Résultats: $DAMAGE_UNMERGED_BASE" | \
    mail -s "Pipeline Coprolites MapDamage - Terminé" pierrelouis.stenger@gmail.com 2>/dev/null || true

exit 0
