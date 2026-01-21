#!/bin/bash

################################################################################
# ÉTAPES 9 & 10: PHYLOGÉNIE ADN ANCIEN - MOUTON (VERSION COMPLÈTE FIXÉE)
# Récupération des reads mappés → Consensus mtDNA → Phylogénie (RAxML)
# 
# VERSION: Utilise consensus poolé + références (même incomplètes)
# AUTEUR: Bioinformatique Paléo
# DATE: 2026-01-20
################################################################################

#SBATCH --job-name=coprolites_phylogeny_fixed
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=200G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/phylogeny_fixed.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/phylogeny_fixed.out"

echo ""
echo "=========================================="
echo "PHYLOGÉNIE ADN ANCIEN MOUTON - VERSION FIXÉE"
echo "=========================================="
echo ""
echo "Timestamp: $(date)"
echo ""

# ========== INITIALISATION CONDA ==========
echo "Initialisation des environnements conda..."

module load conda/4.12.0
source ~/.bashrc

# Créer environnement s'il n'existe pas
if ! conda env list | grep -q "phylogeny_env"; then
    echo "Création environnement phylogeny_env..."
    conda create -y -n phylogeny_env \
        -c bioconda -c conda-forge \
        bcftools samtools mafft raxml \
        python=3.10 biopython numpy pandas
fi

conda activate phylogeny_env
echo "✓ Environnement phylogeny_env activé"
echo ""

# ========== CONFIGURATION GLOBALE ==========
BASE_DIR="/home/plstenge/coprolites_comparison"
BAM_EXTRACTION_DIR="${BASE_DIR}/13_phylogeny/01_bam_extraction"
CONSENSUS_DIR="${BASE_DIR}/13_phylogeny/02_consensus"
ALIGN_DIR="${BASE_DIR}/13_phylogeny/03_alignment"
PHYLOGENY_DIR="${BASE_DIR}/13_phylogeny/04_phylogeny"
VISUALIZATION_DIR="${BASE_DIR}/13_phylogeny/05_visualization"
REFERENCES_DIR="${BASE_DIR}/13_phylogeny/00_references"
LOG_DIR="${BASE_DIR}/13_phylogeny/logs"

THREADS=8
MOUFLON_REF="/home/plstenge/genomes/Ovis_aries_mitochondrion/Ovis_aries_mitochondrion.fa"

echo "Configuration:"
echo "  BASE_DIR: $BASE_DIR"
echo "  BAM_EXTRACTION: $BAM_EXTRACTION_DIR"
echo "  CONSENSUS: $CONSENSUS_DIR"
echo "  ALIGNMENT: $ALIGN_DIR"
echo "  PHYLOGENY: $PHYLOGENY_DIR"
echo "  THREADS: $THREADS"
echo ""

# ========== CRÉATION RÉPERTOIRES ==========
mkdir -p "$BAM_EXTRACTION_DIR" "$CONSENSUS_DIR" "$ALIGN_DIR" "$PHYLOGENY_DIR" "$VISUALIZATION_DIR" "$REFERENCES_DIR" "$LOG_DIR"

LOGFILE="${LOG_DIR}/phylogeny_$(date +%Y%m%d_%H%M%S).log"
touch "$LOGFILE"

# ========== VÉRIFICATION DES OUTILS ==========
echo "=========================================="
echo "Vérification des outils"
echo "=========================================="
echo ""

for tool in samtools bcftools mafft raxmlHPC-PTHREADS python3; do
    if command -v "$tool" &>/dev/null; then
        echo "  ✓ $tool"
    else
        echo "  ✗ $tool MANQUANT"
        if [[ "$tool" == "mafft" || "$tool" == "raxmlHPC-PTHREADS" ]]; then
            exit 1
        fi
    fi
done
echo ""

################################################################################
# ÉTAPE 1: CRÉER LE CONSENSUS POOLÉ
################################################################################

echo "=========================================="
echo "ÉTAPE 1: Création consensus poolé"
echo "=========================================="
echo ""

# Chercher tous les BAM mouton individuels
SHEEP_BAMS=()
for bam_file in "${BAM_EXTRACTION_DIR}"/*Ovis_aries*.bam; do
    if [[ -f "$bam_file" ]]; then
        SHEEP_BAMS+=("$bam_file")
    fi
done

if [[ ${#SHEEP_BAMS[@]} -eq 0 ]]; then
    echo "✗ ERREUR: Aucun BAM mouton trouvé dans $BAM_EXTRACTION_DIR"
    exit 1
fi

echo "Nombre de BAM mouton trouvés: ${#SHEEP_BAMS[@]}"
for bam in "${SHEEP_BAMS[@]}"; do
    SAMPLE=$(basename "$bam" .bam)
    READS=$(samtools view -c "$bam" 2>/dev/null || echo "0")
    echo "  → $SAMPLE: $READS reads"
done
echo ""

# Merger tous les BAM
POOLED_BAM="${BAM_EXTRACTION_DIR}/POOLED_all_sheep_merged.bam"

if [[ -f "$POOLED_BAM" ]]; then
    echo "Réutilisation du BAM poolé existant: $POOLED_BAM"
else
    echo "Fusion des BAM en cours..."
    samtools merge -@ $THREADS -f "$POOLED_BAM" "${SHEEP_BAMS[@]}"
    samtools index "$POOLED_BAM"
    echo "✓ BAM poolé créé et indexé"
fi

# Statistiques
TOTAL_READS=$(samtools view -c "$POOLED_BAM")
MAPPED_READS=$(samtools view -c -F 4 "$POOLED_BAM")
COV_MEAN=$(samtools depth "$POOLED_BAM" 2>/dev/null | \
    awk '{sum+=$3; count++} END {if (count>0) printf "%.1f", sum/count; else print "0"}')

echo ""
echo "Statistiques BAM poolé:"
echo "  Total reads: $TOTAL_READS"
echo "  Mapped reads: $MAPPED_READS"
echo "  Couverture moyenne: ${COV_MEAN}X"
echo ""

################################################################################
# ÉTAPE 2: GÉNÉRER CONSENSUS DEPUIS LE BAM POOLÉ
################################################################################

echo "=========================================="
echo "ÉTAPE 2: Génération consensus"
echo "=========================================="
echo ""

CONSENSUS_RAW="${CONSENSUS_DIR}/pooled_consensus_raw.txt"
CONSENSUS_FA="${CONSENSUS_DIR}/POOLED_consensus_final.fa"

if [[ -f "$CONSENSUS_FA" ]]; then
    echo "Réutilisation du consensus existant: $CONSENSUS_FA"
else
    echo "Génération du consensus depuis BAM..."
    
    # Générer consensus en appelant la majorité
    samtools mpileup -f "$MOUFLON_REF" "$POOLED_BAM" 2>/dev/null | \
        awk '{
            depth=$4
            bases=$5
            if (depth>=1) {
                # Nettoyer bases
                gsub(/[^ACGTacgt]/,"",bases)
                
                # Compter bases
                a=gsub(/[aA]/,"",bases)
                c=gsub(/[cC]/,"",bases)  
                g=gsub(/[gG]/,"",bases)
                t=gsub(/[tT]/,"",bases)
                
                # Appeler majorité
                if (a>=c && a>=g && a>=t) printf "A"
                else if (c>=a && c>=g && c>=t) printf "C"
                else if (g>=a && g>=c && g>=t) printf "G"
                else printf "T"
            } else {
                printf "N"
            }
        }' > "$CONSENSUS_RAW"
    
    echo "✓ Consensus brut généré"
    
    # Créer FASTA formaté
    cat > "$CONSENSUS_FA" << EOF
>POOLED_Ancient_Sheep
EOF
    
    cat "$CONSENSUS_RAW" | sed 's/\(.\{80\}\)/\1\n/g' >> "$CONSENSUS_FA"
    
    echo "✓ FASTA généré"
fi

# Statistiques consensus
CONS_LENGTH=$(awk '/^>/ {next} {total+=length($0)} END {print total}' "$CONSENSUS_FA")
N_COUNT=$(grep -v "^>" "$CONSENSUS_FA" | grep -o 'N' | wc -l)
PERCENT_N=$(echo "scale=2; $N_COUNT * 100 / $CONS_LENGTH" | bc 2>/dev/null || echo "0")

echo ""
echo "Statistiques consensus:"
echo "  Longueur: $CONS_LENGTH bp"
echo "  Bases N: $N_COUNT ($PERCENT_N%)"
echo "  Déterminées: $((CONS_LENGTH - N_COUNT)) bp"
echo ""

# Sauvegarder dans all_samples
cat "$CONSENSUS_FA" >> "${CONSENSUS_DIR}/all_samples_consensus_mtDNA.fasta"

################################################################################
# ÉTAPE 3: CRÉER L'ALIGNEMENT MULTIPLE
################################################################################

echo "=========================================="
echo "ÉTAPE 3: Alignement multiple MAFFT"
echo "=========================================="
echo ""

COMBINED_FA="${ALIGN_DIR}/combined_ancient_modern_mtDNA.fasta"
ALIGNED_FA="${ALIGN_DIR}/aligned_mtDNA.fasta"

if [[ ! -f "$COMBINED_FA" ]]; then
    echo "Création fichier combiné..."
    
    # Combiner consensus + références
    if [[ -f "${REFERENCES_DIR}/sheep_reference_breeds_mtDNA.fasta" ]]; then
        cat "$CONSENSUS_FA" "${REFERENCES_DIR}/sheep_reference_breeds_mtDNA.fasta" > "$COMBINED_FA"
        echo "✓ Combiné: $(grep -c '^>' "$COMBINED_FA") séquences"
    else
        echo "⚠ Pas de références, utilisation du consensus seul"
        cp "$CONSENSUS_FA" "$COMBINED_FA"
    fi
    
    # Statistiques
    echo ""
    echo "Composition du fichier combiné:"
    echo "  Nombre de séquences: $(grep -c '^>' "$COMBINED_FA")"
    echo ""
    echo "Distribution des longueurs:"
    grep -v "^>" "$COMBINED_FA" | awk '{print length}' | sort -rn | uniq -c | head -10
    echo ""
fi

# Aligner avec MAFFT
if [[ ! -f "$ALIGNED_FA" ]]; then
    echo "Alignement MAFFT (--auto)..."
    
    mafft --auto --thread $THREADS "$COMBINED_FA" > "$ALIGNED_FA"
    
    echo "✓ Alignement généré"
    
    # Statistiques alignement
    ALIGN_LEN=$(head -n 2 "$ALIGNED_FA" | tail -n 1 | wc -c)
    ALIGN_SEQS=$(grep -c "^>" "$ALIGNED_FA")
    
    echo ""
    echo "Statistiques alignement:"
    echo "  Séquences: $ALIGN_SEQS"
    echo "  Longueur: $ALIGN_LEN bp"
else
    echo "Réutilisation de l'alignement existant: $ALIGNED_FA"
fi

echo ""

################################################################################
# ÉTAPE 4: CONSTRUIRE L'ARBRE PHYLOGÉNÉTIQUE AVEC RAxML
################################################################################

echo "=========================================="
echo "ÉTAPE 4: Phylogénie RAxML"
echo "=========================================="
echo ""

cd "$PHYLOGENY_DIR"

# Nettoyer anciens fichiers RAxML (optionnel)
# rm -f RAxML_* 2>/dev/null

echo "Construction arbre phylogénétique..."
echo "  Modèle: GTRGAMMA"
echo "  Bootstraps: 100"
echo "  Threads: $THREADS"
echo ""

RUN_ID="sheep_mtDNA_ml_$(date +%Y%m%d_%H%M%S)"

# Lancer RAxML
raxmlHPC-PTHREADS \
    -f a \
    -x 12345 \
    -p 12345 \
    -# 100 \
    -m GTRGAMMA \
    -s "$ALIGNED_FA" \
    -n "$RUN_ID" \
    -T $THREADS \
    2>&1 | tee "${LOG_DIR}/raxml_${RUN_ID}.log"

echo ""
echo "✓ RAxML terminé"
echo ""

# Vérifier résultats
BEST_TREE="RAxML_bipartitions.${RUN_ID}"
if [[ -f "$BEST_TREE" ]]; then
    echo "Fichiers générés:"
    ls -lh RAxML_* | awk '{print "  " $9 " (" $5 ")"}'
    echo ""
    echo "✓ Arbre phylogénétique: $BEST_TREE"
else
    echo "✗ ERREUR: Arbre RAxML non généré"
    exit 1
fi

################################################################################
# ÉTAPE 5: PRÉPARATION VISUALISATION
################################################################################

echo ""
echo "=========================================="
echo "ÉTAPE 5: Préparation visualisation"
echo "=========================================="
echo ""

mkdir -p "${VISUALIZATION_DIR}/figtree"

# Copier l'arbre pour visualisation
cp "$BEST_TREE" "${VISUALIZATION_DIR}/figtree/tree.nwk"

echo "Pour visualiser l'arbre:"
echo ""
echo "  1. Installer FigTree (si pas déjà installé):"
echo "     brew install figtree  # macOS"
echo "     apt-get install figtree  # Linux"
echo ""
echo "  2. Ouvrir l'arbre:"
echo "     figtree ${VISUALIZATION_DIR}/figtree/tree.nwk"
echo ""
echo "  3. Annotations recommandées dans FigTree:"
echo "     - Node Labels → Afficher bootstrap (label)"
echo "     - Tip Labels → Gras, couleurs par groupe"
echo "     - Branch Width → Par support (proportional)"
echo "     - Scale Bar → Ajouter pour montrer distances"
echo ""

################################################################################
# RÉSUMÉ FINAL
################################################################################

echo "=========================================="
echo "✓ PIPELINE PHYLOGÉNIE TERMINÉ"
echo "=========================================="
echo ""

echo "📊 RÉSUMÉS:"
echo "  Consensus: $CONSENSUS_FA"
echo "  Alignement: $ALIGNED_FA"
echo "  Arbre RAxML: ${PHYLOGENY_DIR}/${BEST_TREE}"
echo ""

echo "📁 FICHIERS CLÉS:"
echo "  - $BEST_TREE → Arbre avec support bootstrap"
echo "  - ${ALIGN_DIR}/aligned_mtDNA.fasta → Alignement multiple"
echo "  - ${LOG_DIR}/ → Tous les logs détaillés"
echo ""

echo "📈 ÉTAPES SUIVANTES:"
echo ""
echo "  1. VISUALISATION:"
echo "     figtree ${PHYLOGENY_DIR}/${BEST_TREE}"
echo ""
echo "  2. ANALYSER LES RÉSULTATS:"
echo "     - Identifier le position du POOLED_Ancient_Sheep dans l'arbre"
echo "     - Vérifier la profondeur des branches (bootstrap ≥70 pour confiance)"
echo "     - Comparer avec races de référence"
echo ""
echo "  3. INTERPRÉTATION BIOLOGIQUE:"
echo "     - Quelle race/groupe est le plus proche?"
echo "     - Confiance bootstrap?"
echo "     - Implications pour l'archéologie?"
echo ""

echo "📝 LOG COMPLET: $LOGFILE"
echo ""
echo "Timestamp: $(date)"
echo "=========================================="
echo ""

conda deactivate

exit 0
