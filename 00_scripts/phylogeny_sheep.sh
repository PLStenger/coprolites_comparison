#!/bin/bash
#SBATCH --job-name=coprolites_phylogeny
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=200G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/phylogeny.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/phylogeny.out"

################################################################################
# ÉTAPE 9 & 10: PHYLOGÉNIE ADN ANCIEN - MOUTON
# Récupération des reads mappés → Consensus mtDNA → Phylogénie (RAxML + BEAST)
# 
# Dépend de: ÉTAPE 8 (MapDamage, fichiers BAM mappés)
# Génère: Alignements FASTA → Arbres phylogénétiques
################################################################################

echo ""
echo "=========================================="
echo "ÉTAPES 9 & 10: PHYLOGÉNIE ADN ANCIEN MOUTON"
echo "De l'ADN fragmenté à l'arbre phylogénétique"
echo "=========================================="
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
        bcftools samtools mafft raxml figtree \
        python=3.10 biopython numpy pandas
fi

conda activate phylogeny_env
echo "✓ Environnement phylogeny_env activé"

# ========== CONFIGURATION GLOBALE ==========
BASE_DIR="/home/plstenge/coprolites_comparison"
DAMAGE_DIR="${BASE_DIR}/12_mapdamage/unmerged_reads"
PHYLOGENY_BASE="${BASE_DIR}/13_phylogeny"
THREADS=36

echo ""
echo "Configuration:"
echo "  BASE_DIR: $BASE_DIR"
echo "  DAMAGE_DIR: $DAMAGE_DIR"
echo "  PHYLOGENY_BASE: $PHYLOGENY_BASE"
echo "  THREADS: $THREADS"

# ========== VÉRIFICATION DES OUTILS ==========
echo ""
echo "Vérification des outils:"

for tool in samtools bcftools mafft raxml python3; do
    if command -v "$tool" &>/dev/null; then
        echo "  ✓ $tool"
    else
        echo "  ✗ $tool MANQUANT"
        [[ "$tool" != "raxml" ]] && exit 1
    fi
done

# ========== CRÉATION RÉPERTOIRES ==========
mkdir -p "${PHYLOGENY_BASE}"/{01_bam_extraction,02_consensus,03_alignment,04_phylogeny,05_visualization,logs}

LOGFILE="${PHYLOGENY_BASE}/logs/phylogeny_$(date +%Y%m%d_%H%M%S).log"
touch "$LOGFILE"

echo "$(date): Pipeline phylogénie démarré" | tee -a "$LOGFILE"

# ========== GÉNOMES DE RÉFÉRENCE ==========

# Référence mouflon pour consensus (TRÈS IMPORTANT pour éviter biais)
MOUFLON_REF="/home/plstenge/genomes/Ovis_aries_mitochondrion/Ovis_aries_mitochondrion.fa"
MOUFLON_GENBANK="NC_001941"  # Accession GenBank si besoin de télécharger

# Si références mouton moderne pour phylogénie
SHEEP_REFS_DIR="${BASE_DIR}/13_phylogeny/00_references"
mkdir -p "$SHEEP_REFS_DIR"

echo ""
echo "Configuration génomes de référence:"
echo "  Mouflon (consensus): $MOUFLON_REF"
echo "  Dossier références: $SHEEP_REFS_DIR"

# Vérifier/télécharger référence mouflon si manquante
if [[ ! -f "$MOUFLON_REF" ]]; then
    echo ""
    echo "⚠ Référence mouflon manquante, téléchargement depuis NCBI..."
    mkdir -p "$(dirname "$MOUFLON_REF")"
    
    if command -v efetch &>/dev/null; then
        efetch -db nucleotide -id "$MOUFLON_GENBANK" -format fasta > "$MOUFLON_REF"
        echo "✓ Référence mouflon téléchargée"
    else
        echo "✗ Impossible de télécharger sans Entrez Direct"
        echo "  À faire manuellement: https://www.ncbi.nlm.nih.gov/nuccore/$MOUFLON_GENBANK"
        exit 1
    fi
fi

# Indexer référence
if [[ ! -f "${MOUFLON_REF}.fai" ]]; then
    echo "Indexation référence mouflon..."
    samtools faidx "$MOUFLON_REF"
fi

################################################################################
# ÉTAPE 9.1: EXTRACTION DES READS MAPPÉS DE MOUTON
################################################################################

echo ""
echo "=========================================="
echo "ÉTAPE 9.1: Extraction des reads mouton"
echo "=========================================="
echo ""

BAM_EXTRACTION_DIR="${PHYLOGENY_BASE}/01_bam_extraction"
SHEEP_BAMS=()

if [[ -d "$DAMAGE_DIR" ]]; then
    # Récupérer tous les BAM mouton depuis étape 8
    for BAM_FILE in "$DAMAGE_DIR"/Ovis_aries/*.sorted.bam; do
        if [[ -f "$BAM_FILE" ]]; then
            SAMPLE_NAME=$(basename "$BAM_FILE" .sorted.bam)
            
            echo "→ $SAMPLE_NAME"
            
            # Vérifier index
            if [[ ! -f "${BAM_FILE}.bai" ]]; then
                echo "  Création index..."
                samtools index "$BAM_FILE"
            fi
            
            # Statistiques
            TOTAL_READS=$(samtools view -c "$BAM_FILE")
            MAPPED_READS=$(samtools view -c -F 4 "$BAM_FILE")
            
            if [[ $TOTAL_READS -gt 0 ]]; then
                MAPPING_PCT=$(echo "scale=1; $MAPPED_READS * 100 / $TOTAL_READS" | bc)
            else
                MAPPING_PCT=0
            fi
            
            echo "  Statistiques: $MAPPED_READS/$TOTAL_READS lus ($MAPPING_PCT%)"
            
            if [[ $MAPPED_READS -lt 100 ]]; then
                echo "  ⚠ Trop peu de reads mappés, skipped"
                continue
            fi
            
            # Vérifier couverture mitochondriale
            COVERAGE=$(samtools depth "$BAM_FILE" 2>/dev/null | awk '{sum+=$3; count++} END {if (count>0) print int(sum/count); else print 0}')
            echo "  Couverture mitochondriale moyenne: ${COVERAGE}X"
            
            if [[ $COVERAGE -lt 3 ]]; then
                echo "  ⚠ Couverture trop basse (<3X), peut créer consensus partiel"
            fi
            
            # Copier BAM vers répertoire phylogénie
            cp "$BAM_FILE" "${BAM_EXTRACTION_DIR}/${SAMPLE_NAME}.bam"
            cp "${BAM_FILE}.bai" "${BAM_EXTRACTION_DIR}/${SAMPLE_NAME}.bam.bai"
            
            SHEEP_BAMS+=("${BAM_EXTRACTION_DIR}/${SAMPLE_NAME}.bam")
            
            echo "  ✓ BAM préparé"
        fi
    done
else
    echo "✗ Répertoire mapDamage non trouvé: $DAMAGE_DIR/Ovis_aries"
    exit 1
fi

N_SAMPLES=${#SHEEP_BAMS[@]}
echo ""
echo "✓ $N_SAMPLES BAM mouton préparés"

if [[ $N_SAMPLES -eq 0 ]]; then
    echo "ERREUR: Aucun BAM mouton valide trouvé"
    exit 1
fi

################################################################################
# ÉTAPE 9.1b: OPTION POOLING (si couverture faible)
################################################################################

ENABLE_POOLING=true  # Mettre à false pour désactiver

if [[ "$ENABLE_POOLING" == "true" && ${#SHEEP_BAMS[@]} -gt 1 ]]; then
    echo ""
    echo "=========================================="
    echo "ÉTAPE 9.1b: Création échantillon poolé"
    echo "=========================================="
    echo ""
    
    POOLED_BAM="${BAM_EXTRACTION_DIR}/POOLED_all_sheep.bam"
    
    echo "Pooling de ${#SHEEP_BAMS[@]} échantillons..."
    samtools merge -@ 8 -f "$POOLED_BAM" "${SHEEP_BAMS[@]}"
    samtools index "$POOLED_BAM"
    
    # Vérifier couverture du pooled
    POOLED_COV=$(samtools depth "$POOLED_BAM" | awk '{sum+=$3; count++} END {printf "%.1f", sum/count}')
    echo "✓ BAM poolé créé: couverture ${POOLED_COV}X"
    
    # Ajouter à la liste
    SHEEP_BAMS=("$POOLED_BAM")
    
    echo ""
fi


################################################################################
# ÉTAPE 9.2: GÉNÉRATION DES CONSENSUS MITOCHONDRIAUX
################################################################################

echo ""
echo "=========================================="
echo "ÉTAPE 9.2: Génération consensus mtDNA"
echo "=========================================="
echo ""

CONSENSUS_DIR="${PHYLOGENY_BASE}/02_consensus"
CONSENSUS_MULTIFASTA="${CONSENSUS_DIR}/all_samples_consensus_mtDNA.fasta"
CONSENSUS_STATS="${PHYLOGENY_BASE}/logs/consensus_stats.tsv"

echo -e "Sample\tConsensus_length\tMissing_bases\tN_count\tCoverage_mean" > "$CONSENSUS_STATS"

for BAM_FILE in "${SHEEP_BAMS[@]}"; do
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)
    
    echo "→ $SAMPLE_NAME"
    echo "  Appel de variants avec bcftools..."
    
    # Appeler variants avec mpileup
    # Paramètres: Q30 base quality, q30 mapping quality
    samtools view -b -F 4 "$BAM_FILE" | \
    bcftools mpileup -f "$MOUFLON_REF" \
        -Q 30 -q 30 \
        -d 10000 \
        - 2>"${PHYLOGENY_BASE}/logs/${SAMPLE_NAME}_mpileup.log" | \
    bcftools call -mv \
        -Oz \
        -o "${CONSENSUS_DIR}/${SAMPLE_NAME}_variants.vcf.gz" 2>>"${PHYLOGENY_BASE}/logs/${SAMPLE_NAME}_bcftools.log"
    
    # Indexer VCF
    bcftools index "${CONSENSUS_DIR}/${SAMPLE_NAME}_variants.vcf.gz"
    echo "  ✓ Variants appelés"
    
    # Générer consensus
    echo "  Génération du consensus..."
    bcftools consensus -f "$MOUFLON_REF" \
        "${CONSENSUS_DIR}/${SAMPLE_NAME}_variants.vcf.gz" \
        > "${CONSENSUS_DIR}/${SAMPLE_NAME}_consensus_tmp.fa"
    
    # Renommer header avec nom d'échantillon
    sed -i "s/>.*/>$SAMPLE_NAME/" "${CONSENSUS_DIR}/${SAMPLE_NAME}_consensus_tmp.fa"
    
    # Analyser consensus
    CONS_LENGTH=$(awk '/^>/ {next} {print length}' "${CONSENSUS_DIR}/${SAMPLE_NAME}_consensus_tmp.fa" | head -1)
    N_COUNT=$(grep -o 'N' "${CONSENSUS_DIR}/${SAMPLE_NAME}_consensus_tmp.fa" | wc -l)
    MISSING_BASES=$(awk '/^>/ {next} {gsub(/[atgcn]/,""); print length}' "${CONSENSUS_DIR}/${SAMPLE_NAME}_consensus_tmp.fa" | head -1)
    
    # Calculer couverture moyenne depuis depth
    COV_MEAN=$(samtools depth "$BAM_FILE" 2>/dev/null | \
        awk '{sum+=$3; count++} END {if (count>0) printf "%.1f", sum/count; else print "0"}')
    
    echo "  Consensus length: $CONS_LENGTH bp (N: $N_COUNT, missing: $MISSING_BASES)"
    echo "  Couverture moyenne: ${COV_MEAN}X"
    echo -e "${SAMPLE_NAME}\t${CONS_LENGTH}\t${MISSING_BASES}\t${N_COUNT}\t${COV_MEAN}" >> "$CONSENSUS_STATS"
    
    # Copier vers fichier final
    cat "${CONSENSUS_DIR}/${SAMPLE_NAME}_consensus_tmp.fa" >> "$CONSENSUS_MULTIFASTA"
    rm "${CONSENSUS_DIR}/${SAMPLE_NAME}_consensus_tmp.fa"
    
    echo "  ✓ Consensus généré"
done

echo ""
echo "✓ Fichier multi-FASTA consensus: $CONSENSUS_MULTIFASTA"
echo "  Statistiques: $CONSENSUS_STATS"
echo ""
cat "$CONSENSUS_STATS"

################################################################################
# ÉTAPE 9.3: TÉLÉCHARGEMENT RÉFÉRENCES MODERNES
################################################################################

echo ""
echo "=========================================="
echo "ÉTAPE 9.3: Références races mouton modernes"
echo "=========================================="
echo ""

REFERENCES_FASTA="${SHEEP_REFS_DIR}/sheep_reference_breeds_mtDNA.fasta"

# Vérifier si références existent déjà
if [[ ! -f "$REFERENCES_FASTA" ]]; then
    echo "Téléchargement génomes mitochondriaux mouton depuis NCBI..."
    echo "(Cette opération peut prendre plusieurs minutes)"
    
    # Télécharger via Entrez Direct
    if command -v esearch &>/dev/null; then
        esearch -db nucleotide \
            -query '"Ovis aries"[Organism] AND mitochondrion[Filter] AND complete genome' | \
            efetch -format fasta > "$REFERENCES_FASTA" 2>>"${PHYLOGENY_BASE}/logs/ncbi_download.log"
        
        N_REFS=$(grep -c "^>" "$REFERENCES_FASTA" 2>/dev/null || echo 0)
        echo "✓ $N_REFS séquences de référence téléchargées"
    else
        echo "⚠ Entrez Direct non disponible"
        echo "  Références doivent être téléchargées manuellement et placées dans:"
        echo "  $REFERENCES_FASTA"
        echo ""
        echo "  Instructions:"
        echo "  1. Aller sur: https://www.ncbi.nlm.nih.gov/nuccore/"
        echo "  2. Chercher: \"Ovis aries\"[Organism] AND mitochondrion[Filter] AND complete genome"
        echo "  3. Télécharger en format FASTA"
        echo ""
        
        # Créer fichier placeholder
        touch "$REFERENCES_FASTA"
        echo "  Placeholder créé: $REFERENCES_FASTA"
    fi
else
    N_REFS=$(grep -c "^>" "$REFERENCES_FASTA" 2>/dev/null || echo 0)
    echo "✓ Références existantes: $N_REFS séquences"
fi

# Vérifier qu'on a au minimum quelques références
if [[ ! -s "$REFERENCES_FASTA" ]]; then
    echo ""
    echo "⚠ AVERTISSEMENT: Pas de références modernes"
    echo "  Vous devez télécharger des séquences mouton de référence"
    echo "  depuis NCBI GenBank pour une phylogénie significative"
    echo ""
    echo "  Format attendu: FASTA avec > 10 séquences représentant"
    echo "  différentes races (européennes, asiatiques, moyen-oriental)"
fi

################################################################################
# ÉTAPE 10.1: ALIGNEMENT MULTIPLE
################################################################################

echo ""
echo "=========================================="
echo "ÉTAPE 10.1: Alignement multiple MAFFT"
echo "=========================================="
echo ""

ALIGNMENT_DIR="${PHYLOGENY_BASE}/03_alignment"
COMBINED_FASTA="${ALIGNMENT_DIR}/combined_ancient_modern_mtDNA.fasta"
ALIGNED_FASTA="${ALIGNMENT_DIR}/aligned_mtDNA.fasta"

# Combiner anciens + modernes
echo "Combinaison des séquences..."
cat "$CONSENSUS_MULTIFASTA" > "$COMBINED_FASTA"

if [[ -s "$REFERENCES_FASTA" ]]; then
    cat "$REFERENCES_FASTA" >> "$COMBINED_FASTA"
    TOTAL_SEQS=$(grep -c "^>" "$COMBINED_FASTA")
    echo "✓ Fichier combiné: $TOTAL_SEQS séquences"
else
    echo "⚠ Fichier références vide, alignement avec échantillons anciens uniquement"
fi

# Vérifier qu'on a au moins 2 séquences
NSEQS=$(grep -c "^>" "$COMBINED_FASTA")
if [[ $NSEQS -lt 2 ]]; then
    echo "✗ ERREUR: Besoin d'au moins 2 séquences pour alignement"
    exit 1
fi

# Alignement MAFFT
echo ""
echo "Alignement MAFFT (--auto)..."
echo "  Nombre de séquences: $NSEQS"

if [[ $NSEQS -gt 200 ]]; then
    MAFFT_MODE="--parttree"
    echo "  Mode: parttree (grand nombre de séquences)"
else
    MAFFT_MODE="--auto"
    echo "  Mode: auto"
fi

time mafft $MAFFT_MODE \
    --thread "$THREADS" \
    "$COMBINED_FASTA" > "$ALIGNED_FASTA" 2>>"${PHYLOGENY_BASE}/logs/mafft.log"

ALIGN_LENGTH=$(head -n 2 "$ALIGNED_FASTA" | tail -n 1 | wc -c)
echo "✓ Alignement généré: ${ALIGN_LENGTH} bp"
echo "  Fichier: $ALIGNED_FASTA"

################################################################################
# ÉTAPE 10.2: PHYLOGÉNIE RAXML
################################################################################

echo ""
echo "=========================================="
echo "ÉTAPE 10.2: Phylogénie RAxML (Maximum Likelihood)"
echo "=========================================="
echo ""

PHYLOGENY_DIR="${PHYLOGENY_BASE}/04_phylogeny"
RAXML_RUN_NAME="sheep_mtDNA_ml_$(date +%Y%m%d)"

echo "Construction arbre phylogénétique..."
echo "  Modèle: GTRGAMMA"
echo "  Bootstraps: 100"
echo "  Threads: $THREADS"
echo ""

# Vérifier RAxML disponible
if ! command -v raxmlHPC-PTHREADS &>/dev/null && ! command -v raxml &>/dev/null; then
    echo "⚠ RAxML non trouvé"
    echo "  Installation depuis bioconda:"
    echo "    conda install -c bioconda raxml"
    echo ""
else
    # Rechercher la bonne version
    RAXML_CMD="raxmlHPC-PTHREADS"
    if ! command -v "$RAXML_CMD" &>/dev/null; then
        RAXML_CMD="raxml"
    fi
    
    cd "$PHYLOGENY_DIR"
    
    echo "Exécution: $RAXML_CMD ..."
    time "$RAXML_CMD" \
        -f a \
        -x 12345 \
        -p 12345 \
        -# 100 \
        -m GTRGAMMA \
        -s "$ALIGNED_FASTA" \
        -n "$RAXML_RUN_NAME" \
        -T "$THREADS" \
        2>&1 | tee -a "${PHYLOGENY_BASE}/logs/raxml.log"
    
    if [[ -f "RAxML_bipartitions.${RAXML_RUN_NAME}" ]]; then
        echo "✓ Arbre RAxML généré: RAxML_bipartitions.${RAXML_RUN_NAME}"
        echo "  Bootstrap support calculé et appliqué aux branches"
    else
        echo "✗ Erreur lors de la construction arbre RAxML"
    fi
    
    cd - > /dev/null
fi

################################################################################
# ÉTAPE 10.3: SETUP PHYLOGÉNIE BEAST (optionnel)
################################################################################

echo ""
echo "=========================================="
echo "ÉTAPE 10.3: Préparation BEAST (optionnel)"
echo "=========================================="
echo ""

BEAST_DIR="${PHYLOGENY_BASE}/04_phylogeny/beast_input"
mkdir -p "$BEAST_DIR"

echo "Pour exécuter BEAST (phylogénie Bayésienne avec dates):"
echo ""
echo "1. Préparer le fichier d'alignement:"
cp "$ALIGNED_FASTA" "${BEAST_DIR}/aligned_mtDNA.fasta"
echo "   ✓ Copié: ${BEAST_DIR}/aligned_mtDNA.fasta"
echo ""

echo "2. Ouvrir BEAUti (interface graphique):"
echo "   beauti"
echo ""

echo "3. Étapes dans BEAUti:"
echo "   a) Charger l'alignement FASTA"
echo "   b) Tips panel: Si dates disponibles, activer 'Use tip dates'"
echo "   c) Site Model: HKY ou GTR + Gamma"
echo "   d) Clock Model: Strict clock (ou Relaxed si variation suspectée)"
echo "   e) Tree Prior: Coalescent Constant Size"
echo "   f) MCMC: Chain length 10,000,000"
echo "   g) Générer XML (ex: sheep_mtDNA.xml)"
echo ""

echo "4. Sauvegarder XML dans: $BEAST_DIR/"
echo ""

echo "5. Exécuter BEAST:"
echo "   beast -beagle_off -threads $THREADS sheep_mtDNA.xml"
echo ""

# Créer template de fichier de métadonnées pour BEAST
cat > "${BEAST_DIR}/README_BEAST.txt" <<'EOF'
GUIDE D'UTILISATION DE BEAST POUR ADN ANCIEN

=== PRÉPARATION ===
1. Alignement: aligned_mtDNA.fasta (créé automatiquement)

2. Métadonnées (si applicable): Créer fichier CSV avec dates
   Format: Sample,Date_BP (avant présent)
   Exemple:
   sample_01,3200
   sample_02,3150
   Gotland,0  # Moderne = 0 BP

=== CONFIGURATION BEAST ===
Site Model:
  - Substitution: HKY ou GTR (mouton mtDNA = GTR recommandé)
  - Site heterogeneity: Gamma, 4 catégories
  - Codon positions: 1,2,3

Clock Model:
  - Strict clock: Si taux évolution constant
  - Relaxed clock (lognormal): Si variation des taux
    Recommandé pour ADN ancien + références modernes

Tree Prior:
  - Coalescent Constant Size: Populationsimple
  - Coalescent Exponential Growth: Si expansion suspectée
  
MCMC:
  - Chain length: 10,000,000 minimum (100M pour publication)
  - Log parameters: tous les 1000
  - Log trees: tous les 10,000
  
=== VÉRIFICATION CONVERGENCE ===
Avec Tracer:
  tracer sheep_mtDNA.log
  
Critères OK:
  - ESS (Effective Sample Size) ≥ 200 pour tous paramètres
  - Traces visuellement stables

=== RÉSULTATS ===
TreeAnnotator (burnin 10%):
  treeannotator -burnin 10 -heights median sheep_mtDNA.trees sheep_mtDNA_MCC.tree

=== VISUALISATION ===
FigTree:
  figtree sheep_mtDNA_MCC.tree
EOF

echo "✓ Template README créé: ${BEAST_DIR}/README_BEAST.txt"

################################################################################
# ÉTAPE 10.4: PRÉPARATION MITOTOOLPY (HAPLOGROUPE)
################################################################################

echo ""
echo "=========================================="
echo "ÉTAPE 10.4: Assignation haplogroupes (optionnel)"
echo "=========================================="
echo ""

HAPLOGROUP_DIR="${PHYLOGENY_BASE}/05_visualization/haplogroups"
mkdir -p "$HAPLOGROUP_DIR"

echo "Pour assigner les haplogroupes mitochondriaux:"
echo ""
echo "1. Installer MitoToolPy:"
echo "   git clone https://github.com/kizcas/MitoToolPy.git"
echo "   cd MitoToolPy"
echo ""

echo "2. Exécuter pour chaque consensus:"
echo "   python mitotoolpy-seq.py -s sheep -r whole -i consensus.fasta -o results.txt"
echo ""

# Créer script helper
cat > "${HAPLOGROUP_DIR}/run_mitotoolpy.sh" <<'EOFSCRIPT'
#!/bin/bash

# Script helper pour assignation haplogroupes

CONSENSUS_DIR="../../../02_consensus"
RESULTS_DIR="."

# Installer MitoToolPy si manquant
if [[ ! -d "MitoToolPy" ]]; then
    echo "Clonage MitoToolPy..."
    git clone https://github.com/kizcas/MitoToolPy.git
fi

cd MitoToolPy

echo "Assignation haplogroupes..."

for CONS_FILE in $CONSENSUS_DIR/*_consensus_tmp.fa; do
    SAMPLE=$(basename "$CONS_FILE" _consensus_tmp.fa)
    echo "  → $SAMPLE"
    
    python mitotoolpy-seq.py \
        -s sheep \
        -r whole \
        -i "$CONS_FILE" \
        -o "${RESULTS_DIR}/${SAMPLE}_haplogroup.txt"
done

echo "✓ Terminé"
echo "Résultats dans: $RESULTS_DIR/"
EOFSCRIPT

chmod +x "${HAPLOGROUP_DIR}/run_mitotoolpy.sh"
echo "✓ Script helper créé: ${HAPLOGROUP_DIR}/run_mitotoolpy.sh"

################################################################################
# ÉTAPE 10.5: VISUALISATION FIGTREE
################################################################################

echo ""
echo "=========================================="
echo "ÉTAPE 10.5: Préparation visualisation"
echo "=========================================="
echo ""

VIZ_DIR="${PHYLOGENY_BASE}/05_visualization"
TREE_FILES=(
    "${PHYLOGENY_DIR}/RAxML_bipartitions.${RAXML_RUN_NAME}"
)

echo "Pour visualiser l'arbre phylogénétique avec FigTree:"
echo ""
echo "  figtree ${PHYLOGENY_DIR}/RAxML_bipartitions.${RAXML_RUN_NAME}"
echo ""

echo "Annotations recommandées dans FigTree:"
echo "  1. Node Labels → Afficher bootstrap (label)"
echo "  2. Tip Labels → Gras, couleur par origine"
echo "  3. Branch Width → Par bootstrap value"
echo "  4. Scale Bar → Ajouter"
echo "  5. Tree Root → Ovis_orientalis (outgroup)"
echo ""

################################################################################
# RÉSUMÉ FINAL
################################################################################

echo ""
echo "=========================================="
echo "PIPELINE PHYLOGÉNIE TERMINÉ"
echo "=========================================="
echo ""

echo "📊 RÉSUMÉS:"
echo "  Statistiques consensus: ${PHYLOGENY_BASE}/logs/consensus_stats.tsv"
echo ""

echo "📁 FICHIERS GÉNÉRÉS:"
echo "  Consensus FASTA: $CONSENSUS_MULTIFASTA"
echo "  Alignement: $ALIGNED_FASTA"
echo "  Arbre RAxML: ${PHYLOGENY_DIR}/RAxML_bipartitions.${RAXML_RUN_NAME}"
echo ""

echo "📈 ÉTAPES SUIVANTES:"
echo ""
echo "1. HAPLOGROUPE ASSIGNMENT (optionnel mais recommandé):"
echo "   cd ${HAPLOGROUP_DIR}"
echo "   bash run_mitotoolpy.sh"
echo ""

echo "2. PHYLOGÉNIE BAYÉSIENNE BEAST (si dates disponibles):"
echo "   cd ${BEAST_DIR}"
echo "   # 1. Ouvrir BEAUti, créer XML"
echo "   # 2. Lancer: beast -beagle_off -threads 8 sheep_mtDNA.xml"
echo "   # 3. Vérifier: tracer sheep_mtDNA.log"
echo "   # 4. Consensus: treeannotator -burnin 10 sheep_mtDNA.trees sheep_mtDNA_MCC.tree"
echo ""

echo "3. VISUALISATION FIGTREE:"
echo "   figtree ${PHYLOGENY_DIR}/RAxML_bipartitions.${RAXML_RUN_NAME}"
echo ""

echo "📝 LOG COMPLET: $LOGFILE"
echo ""

# Afficher statistiques consensus
echo "STATISTIQUES CONSENSUS:"
echo ""
column -t "${PHYLOGENY_BASE}/logs/consensus_stats.tsv"
echo ""

echo "=========================================="
echo "✓ Pipeline réussi"
echo "=========================================="
echo ""

conda deactivate

# Notification email
{
    echo "Pipeline phylogénie terminé avec succès"
    echo ""
    echo "Résultats:"
    echo "  Consensus: $CONSENSUS_MULTIFASTA"
    echo "  Alignement: $ALIGNED_FASTA"
    echo "  Arbre RAxML: ${PHYLOGENY_DIR}/RAxML_bipartitions.${RAXML_RUN_NAME}"
    echo ""
    echo "Visualiser: figtree ${PHYLOGENY_DIR}/RAxML_bipartitions.${RAXML_RUN_NAME}"
} | mail -s "Pipeline Phylogénie Mouton - Terminé" pierrelouis.stenger@gmail.com 2>/dev/null || true

exit 0
