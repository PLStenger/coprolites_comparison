#!/bin/bash

#SBATCH --job-name=genome_sheep_race_identification
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=200G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/phylogeny_identification_genomic.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/phylogeny_identification_genomic.out"

# ========== INITIALISATION CONDA POUR SLURM ==========
echo "Initialisation de conda..."
module load conda/4.12.0
source ~/.bashrc
conda activate genomics_analysis

# Charger R si nécessaire (vérifie si Rscript est disponible)
which Rscript > /dev/null 2>&1 || module load R/4.3.0

# ========== CONFIG ==========
# Utiliser les BAMs de paired_end (choix arbitraire, tu peux changer pour single_end_R1)
BAM_DIR="/home/plstenge/coprolites_comparison/12_mapdamage/unmerged_reads/Ovis_aries/paired_end"
REF="/home/plstenge/genomes/Ovis_aries/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa" 
OUT_DIR="/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis"

SAMPLES=(cop408 cop410 cop412 cop414)

mkdir -p $OUT_DIR/{01_snps,02_merge,03_analysis}

# ========== VÉRIFICATION PRÉALABLE ==========
echo "=========================================="
echo "VÉRIFICATION DES FICHIERS"
echo "=========================================="

echo "Référence: $REF"
if [[ ! -f "$REF" ]]; then
    echo "⚠️  ERREUR: Génome de référence introuvable !"
    exit 1
fi

echo "BAMs à traiter:"
for SAMPLE in "${SAMPLES[@]}"; do
    BAM_FILE="${BAM_DIR}/${SAMPLE}_unmerged_Ovis_aries.sorted.bam"
    if [[ -f "$BAM_FILE" ]]; then
        echo "  ✓ $SAMPLE: $(ls -lh $BAM_FILE | awk '{print $5}')"
    else
        echo "  ✗ $SAMPLE: INTROUVABLE à $BAM_FILE"
        exit 1
    fi
done

echo ""
echo "=========================================="
echo "DÉBUT DU PIPELINE"
echo "=========================================="

# ========== ÉTAPE 1 : SNP calling par échantillon ==========
echo ""
echo "ÉTAPE 1: SNP calling (bcftools)"
echo "----------------------------------------"

for SAMPLE in "${SAMPLES[@]}"; do
    BAM_FILE="${BAM_DIR}/${SAMPLE}_unmerged_Ovis_aries.sorted.bam"
    VCF_OUT="${OUT_DIR}/01_snps/${SAMPLE}.vcf.gz"
    
    echo "  Processing $SAMPLE..."
    
    bcftools mpileup \
        -f $REF \
        -Q 30 -q 30 \
        -a FORMAT/AD,FORMAT/DP \
        $BAM_FILE | \
    bcftools call \
        -m -v \
        -Oz -o $VCF_OUT
    
    if [[ $? -ne 0 ]]; then
        echo "  ✗ Échec du SNP calling pour $SAMPLE"
        exit 1
    fi
    
    bcftools index $VCF_OUT
    
    # Compter les SNPs
    N_SNPS=$(bcftools view -H $VCF_OUT | wc -l)
    echo "    → $N_SNPS SNPs détectés"
done

echo "  ✓ SNP calling terminé pour tous les échantillons"

# ========== ÉTAPE 2 : Pooler les échantillons ==========
echo ""
echo "ÉTAPE 2: Pooling samples"
echo "----------------------------------------"

# Créer liste des VCFs
VCF_LIST=""
for SAMPLE in "${SAMPLES[@]}"; do
    VCF_LIST="$VCF_LIST ${OUT_DIR}/01_snps/${SAMPLE}.vcf.gz"
done

bcftools merge \
    $VCF_LIST \
    -Oz -o ${OUT_DIR}/02_merge/POOLED_ancient.vcf.gz

if [[ $? -ne 0 ]]; then
    echo "  ✗ Échec du merge"
    exit 1
fi

bcftools index ${OUT_DIR}/02_merge/POOLED_ancient.vcf.gz

# Filtrer SNPs haute qualité
echo "  Filtrage des SNPs (QUAL>30, DP>5)..."
bcftools view \
    -i 'QUAL>30 && DP>5' \
    ${OUT_DIR}/02_merge/POOLED_ancient.vcf.gz \
    -Oz -o ${OUT_DIR}/02_merge/POOLED_ancient_filtered.vcf.gz

bcftools index ${OUT_DIR}/02_merge/POOLED_ancient_filtered.vcf.gz

# Statistiques
bcftools stats ${OUT_DIR}/02_merge/POOLED_ancient_filtered.vcf.gz \
    > ${OUT_DIR}/02_merge/stats.txt

N_SNPS_POOLED=$(bcftools view -H ${OUT_DIR}/02_merge/POOLED_ancient_filtered.vcf.gz | wc -l)
echo "  ✓ SNPs poolés après filtrage: $N_SNPS_POOLED"

if [[ $N_SNPS_POOLED -eq 0 ]]; then
    echo "  ⚠️  ATTENTION: Aucun SNP détecté après filtrage !"
    echo "  Vérifier la profondeur de séquençage et la qualité du mapping"
    exit 1
fi

# ========== ÉTAPE 3 : Préparer références modernes ==========
echo ""
echo "ÉTAPE 3: Préparation des données modernes"
echo "----------------------------------------"

# Copier et renommer le VCF moderne
if [[ ! -f "${OUT_DIR}/02_merge/modern_sheep_1000genomes.vcf.gz" ]]; then
    echo "  Copie du VCF moderne..."
    cp /home/plstenge/coprolites_comparison/15_modern_sheep_vcf/ISGC_SNP50_Breedv2.vcf.gz \
       ${OUT_DIR}/02_merge/modern_sheep_1000genomes.vcf.gz
    
    # Créer l'index
    bcftools index ${OUT_DIR}/02_merge/modern_sheep_1000genomes.vcf.gz
else
    echo "  ✓ VCF moderne déjà présent"
fi

# Vérifier le nombre de SNPs modernes
N_SNPS_MODERN=$(bcftools view -H ${OUT_DIR}/02_merge/modern_sheep_1000genomes.vcf.gz | wc -l)
echo "  → SNPs modernes: $N_SNPS_MODERN"

# ========== ÉTAPE 4 : Merge ancien + moderne ==========
echo ""
echo "ÉTAPE 4: Merge ancient + modern"
echo "----------------------------------------"

# Merger uniquement sur les positions communes
bcftools isec \
    -p ${OUT_DIR}/02_merge/isec_temp \
    ${OUT_DIR}/02_merge/POOLED_ancient_filtered.vcf.gz \
    ${OUT_DIR}/02_merge/modern_sheep_1000genomes.vcf.gz

echo "  Positions communes détectées, merge en cours..."

bcftools merge \
    --force-samples \
    ${OUT_DIR}/02_merge/POOLED_ancient_filtered.vcf.gz \
    ${OUT_DIR}/02_merge/modern_sheep_1000genomes.vcf.gz \
    -Oz -o ${OUT_DIR}/02_merge/combined_ancient_modern.vcf.gz

bcftools index ${OUT_DIR}/02_merge/combined_ancient_modern.vcf.gz

N_SNPS_COMBINED=$(bcftools view -H ${OUT_DIR}/02_merge/combined_ancient_modern.vcf.gz | wc -l)
echo "  ✓ SNPs dans le merge: $N_SNPS_COMBINED"

# Convertir en PLINK format
echo "  Conversion en format PLINK..."
plink --vcf ${OUT_DIR}/02_merge/combined_ancient_modern.vcf.gz \
      --sheep \
      --make-bed \
      --allow-extra-chr \
      --out ${OUT_DIR}/03_analysis/combined

if [[ $? -ne 0 ]]; then
    echo "  ✗ Échec de la conversion PLINK"
    exit 1
fi

echo "  ✓ Conversion PLINK terminée"

# ========== ÉTAPE 5 : PCA ==========
echo ""
echo "ÉTAPE 5: PCA analysis"
echo "----------------------------------------"

plink --bfile ${OUT_DIR}/03_analysis/combined \
      --pca 10 \
      --allow-extra-chr \
      --out ${OUT_DIR}/03_analysis/pca

if [[ $? -ne 0 ]]; then
    echo "  ✗ Échec de la PCA"
    exit 1
fi

# Plot PCA (R) - SANS guillemets simples autour de EOF
echo "  Génération du plot PCA..."
Rscript <<EOF
library(ggplot2)

# Lire les données PCA
pca <- read.table("${OUT_DIR}/03_analysis/pca.eigenvec", header=FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:10))

# Identifier les échantillons anciens
pca\$Type <- ifelse(grepl("cop", pca\$IID), "Ancient", "Modern")

# Plot
p <- ggplot(pca, aes(PC1, PC2, color=Type, shape=Type)) +
  geom_point(size=3, alpha=0.7) +
  scale_color_manual(values=c("Ancient"="red", "Modern"="blue")) +
  theme_bw(base_size=14) +
  labs(title="PCA: Ancient vs Modern Sheep",
       x="PC1", y="PC2") +
  theme(legend.position="top")

ggsave("${OUT_DIR}/03_analysis/pca_plot.pdf", plot=p, width=10, height=8)
ggsave("${OUT_DIR}/03_analysis/pca_plot.png", plot=p, width=10, height=8, dpi=300)

cat("✓ PCA plot sauvegardé\n")
EOF

if [[ $? -eq 0 ]]; then
    echo "  ✓ PCA terminée"
else
    echo "  ⚠️  Erreur dans le script R (plot non généré)"
fi

# ========== ÉTAPE 6 : ADMIXTURE ==========
echo ""
echo "ÉTAPE 6: ADMIXTURE analysis"
echo "----------------------------------------"

# Vérifier si admixture est disponible
which admixture > /dev/null 2>&1
if [[ $? -ne 0 ]]; then
    echo "  ⚠️  ADMIXTURE non disponible, étape ignorée"
else
    cd ${OUT_DIR}/03_analysis
    
    for K in {2..10}; do
        echo "  Running K=$K..."
        admixture --cv \
            combined.bed \
            $K \
            -j8 \
            > admixture_K${K}.log 2>&1
    done
    
    # Choisir meilleur K (lowest CV error)
    echo ""
    echo "  Cross-validation errors:"
    grep "CV error" admixture_K*.log | sort -k4 -n
    
    cd - > /dev/null
    echo "  ✓ ADMIXTURE terminé"
fi

# ========== ÉTAPE 7 : f-statistics (ADMIXTOOLS) ==========
echo ""
echo "ÉTAPE 7: f-statistics (ADMIXTOOLS)"
echo "----------------------------------------"

# Vérifier si ADMIXTOOLS est disponible
which qp3Pop > /dev/null 2>&1
if [[ $? -ne 0 ]]; then
    echo "  ⚠️  ADMIXTOOLS non disponible, étape ignorée"
    echo "  (Installation: conda install -c bioconda admixtools)"
else
    echo "  ⚠️  Configuration manuelle requise pour ADMIXTOOLS"
    echo "  Fichiers de configuration (.par) à créer selon ton design expérimental"
    # qp3Pop -p qp3Pop.par > ${OUT_DIR}/03_analysis/f3_results.txt
    # qpDstat -p qpDstat.par > ${OUT_DIR}/03_analysis/D_statistics.txt
fi

# ========== RÉSUMÉ ==========
echo ""
echo "=========================================="
echo "✓ PIPELINE TERMINÉ"
echo "=========================================="
echo ""
echo "Résumé des résultats:"
echo "  • SNPs anciens (après filtrage): $N_SNPS_POOLED"
echo "  • SNPs modernes: $N_SNPS_MODERN"
echo "  • SNPs combinés: $N_SNPS_COMBINED"
echo ""
echo "Fichiers générés:"
echo "  • VCFs: ${OUT_DIR}/02_merge/"
echo "  • PLINK: ${OUT_DIR}/03_analysis/combined.{bed,bim,fam}"
echo "  • PCA: ${OUT_DIR}/03_analysis/pca_plot.{pdf,png}"
echo "  • ADMIXTURE: ${OUT_DIR}/03_analysis/admixture_K*.log"
echo ""
echo "=========================================="
