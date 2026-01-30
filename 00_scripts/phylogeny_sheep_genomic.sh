#!/bin/bash

#SBATCH --job-name=genome_sheep_race_identification
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=200G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/phylogeny_identification_genomic.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/phylogeny_identification_genomic.out"


# ========== CONFIG ==========
BAM_DIR="/home/plstenge/coprolites_comparison/11_map_to_sheep_genome"
REF="/home/plstenge/genomes/Ovis_aries/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa" 
OUT_DIR="/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis"

SAMPLES=(cop408 cop410 cop412 cop414)

mkdir -p $OUT_DIR/{01_snps,02_merge,03_analysis}

# ========== ÉTAPE 1 : SNP calling par échantillon ==========
echo "ÉTAPE 1: SNP calling (bcftools)"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "  Processing $SAMPLE..."
    
    bcftools mpileup \
        -f $REF \
        -Q 30 -q 30 \
        -a FORMAT/AD,FORMAT/DP \
        ${BAM_DIR}/${SAMPLE}_sorted_RG_dedup_realigned_q30.bam | \
    bcftools call \
        -m -v \
        -Oz -o ${OUT_DIR}/01_snps/${SAMPLE}.vcf.gz
    
    bcftools index ${OUT_DIR}/01_snps/${SAMPLE}.vcf.gz
done

# ========== ÉTAPE 2 : Pooler les échantillons ==========
echo "ÉTAPE 2: Pooling samples"

bcftools merge \
    ${OUT_DIR}/01_snps/*.vcf.gz \
    -Oz -o ${OUT_DIR}/02_merge/POOLED_ancient.vcf.gz

bcftools index ${OUT_DIR}/02_merge/POOLED_ancient.vcf.gz

# Filtrer SNPs haute qualité
bcftools view \
    -i 'QUAL>30 && DP>5' \
    ${OUT_DIR}/02_merge/POOLED_ancient.vcf.gz \
    -Oz -o ${OUT_DIR}/02_merge/POOLED_ancient_filtered.vcf.gz

# Statistiques
bcftools stats ${OUT_DIR}/02_merge/POOLED_ancient_filtered.vcf.gz \
    > ${OUT_DIR}/02_merge/stats.txt

echo "✓ SNPs poolés: $(bcftools view -H ${OUT_DIR}/02_merge/POOLED_ancient_filtered.vcf.gz | wc -l)"

# ========== ÉTAPE 3 : Télécharger références modernes ==========
echo "ÉTAPE 3: Download modern reference genomes"

# Télécharger VCF moutons modernes (ex: 1000 Sheep Genomes Project)
# http://www.sheephapmap.org/

scp -r /home/plstenge/coprolites_comparison/15_modern_sheep_vcf/ISGC_SNP50_Breedv2.vcf.gz ${OUT_DIR}/02_merge/ 
mv ${OUT_DIR}/02_merge/ISGC_SNP50_Breedv2.vcf.gz ${OUT_DIR}/02_merge/modern_sheep_1000genomes.vcf.gz 

# ========== ÉTAPE 4 : Merge ancien + moderne ==========
echo "ÉTAPE 4: Merge ancient + modern"

bcftools merge \
    ${OUT_DIR}/02_merge/POOLED_ancient_filtered.vcf.gz \
    ${OUT_DIR}/02_merge/modern_sheep_1000genomes.vcf.gz \
    -Oz -o ${OUT_DIR}/02_merge/combined_ancient_modern.vcf.gz

# Convertir en PLINK format
plink --vcf ${OUT_DIR}/02_merge/combined_ancient_modern.vcf.gz \
      --sheep \
      --make-bed \
      --out ${OUT_DIR}/03_analysis/combined

# ========== ÉTAPE 5 : PCA ==========
echo "ÉTAPE 5: PCA analysis"

plink --bfile ${OUT_DIR}/03_analysis/combined \
      --pca 10 \
      --out ${OUT_DIR}/03_analysis/pca

# Plot PCA (R)
Rscript << 'EOF'
library(ggplot2)
pca <- read.table("${OUT_DIR}/03_analysis/pca.eigenvec", header=F)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:10))

# Identifier POOLED
pca$Type <- ifelse(grepl("POOLED", pca$IID), "Ancient", "Modern")

ggplot(pca, aes(PC1, PC2, color=Type, shape=Type)) +
  geom_point(size=3) +
  scale_color_manual(values=c("Ancient"="red", "Modern"="blue")) +
  theme_bw() +
  ggtitle("PCA: Ancient vs Modern Sheep")
ggsave("${OUT_DIR}/03_analysis/pca_plot.pdf")
EOF

# ========== ÉTAPE 6 : ADMIXTURE ==========
echo "ÉTAPE 6: ADMIXTURE analysis"

for K in {2..10}; do
    admixture --cv \
        ${OUT_DIR}/03_analysis/combined.bed \
        $K \
        -j8 \
        > ${OUT_DIR}/03_analysis/admixture_K${K}.log
done

# Choisir meilleur K (lowest CV error)
grep "CV error" ${OUT_DIR}/03_analysis/admixture_K*.log

# ========== ÉTAPE 7 : f-statistics (introgression) ==========
echo "ÉTAPE 7: f-statistics (ADMIXTOOLS)"

# Préparer fichiers pour ADMIXTOOLS
# ... (format spécifique)

qp3Pop -p qp3Pop.par > f3_results.txt
qpDstat -p qpDstat.par > D_statistics.txt

echo "✓ PIPELINE TERMINÉ"
