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

# Charger R si nécessaire
which Rscript > /dev/null 2>&1 || module load R/4.3.0

# ========== CONFIG ==========
BAM_DIR="/home/plstenge/coprolites_comparison/12_mapdamage/unmerged_reads/Ovis_aries/paired_end"
REF="/home/plstenge/genomes/Ovis_aries/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa" 
OUT_DIR="/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis"

SAMPLES=(cop408 cop410 cop412 cop414)

mkdir -p $OUT_DIR/{01_snps,02_merge,03_analysis,04_logs}

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
        N_READS=$(samtools view -c $BAM_FILE)
        echo "  ✓ $SAMPLE: $(ls -lh $BAM_FILE | awk '{print $5}') ($N_READS reads)"
    else
        echo "  ✗ $SAMPLE: INTROUVABLE"
        exit 1
    fi
done

echo ""
echo "=========================================="
echo "DÉBUT DU PIPELINE - ADN ANCIEN"
echo "=========================================="

# ========== ÉTAPE 1 : SNP calling ULTRA-PERMISSIF pour ADN ancien ==========
echo ""
echo "ÉTAPE 1: SNP calling avec paramètres ADN ancien"
echo "----------------------------------------"

for SAMPLE in "${SAMPLES[@]}"; do
    BAM_FILE="${BAM_DIR}/${SAMPLE}_unmerged_Ovis_aries.sorted.bam"
    VCF_OUT="${OUT_DIR}/01_snps/${SAMPLE}.vcf.gz"
    
    echo "  Processing $SAMPLE..."
    
    # PARAMÈTRES ADAPTÉS ADN ANCIEN:
    # -Q 20 (au lieu de 30) : qualité base moins stricte
    # -q 20 (au lieu de 30) : qualité mapping moins stricte  
    # --max-depth 1000000 : pas de limite de profondeur
    # -A : ne pas ignorer les anomalies de reads
    bcftools mpileup \
        -f $REF \
        -Q 20 \
        -q 20 \
        -A \
        --max-depth 1000000 \
        -a FORMAT/AD,FORMAT/DP,INFO/AD \
        $BAM_FILE | \
    bcftools call \
        -m \
        -v \
        --ploidy 2 \
        -Oz -o $VCF_OUT
    
    if [[ $? -ne 0 ]]; then
        echo "  ✗ Échec du SNP calling pour $SAMPLE"
        exit 1
    fi
    
    bcftools index -f $VCF_OUT
    
    # Statistiques détaillées
    N_SNPS_RAW=$(bcftools view -H $VCF_OUT | wc -l)
    echo "    → $N_SNPS_RAW SNPs bruts détectés"
    
    # Stats de profondeur
    bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' $VCF_OUT | \
        awk '{sum+=$3; n++} END {if(n>0) print "    → Profondeur moyenne: "sum/n}' 
done

echo "  ✓ SNP calling terminé pour tous les échantillons"

# ========== ÉTAPE 2 : Pooler avec filtres TRÈS PERMISSIFS ==========
echo ""
echo "ÉTAPE 2: Pooling avec filtres ADN ancien"
echo "----------------------------------------"

# Créer liste des VCFs
VCF_LIST=""
for SAMPLE in "${SAMPLES[@]}"; do
    VCF_LIST="$VCF_LIST ${OUT_DIR}/01_snps/${SAMPLE}.vcf.gz"
done

bcftools merge \
    --force-samples \
    $VCF_LIST \
    -Oz -o ${OUT_DIR}/02_merge/POOLED_ancient_raw.vcf.gz

bcftools index -f ${OUT_DIR}/02_merge/POOLED_ancient_raw.vcf.gz

N_SNPS_RAW=$(bcftools view -H ${OUT_DIR}/02_merge/POOLED_ancient_raw.vcf.gz | wc -l)
echo "  → SNPs bruts poolés: $N_SNPS_RAW"

# Filtrage PERMISSIF pour ADN ancien:
# QUAL>20 (au lieu de 30)
# FORMAT/DP>2 (au lieu de 5) - accepter profondeur très faible
# Exclure indels si trop bruités
echo "  Filtrage permissif (QUAL>20, FORMAT/DP>2)..."

bcftools view \
    -i 'QUAL>20 && FORMAT/DP>2' \
    -v snps \
    ${OUT_DIR}/02_merge/POOLED_ancient_raw.vcf.gz \
    -Oz -o ${OUT_DIR}/02_merge/POOLED_ancient_filtered.vcf.gz

bcftools index -f ${OUT_DIR}/02_merge/POOLED_ancient_filtered.vcf.gz

# Statistiques
bcftools stats ${OUT_DIR}/02_merge/POOLED_ancient_filtered.vcf.gz \
    > ${OUT_DIR}/02_merge/stats_filtered.txt

N_SNPS_FILTERED=$(bcftools view -H ${OUT_DIR}/02_merge/POOLED_ancient_filtered.vcf.gz | wc -l)
echo "  ✓ SNPs après filtrage: $N_SNPS_FILTERED"

if [[ $N_SNPS_FILTERED -eq 0 ]]; then
    echo ""
    echo "  ⚠️  ATTENTION: Toujours 0 SNPs !"
    echo ""
    echo "  Diagnostic de la profondeur par échantillon:"
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF="${OUT_DIR}/01_snps/${SAMPLE}.vcf.gz"
        echo "  $SAMPLE:"
        bcftools query -f '[%SAMPLE\t%DP\n]' $VCF | \
            awk '{sum+=$2; n++; if($2>0) npos++} END {
                print "    - Positions totales: "n
                print "    - Positions avec couverture >0: "npos
                if(n>0) print "    - Profondeur moyenne: "sum/n
            }'
    done
    
    echo ""
    echo "  Solutions possibles:"
    echo "  1. Merger paired_end + single_end_R1 pour augmenter la couverture"
    echo "  2. Ne garder que cop414 (seul avec SNPs)"
    echo "  3. Utiliser une approche pseudohaploid (sampling aléatoire)"
    echo "  4. Travailler sur génome mitochondrial au lieu du nucléaire"
    echo ""
    
    # Créer un rapport diagnostic
    cat > ${OUT_DIR}/04_logs/diagnostic.txt <<DIAG
DIAGNOSTIC - Couverture insuffisante pour ADN nucléaire ancien
================================================================

Reads mappés:
$(for S in "${SAMPLES[@]}"; do
    BAM="${BAM_DIR}/${S}_unmerged_Ovis_aries.sorted.bam"
    echo "  $S: $(samtools view -c $BAM) reads"
done)

SNPs détectés (bruts):
$(for S in "${SAMPLES[@]}"; do
    VCF="${OUT_DIR}/01_snps/${S}.vcf.gz"
    echo "  $S: $(bcftools view -H $VCF 2>/dev/null | wc -l) SNPs"
done)

Recommandations:
- Les BAMs sont très petits (34K à 1.5M) → couverture <0.001X
- Pour phylogénie nucléaire, il faut idéalement >1X de couverture
- Options:
  1. Analyser l'ADN mitochondrial (meilleure préservation)
  2. Utiliser des captures ciblées (si disponibles)
  3. Approche pseudohaploid + ADMIXTOOLS pour faible couverture
DIAG
    
    cat ${OUT_DIR}/04_logs/diagnostic.txt
    
    echo ""
    echo "  Veux-tu continuer avec cop414 uniquement ? (y/n)"
    echo "  Ou lancer l'analyse mitochondriale ? (Dis-moi)"
    
    exit 1
fi

# ========== ÉTAPE 3 : Intersection avec données modernes ==========
echo ""
echo "ÉTAPE 3: Préparation et intersection avec données modernes"
echo "----------------------------------------"

# Copier VCF moderne si nécessaire
if [[ ! -f "${OUT_DIR}/02_merge/modern_sheep.vcf.gz" ]]; then
    echo "  Copie du VCF moderne..."
    cp /home/plstenge/coprolites_comparison/15_modern_sheep_vcf/ISGC_SNP50_Breedv2.vcf.gz \
       ${OUT_DIR}/02_merge/modern_sheep.vcf.gz
    bcftools index -f ${OUT_DIR}/02_merge/modern_sheep.vcf.gz
fi

N_SNPS_MODERN=$(bcftools view -H ${OUT_DIR}/02_merge/modern_sheep.vcf.gz | wc -l)
echo "  → SNPs modernes: $N_SNPS_MODERN"

# Trouver positions communes (CRITIQUE pour ADN ancien)
echo "  Intersection des positions..."
bcftools isec \
    -n=2 \
    -w1 \
    ${OUT_DIR}/02_merge/POOLED_ancient_filtered.vcf.gz \
    ${OUT_DIR}/02_merge/modern_sheep.vcf.gz \
    -Oz -o ${OUT_DIR}/02_merge/ancient_common_pos.vcf.gz

bcftools isec \
    -n=2 \
    -w2 \
    ${OUT_DIR}/02_merge/POOLED_ancient_filtered.vcf.gz \
    ${OUT_DIR}/02_merge/modern_sheep.vcf.gz \
    -Oz -o ${OUT_DIR}/02_merge/modern_common_pos.vcf.gz

bcftools index -f ${OUT_DIR}/02_merge/ancient_common_pos.vcf.gz
bcftools index -f ${OUT_DIR}/02_merge/modern_common_pos.vcf.gz

N_COMMON=$(bcftools view -H ${OUT_DIR}/02_merge/ancient_common_pos.vcf.gz | wc -l)
echo "  → Positions communes: $N_COMMON"

if [[ $N_COMMON -eq 0 ]]; then
    echo "  ⚠️  Aucune position commune ! Arrays SNP incompatibles."
    echo "  → Le VCF moderne est probablement sur un array de SNPs fixes"
    echo "  → Les échantillons anciens ont trop peu de couverture pour ces positions"
    exit 1
fi

# Merger
bcftools merge \
    --force-samples \
    ${OUT_DIR}/02_merge/ancient_common_pos.vcf.gz \
    ${OUT_DIR}/02_merge/modern_common_pos.vcf.gz \
    -Oz -o ${OUT_DIR}/02_merge/combined_ancient_modern.vcf.gz

bcftools index -f ${OUT_DIR}/02_merge/combined_ancient_modern.vcf.gz

echo "  ✓ Merge terminé: $N_COMMON SNPs"

# ========== ÉTAPE 4 : Conversion PLINK ==========
echo ""
echo "ÉTAPE 4: Conversion PLINK"
echo "----------------------------------------"

plink --vcf ${OUT_DIR}/02_merge/combined_ancient_modern.vcf.gz \
      --sheep \
      --make-bed \
      --allow-extra-chr \
      --set-missing-var-ids @:# \
      --out ${OUT_DIR}/03_analysis/combined

if [[ $? -ne 0 ]]; then
    echo "  ✗ Échec PLINK"
    exit 1
fi

echo "  ✓ Conversion terminée"

# ========== ÉTAPE 5 : PCA ==========
echo ""
echo "ÉTAPE 5: PCA"
echo "----------------------------------------"

plink --bfile ${OUT_DIR}/03_analysis/combined \
      --pca 10 \
      --allow-extra-chr \
      --out ${OUT_DIR}/03_analysis/pca

Rscript <<EOF
library(ggplot2)

pca <- read.table("${OUT_DIR}/03_analysis/pca.eigenvec", header=FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:10))

pca\$Type <- ifelse(grepl("cop", pca\$IID), "Ancient", "Modern")

p <- ggplot(pca, aes(PC1, PC2, color=Type, shape=Type)) +
  geom_point(size=3, alpha=0.7) +
  scale_color_manual(values=c("Ancient"="red", "Modern"="blue")) +
  theme_bw(base_size=14) +
  labs(title="PCA: Ancient vs Modern Sheep", x="PC1", y="PC2")

ggsave("${OUT_DIR}/03_analysis/pca_plot.pdf", plot=p, width=10, height=8)
ggsave("${OUT_DIR}/03_analysis/pca_plot.png", plot=p, width=10, height=8, dpi=300)
EOF

echo "  ✓ PCA terminée"

# ========== ÉTAPE 6 : ADMIXTURE ==========
echo ""
echo "ÉTAPE 6: ADMIXTURE"
echo "----------------------------------------"

which admixture > /dev/null 2>&1
if [[ $? -ne 0 ]]; then
    echo "  ⚠️  ADMIXTURE non disponible"
else
    cd ${OUT_DIR}/03_analysis
    for K in {2..6}; do
        echo "  K=$K..."
        admixture --cv combined.bed $K -j8 > admixture_K${K}.log 2>&1
    done
    grep "CV error" admixture_K*.log | sort -k4 -n
    cd - > /dev/null
fi

echo ""
echo "=========================================="
echo "✓ PIPELINE TERMINÉ"
echo "=========================================="
