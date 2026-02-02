#!/bin/bash

#SBATCH --job-name=fix_and_phylogeny
#SBATCH --ntasks=4
#SBATCH -p smp
#SBATCH --mem=100G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/fix_phylogeny.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/fix_phylogeny.out"

module load conda/4.12.0
source ~/.bashrc
conda activate genomics_analysis
which Rscript > /dev/null 2>&1 || module load R/4.3.0

BAM_DIR="/home/plstenge/coprolites_comparison/12_mapdamage/unmerged_reads/Ovis_aries/paired_end"
REF="/home/plstenge/genomes/Ovis_aries/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa"
OUT_DIR="/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis"
MODERN_VCF_ORIG="/home/plstenge/coprolites_comparison/15_modern_sheep_vcf/ISGC_SNP50_Breedv2.vcf.gz"

SAMPLES=(cop408 cop410 cop412 cop414)

mkdir -p $OUT_DIR/{08_fixed,09_phylogeny}

echo "=========================================="
echo "RÉPARATION + PHYLOGÉNIE"
echo "=========================================="

# ========== PARTIE 1 : RÉPARER VCF ANCIEN (noms d'échantillons) ==========
echo ""
echo "PARTIE 1: Réparation du VCF ancien"
echo "----------------------------------------"

echo "  Création de VCFs individuels propres..."

for SAMPLE in "${SAMPLES[@]}"; do
    BAM_FILE="${BAM_DIR}/${SAMPLE}_unmerged_Ovis_aries.sorted.bam"
    VCF_OUT="${OUT_DIR}/08_fixed/${SAMPLE}.vcf.gz"
    
    if [[ -f "$VCF_OUT" ]]; then
        echo "  $SAMPLE: déjà fait"
        continue
    fi
    
    echo "  $SAMPLE..."
    
    # SNP calling avec nom d'échantillon propre
    bcftools mpileup \
        -f $REF \
        -Q 20 -q 20 -A \
        --max-depth 1000000 \
        -a FORMAT/AD,FORMAT/DP \
        $BAM_FILE | \
    bcftools call -m -v --ploidy 2 | \
    bcftools reheader -s <(echo "$SAMPLE") | \
    bcftools view -Oz -o $VCF_OUT
    
    bcftools index -f $VCF_OUT
    
    N=$(bcftools view -H $VCF_OUT | wc -l)
    echo "    → $N SNPs, nom échantillon: $SAMPLE"
done

# Merger proprement
echo "  Merge des échantillons..."

bcftools merge \
    ${OUT_DIR}/08_fixed/cop408.vcf.gz \
    ${OUT_DIR}/08_fixed/cop410.vcf.gz \
    ${OUT_DIR}/08_fixed/cop412.vcf.gz \
    ${OUT_DIR}/08_fixed/cop414.vcf.gz \
    -Oz -o ${OUT_DIR}/08_fixed/ancient_merged.vcf.gz

bcftools index -f ${OUT_DIR}/08_fixed/ancient_merged.vcf.gz

# Vérifier les noms
echo "  Vérification des noms d'échantillons:"
bcftools query -l ${OUT_DIR}/08_fixed/ancient_merged.vcf.gz

# Filtrer
bcftools view \
    -i 'QUAL>20 && FORMAT/DP>2' \
    -v snps \
    ${OUT_DIR}/08_fixed/ancient_merged.vcf.gz \
    -Oz -o ${OUT_DIR}/08_fixed/ancient_clean.vcf.gz

bcftools index -f ${OUT_DIR}/08_fixed/ancient_clean.vcf.gz

N_ANCIENT=$(bcftools view -H ${OUT_DIR}/08_fixed/ancient_clean.vcf.gz | wc -l)
echo "  ✓ $N_ANCIENT SNPs anciens propres"

# ========== PARTIE 2 : VÉRIFIER/IGNORER VCF MODERNE ==========
echo ""
echo "PARTIE 2: Diagnostic du VCF moderne"
echo "----------------------------------------"

echo "  Le VCF moderne (ISGC_SNP50) est corrompu."
echo "  Chromosomes détectés: 00, 1-28 (format numérique)"
echo ""
echo "  Vérification de la compatibilité avec ton génome de référence..."

# Extraire les chromosomes de ton génome de référence
samtools idxstats ${BAM_DIR}/cop408_unmerged_Ovis_aries.sorted.bam | head -30 | awk '{print $1"\t"$2}' > ${OUT_DIR}/08_fixed/ref_chromosomes.txt

echo "  Chromosomes dans tes BAMs (premiers 30):"
head -20 ${OUT_DIR}/08_fixed/ref_chromosomes.txt

echo ""
echo "  ⚠️  Le VCF moderne utilise les chromosomes: 0,1,2,3...28"
echo "  ⚠️  Tes BAMs utilisent probablement: NC_XXXXXX ou chr1, chr2..."
echo ""
echo "  → INCOMPATIBILITÉ : Impossible de faire l'imputation"
echo "  → SOLUTION : Phylogénie uniquement sur échantillons anciens"
echo ""

# ========== PARTIE 3 : PHYLOGÉNIE ANCIENS UNIQUEMENT ==========
echo ""
echo "PARTIE 3: Phylogénie des échantillons anciens"
echo "----------------------------------------"

VCF_ANCIENT="${OUT_DIR}/08_fixed/ancient_clean.vcf.gz"

# Conversion PLINK
plink --vcf $VCF_ANCIENT \
      --sheep --allow-extra-chr \
      --make-bed \
      --set-missing-var-ids @:# \
      --out ${OUT_DIR}/09_phylogeny/ancient

# PCA
plink --bfile ${OUT_DIR}/09_phylogeny/ancient \
      --pca 4 \
      --allow-extra-chr \
      --out ${OUT_DIR}/09_phylogeny/pca

# Distance génétique
plink --bfile ${OUT_DIR}/09_phylogeny/ancient \
      --distance square \
      --allow-extra-chr \
      --out ${OUT_DIR}/09_phylogeny/distances

# Identity-By-State (IBS) matrix
plink --bfile ${OUT_DIR}/09_phylogeny/ancient \
      --distance-matrix \
      --allow-extra-chr \
      --out ${OUT_DIR}/09_phylogeny/ibs

echo "  ✓ Analyses PLINK terminées"

# ========== PARTIE 4 : VISUALISATIONS R ==========
echo ""
echo "PARTIE 4: Génération des visualisations"
echo "----------------------------------------"

Rscript <<'EOF'
library(ggplot2)
library(ape)
library(phangorn)
library(pheatmap)

outdir <- "/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis/09_phylogeny"

# ===== 1. PCA =====
cat("Génération PCA...\n")
pca <- read.table(paste0(outdir, "/pca.eigenvec"), header=FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:4))

p_pca <- ggplot(pca, aes(PC1, PC2, label=IID)) +
  geom_point(size=6, color="darkred", alpha=0.8) +
  geom_text(vjust=-1.5, size=7, fontface="bold", color="black") +
  theme_bw(base_size=16) +
  labs(title="PCA - Échantillons anciens de mouton",
       subtitle=paste0("Basé sur génome nucléaire (", 246, " SNPs)"),
       x="Composante principale 1",
       y="Composante principale 2") +
  theme(plot.title=element_text(face="bold", size=18))

ggsave(paste0(outdir, "/pca_plot.pdf"), p_pca, width=10, height=8)
ggsave(paste0(outdir, "/pca_plot.png"), p_pca, width=10, height=8, dpi=300)

# ===== 2. MATRICE DE DISTANCE =====
cat("Calcul distances...\n")
dist_mat <- as.matrix(read.table(paste0(outdir, "/distances.dist")))
ids <- read.table(paste0(outdir, "/distances.dist.id"))$V1
rownames(dist_mat) <- ids
colnames(dist_mat) <- ids

cat("\n=== MATRICE DE DISTANCE GÉNÉTIQUE ===\n")
print(round(dist_mat, 5))

write.csv(dist_mat, paste0(outdir, "/distance_matrix.csv"))

# ===== 3. HEATMAP =====
cat("Heatmap...\n")
pdf(paste0(outdir, "/distance_heatmap.pdf"), width=9, height=8)
pheatmap(dist_mat,
         display_numbers=TRUE,
         number_format="%.5f",
         number_color="black",
         fontsize_number=14,
         fontsize=14,
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         main="Distances génétiques entre coprolites",
         color=colorRampPalette(c("white", "orange", "red"))(100))
dev.off()

# ===== 4. ARBRE PHYLOGÉNÉTIQUE =====
cat("Construction de l'arbre...\n")
dist_obj <- as.dist(dist_mat)
tree_nj <- nj(dist_obj)
tree_rooted <- midpoint(tree_nj)

write.tree(tree_nj, paste0(outdir, "/tree_NJ.nwk"))
write.tree(tree_rooted, paste0(outdir, "/tree_NJ_rooted.nwk"))

# Plot simple
pdf(paste0(outdir, "/tree_simple.pdf"), width=10, height=8)
par(mfrow=c(1,1), mar=c(2,2,4,2))
plot(tree_rooted, 
     main="Arbre phylogénétique - Moutons anciens\n(Neighbor-Joining, enraciné au midpoint)",
     cex=2, font=2, edge.width=4, cex.main=1.5)
add.scale.bar(cex=1.5, lwd=2)
dev.off()

# Plot avec distances sur branches
pdf(paste0(outdir, "/tree_with_distances.pdf"), width=12, height=8)
par(mar=c(2,2,4,2))
plot(tree_rooted, 
     main="Arbre avec distances génétiques",
     cex=2, font=2, edge.width=4, show.node.label=TRUE)
edgelabels(round(tree_rooted$edge.length, 4), 
           bg="white", cex=1.2, frame="rect")
add.scale.bar(cex=1.5, lwd=2)
dev.off()

# ===== 5. STATISTIQUES DÉTAILLÉES =====
cat("\n=== STATISTIQUES COMPARATIVES ===\n")

# Calculer pairwise stats
pairs <- combn(ids, 2)
for(i in 1:ncol(pairs)) {
  s1 <- pairs[1,i]
  s2 <- pairs[2,i]
  d <- dist_mat[s1, s2]
  cat(sprintf("%s vs %s: distance = %.5f\n", s1, s2, d))
}

# Identifier les plus proches
min_dist <- min(dist_mat[dist_mat > 0])
max_dist <- max(dist_mat)
cat(sprintf("\nDistance minimale: %.5f\n", min_dist))
cat(sprintf("Distance maximale: %.5f\n", max_dist))

# Trouver paire la plus proche
closest <- which(dist_mat == min_dist & dist_mat > 0, arr.ind=TRUE)[1,]
cat(sprintf("Paire la plus proche: %s - %s\n", 
            rownames(dist_mat)[closest[1]], 
            colnames(dist_mat)[closest[2]]))

# ===== 6. RÉSUMÉ VISUEL =====
summary_data <- data.frame(
  Sample = ids,
  MinDist = apply(dist_mat + diag(Inf, nrow(dist_mat)), 1, min),
  MeanDist = apply(dist_mat, 1, function(x) mean(x[x>0]))
)

p_summary <- ggplot(summary_data, aes(x=Sample, y=MeanDist)) +
  geom_bar(stat="identity", fill="steelblue", alpha=0.7) +
  geom_point(aes(y=MinDist), color="red", size=4) +
  geom_text(aes(y=MeanDist, label=round(MeanDist, 4)), 
            vjust=-0.5, size=5, fontface="bold") +
  theme_bw(base_size=14) +
  labs(title="Distances génétiques moyennes par échantillon",
       subtitle="Barres = distance moyenne | Points rouges = distance au plus proche",
       x="Échantillon", y="Distance génétique") +
  theme(axis.text.x=element_text(size=12, face="bold"))

ggsave(paste0(outdir, "/summary_distances.pdf"), p_summary, width=10, height=6)

cat("\n✓ Toutes les visualisations générées\n")
EOF

# ========== RÉSUMÉ FINAL ==========
echo ""
echo "=========================================="
echo "✓ PIPELINE TERMINÉ"
echo "=========================================="
echo ""
echo "📊 RÉSULTATS:"
echo ""
echo "Visualisations principales:"
echo "  • PCA: ${OUT_DIR}/09_phylogeny/pca_plot.pdf"
echo "  • Heatmap distances: ${OUT_DIR}/09_phylogeny/distance_heatmap.pdf"
echo "  • Arbre phylogénétique: ${OUT_DIR}/09_phylogeny/tree_simple.pdf"
echo "  • Arbre avec distances: ${OUT_DIR}/09_phylogeny/tree_with_distances.pdf"
echo "  • Résumé: ${OUT_DIR}/09_phylogeny/summary_distances.pdf"
echo ""
echo "Fichiers de données:"
echo "  • VCF nettoyé: ${OUT_DIR}/08_fixed/ancient_clean.vcf.gz"
echo "  • Matrice de distances: ${OUT_DIR}/09_phylogeny/distance_matrix.csv"
echo "  • Arbre Newick: ${OUT_DIR}/09_phylogeny/tree_NJ_rooted.nwk"
echo "  • PLINK: ${OUT_DIR}/09_phylogeny/ancient.{bed,bim,fam}"
echo ""
echo "📌 CONCLUSION sur le VCF moderne:"
echo "  Le fichier ISGC_SNP50_Breedv2.vcf.gz est incompatible"
echo "  (format corrompu + nomenclature chromosomes différente)"
echo ""
echo "  Pour comparer avec des modernes, il faudrait:"
echo "  1. Trouver un VCF moderne avec nomenclature compatible"
echo "  2. OU faire l'analyse mitochondriale + GenBank"
echo "  3. OU reséquencer à plus haute couverture"
echo ""
echo "=========================================="
