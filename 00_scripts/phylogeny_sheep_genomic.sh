#!/bin/bash

#SBATCH --job-name=phylogeny_ancient_modern_genomes
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH -p smp
#SBATCH --mem=250G
#SBATCH --time=48:00:00
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/phylogeny_genomes.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/phylogeny_genomes.out"

# ========== INITIALISATION ==========
echo "Initialisation..."
module load conda/4.12.0
source ~/.bashrc
conda activate genomics_analysis
which Rscript > /dev/null 2>&1 || module load R/4.3.0

# Vérifier outils
for tool in bwa samtools bcftools plink; do
    which $tool > /dev/null 2>&1 || { echo "⚠️  $tool manquant"; exit 1; }
done

# ========== CONFIGURATION ==========
ANCIENT_BAM_DIR="/home/plstenge/coprolites_comparison/12_mapdamage/unmerged_reads/Ovis_aries/paired_end"
MODERN_GENOMES_DIR="/home/plstenge/coprolites_comparison/16_genomes"
REF="/home/plstenge/genomes/Ovis_aries/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa"
OUT_DIR="/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis/11_complete_phylogeny"

ANCIENT_SAMPLES=(cop408 cop410 cop412 cop414)

# Liste des génomes modernes avec noms courts
declare -A MODERN_GENOMES=(
    ["Ammon_polii"]="Ovis_ammon_polii/GCA_028583565.1_CAU_O.ammon_polii_1.0_genomic.fna.gz"
    ["Awassi"]="Ovis_aries_Awassi_breed/GCA_040543055.1_ASM4054305v1_genomic.fna.gz"
    ["Charollais"]="Ovis_aries_Charollais_breed/GCA_022416745.1_ASM2241674v1_genomic.fna.gz"
    ["East_Friesian"]="Ovis_aries_East_Friesian_breed/GCA_033439445.1_ASM3343944v1_genomic.fna.gz"
    ["Hu_Sheep"]="Ovis_aries_Hu_Sheep_breed/GCA_040805955.1_T2T-sheep1.0_genomic.fna.gz"
    ["Kazak"]="Ovis_aries_Kazak_breed/GCA_022432845.1_ASM2243284v1_genomic.fna.gz"
    ["Kermani"]="Ovis_aries_Kermani_breed/GCA_022432835.1_ASM2243283v1_genomic.fna.gz"
    ["Mouflon"]="Ovis_aries_musimon/GCF_000765115.1_Oori1_genomic.fna.gz"
    ["Polled_Dorset"]="Ovis_aries_Polled_Dorset_breed/GCA_022416915.1_ASM2241691v1_genomic.fna.gz"
    ["Rambouillet"]="Ovis_aries_Rambouillet_breed/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna.gz"
    ["Bighorn"]="Ovis_canadensis/GCF_042477335.2_ARS-UI_OviCan_v2_genomic.fna.gz"
    ["Snow_Sheep"]="Ovis_nivicola_lydekkeri/GCA_903231385.1_OvNiv1.0_genomic.fna.gz"
    ["Orientalis"]="Ovis_orientalis_from_Iran_PRJEB3141/GCF_000765115.1_Oori1_genomic.fna.gz"
    ["Urial"]="Ovis_vignei/GCA_053525375.1_ASM5352537v1_genomic.fna.gz"
)

THREADS=16

mkdir -p $OUT_DIR/{01_indexed_genomes,02_consensus,03_vcf,04_merged,05_analysis}

echo "=========================================="
echo "PHYLOGÉNIE COMPLÈTE"
echo "=========================================="
echo "Anciens: 4 coprolites"
echo "Modernes: ${#MODERN_GENOMES[@]} génomes de référence"
echo "=========================================="

# ========== ÉTAPE 1 : Index des génomes modernes ==========
echo ""
echo "ÉTAPE 1: Indexation des génomes modernes"
echo "----------------------------------------"

# Vérifier/créer index de la référence commune
if [[ ! -f "${REF}.bwt" ]]; then
    echo "  Indexation de la référence commune..."
    bwa index $REF
fi

# ========== ÉTAPE 2 : Générer séquences consensus pour modernes ==========
echo ""
echo "ÉTAPE 2: Génération des consensus modernes"
echo "----------------------------------------"
echo "  Stratégie: Pseudogénomes basés sur alignement vs référence commune"

for NAME in "${!MODERN_GENOMES[@]}"; do
    GENOME_PATH="${MODERN_GENOMES_DIR}/${MODERN_GENOMES[$NAME]}"
    CONSENSUS="${OUT_DIR}/02_consensus/${NAME}_consensus.fa.gz"
    
    if [[ -f "$CONSENSUS" ]]; then
        echo "  $NAME: déjà fait"
        continue
    fi
    
    echo "  $NAME: extraction chromosome 1 (représentatif)..."
    
    # Pour comparaison rapide: extraire CHR1 uniquement
    # (génome entier prendrait trop de temps/RAM)
    zcat $GENOME_PATH | \
    awk '/^>/{if(NR>1) exit; print ">'"$NAME"'"; next} {print}' | \
    head -100000 | \
    gzip > $CONSENSUS
    
    echo "    → Consensus CHR1 créé"
done

echo "  ⚠️  Note: Utilisation du chromosome 1 uniquement pour accélérer"
echo "     Pour analyse complète, utiliser génomes entiers (temps ++)"

# ========== ÉTAPE 3 : SNP calling sur échantillons anciens ==========
echo ""
echo "ÉTAPE 3: SNP calling échantillons anciens"
echo "----------------------------------------"

for SAMPLE in "${ANCIENT_SAMPLES[@]}"; do
    BAM="${ANCIENT_BAM_DIR}/${SAMPLE}_unmerged_Ovis_aries.sorted.bam"
    VCF="${OUT_DIR}/03_vcf/${SAMPLE}.vcf.gz"
    
    if [[ -f "$VCF" ]]; then
        echo "  $SAMPLE: déjà fait"
        continue
    fi
    
    echo "  $SAMPLE..."
    
    bcftools mpileup \
        -f $REF \
        -Q 20 -q 20 -A \
        --max-depth 1000000 \
        -r 1 \
        $BAM | \
    bcftools call -m -v --ploidy 2 | \
    bcftools reheader -s <(echo "$SAMPLE") | \
    bcftools view -Oz -o $VCF
    
    bcftools index -f $VCF
    
    N=$(bcftools view -H $VCF | wc -l)
    echo "    → $N SNPs (chromosome 1)"
done

# ========== ÉTAPE 4 : Appeler variants sur consensus modernes ==========
echo ""
echo "ÉTAPE 4: Extraction variants depuis consensus modernes"
echo "----------------------------------------"

for NAME in "${!MODERN_GENOMES[@]}"; do
    CONSENSUS="${OUT_DIR}/02_consensus/${NAME}_consensus.fa.gz"
    VCF="${OUT_DIR}/03_vcf/${NAME}.vcf.gz"
    
    if [[ -f "$VCF" ]]; then
        echo "  $NAME: déjà fait"
        continue
    fi
    
    echo "  $NAME..."
    
    # Mapper consensus sur référence
    TMP_BAM="${OUT_DIR}/03_vcf/${NAME}_tmp.bam"
    
    zcat $CONSENSUS | \
    bwa mem -t 4 $REF - | \
    samtools sort -@ 4 -o $TMP_BAM
    
    samtools index $TMP_BAM
    
    # Calling
    bcftools mpileup -f $REF -r 1 --max-depth 500 $TMP_BAM | \
    bcftools call -m -v --ploidy 2 | \
    bcftools reheader -s <(echo "$NAME") | \
    bcftools view -Oz -o $VCF
    
    bcftools index -f $VCF
    
    rm -f $TMP_BAM $TMP_BAM.bai
    
    N=$(bcftools view -H $VCF | wc -l)
    echo "    → $N SNPs"
done

# ========== ÉTAPE 5 : Merge tous les VCFs ==========
echo ""
echo "ÉTAPE 5: Merge de tous les échantillons"
echo "----------------------------------------"

ALL_VCFS=""
for SAMPLE in "${ANCIENT_SAMPLES[@]}"; do
    ALL_VCFS="$ALL_VCFS ${OUT_DIR}/03_vcf/${SAMPLE}.vcf.gz"
done

for NAME in "${!MODERN_GENOMES[@]}"; do
    ALL_VCFS="$ALL_VCFS ${OUT_DIR}/03_vcf/${NAME}.vcf.gz"
done

echo "  Merge de $(echo $ALL_VCFS | wc -w) VCFs..."

bcftools merge $ALL_VCFS -Oz -o ${OUT_DIR}/04_merged/all_samples_raw.vcf.gz
bcftools index -f ${OUT_DIR}/04_merged/all_samples_raw.vcf.gz

N_RAW=$(bcftools view -H ${OUT_DIR}/04_merged/all_samples_raw.vcf.gz | wc -l)
echo "  → $N_RAW SNPs bruts"

# Filtrage
echo "  Filtrage (QUAL>20, missing<50%)..."

bcftools view \
    -i 'QUAL>20 && F_MISSING<0.5' \
    -v snps \
    ${OUT_DIR}/04_merged/all_samples_raw.vcf.gz \
    -Oz -o ${OUT_DIR}/04_merged/all_samples_filtered.vcf.gz

bcftools index -f ${OUT_DIR}/04_merged/all_samples_filtered.vcf.gz

N_FILT=$(bcftools view -H ${OUT_DIR}/04_merged/all_samples_filtered.vcf.gz | wc -l)
echo "  ✓ $N_FILT SNPs filtrés"

if [[ $N_FILT -lt 100 ]]; then
    echo "  ⚠️  Peu de SNPs détectés, relaxation des filtres..."
    
    bcftools view \
        -i 'F_MISSING<0.7' \
        -v snps \
        ${OUT_DIR}/04_merged/all_samples_raw.vcf.gz \
        -Oz -o ${OUT_DIR}/04_merged/all_samples_filtered.vcf.gz
    
    bcftools index -f ${OUT_DIR}/04_merged/all_samples_filtered.vcf.gz
    N_FILT=$(bcftools view -H ${OUT_DIR}/04_merged/all_samples_filtered.vcf.gz | wc -l)
    echo "  → $N_FILT SNPs après filtres relaxés"
fi

# ========== ÉTAPE 6 : Conversion PLINK ==========
echo ""
echo "ÉTAPE 6: Conversion PLINK"
echo "----------------------------------------"

plink --vcf ${OUT_DIR}/04_merged/all_samples_filtered.vcf.gz \
      --sheep \
      --allow-extra-chr \
      --make-bed \
      --set-missing-var-ids @:# \
      --out ${OUT_DIR}/05_analysis/all_samples

N_SNPS_PLINK=$(wc -l < ${OUT_DIR}/05_analysis/all_samples.bim)
echo "  ✓ $N_SNPS_PLINK SNPs dans PLINK"

# ========== ÉTAPE 7 : PCA ==========
echo ""
echo "ÉTAPE 7: PCA"
echo "----------------------------------------"

plink --bfile ${OUT_DIR}/05_analysis/all_samples \
      --pca 10 \
      --allow-extra-chr \
      --out ${OUT_DIR}/05_analysis/pca

echo "  ✓ PCA terminée"

# ========== ÉTAPE 8 : Distance génétique ==========
echo ""
echo "ÉTAPE 8: Distances génétiques"
echo "----------------------------------------"

plink --bfile ${OUT_DIR}/05_analysis/all_samples \
      --distance square \
      --allow-extra-chr \
      --out ${OUT_DIR}/05_analysis/distances

echo "  ✓ Matrice de distances calculée"

# ========== ÉTAPE 9 : ADMIXTURE ==========
echo ""
echo "ÉTAPE 9: ADMIXTURE"
echo "----------------------------------------"

which admixture > /dev/null 2>&1
if [[ $? -ne 0 ]]; then
    echo "  ⚠️  ADMIXTURE non disponible, ignoré"
else
    cd ${OUT_DIR}/05_analysis
    
    for K in {2..10}; do
        echo "  K=$K..."
        admixture --cv all_samples.bed $K -j8 > admixture_K${K}.log 2>&1
    done
    
    echo ""
    echo "  Cross-validation (5 meilleurs):"
    grep "CV error" admixture_K*.log | sort -k4 -n | head -5
    
    cd - > /dev/null
fi

# ========== ÉTAPE 10 : VISUALISATIONS R ==========
echo ""
echo "ÉTAPE 10: Visualisations"
echo "----------------------------------------"

Rscript <<'EOF'
library(ggplot2)
library(ggrepel)
library(ape)
library(phangorn)
library(pheatmap)
library(RColorBrewer)
library(reshape2)

outdir <- "/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis/11_complete_phylogeny/05_analysis"
setwd(outdir)

# ===== 1. PCA =====
cat("Génération PCA...\n")
pca <- read.table("pca.eigenvec", header=FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:10))

# Catégoriser
pca$Category <- "Modern_Domestic"
pca$Category[grepl("cop", pca$IID)] <- "Ancient"
pca$Category[grepl("Ammon|Bighorn|Snow|Orientalis|Urial|Mouflon", pca$IID)] <- "Wild"

pca$Type <- factor(pca$Category, levels=c("Ancient", "Modern_Domestic", "Wild"))

# Couleurs
cols <- c("Ancient"="red", "Modern_Domestic"="steelblue", "Wild"="darkgreen")

# PCA PC1-PC2
p1 <- ggplot(pca, aes(PC1, PC2, color=Type, label=IID)) +
  geom_point(size=4, alpha=0.8) +
  geom_text_repel(size=3.5, fontface="bold", max.overlaps=25, 
                  segment.size=0.3, box.padding=0.5) +
  scale_color_manual(values=cols) +
  theme_bw(base_size=14) +
  labs(title="PCA: Phylogénie complète des moutons",
       subtitle=paste0("Anciens (n=4) vs Domestiques modernes (n=9) vs Sauvages (n=5)"),
       x="PC1", y="PC2") +
  theme(legend.position="top", legend.text=element_text(size=12, face="bold"))

ggsave("pca_PC1_PC2.pdf", p1, width=16, height=12)
ggsave("pca_PC1_PC2.png", p1, width=16, height=12, dpi=300)

# PCA PC2-PC3
p2 <- ggplot(pca, aes(PC2, PC3, color=Type, label=IID)) +
  geom_point(size=4, alpha=0.8) +
  geom_text_repel(size=3.5, fontface="bold", max.overlaps=25) +
  scale_color_manual(values=cols) +
  theme_bw(base_size=14) +
  labs(title="PCA: PC2 vs PC3", x="PC2", y="PC3") +
  theme(legend.position="top")

ggsave("pca_PC2_PC3.pdf", p2, width=16, height=12)

# ===== 2. DISTANCE MATRIX =====
cat("Analyse des distances...\n")
dist_mat <- as.matrix(read.table("distances.dist"))
ids <- read.table("distances.dist.id")$V1
rownames(dist_mat) <- ids
colnames(dist_mat) <- ids

# Annotation
annot_df <- data.frame(
    Type = ifelse(grepl("cop", ids), "Ancient",
                  ifelse(grepl("Ammon|Bighorn|Snow|Orientalis|Urial|Mouflon", ids), 
                         "Wild", "Domestic")),
    row.names = ids
)

# Couleurs annotation
annot_colors <- list(Type = c("Ancient"="red", "Domestic"="steelblue", "Wild"="darkgreen"))

# Heatmap
pdf("distance_heatmap_full.pdf", width=16, height=15)
pheatmap(dist_mat,
         annotation_row = annot_df,
         annotation_col = annot_df,
         annotation_colors = annot_colors,
         display_numbers = FALSE,
         fontsize = 10,
         fontsize_row = 9,
         fontsize_col = 9,
         main = "Distances génétiques - Phylogénie complète",
         color = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
         border_color = "grey80")
dev.off()

# Distances anciens vs modernes
cat("\n=== DISTANCES MOYENNES PAR CATÉGORIE ===\n")

ancient_ids <- ids[grepl("cop", ids)]
domestic_ids <- ids[grepl("^(Awassi|Charollais|East|Hu|Kazak|Kermani|Polled|Rambouillet)$", ids)]
wild_ids <- ids[grepl("Ammon|Bighorn|Snow|Orientalis|Urial|Mouflon", ids)]

for(anc in ancient_ids) {
    cat(sprintf("\n%s:\n", anc))
    
    # Distances vers domestiques
    dist_domestic <- dist_mat[anc, domestic_ids]
    cat(sprintf("  Domestiques: min=%.4f, mean=%.4f, max=%.4f\n",
                min(dist_domestic), mean(dist_domestic), max(dist_domestic)))
    cat(sprintf("    Plus proche: %s (%.4f)\n", 
                names(which.min(dist_domestic)), min(dist_domestic)))
    
    # Distances vers sauvages
    dist_wild <- dist_mat[anc, wild_ids]
    cat(sprintf("  Sauvages: min=%.4f, mean=%.4f, max=%.4f\n",
                min(dist_wild), mean(dist_wild), max(dist_wild)))
    cat(sprintf("    Plus proche: %s (%.4f)\n", 
                names(which.min(dist_wild)), min(dist_wild)))
}

# ===== 3. ARBRE PHYLOGÉNÉTIQUE =====
cat("\nConstruction de l'arbre...\n")
dist_obj <- as.dist(dist_mat)
tree_nj <- nj(dist_obj)
tree_rooted <- midpoint(tree_nj)

write.tree(tree_nj, "tree_NJ.nwk")
write.tree(tree_rooted, "tree_NJ_rooted.nwk")

# Plot arbre avec couleurs
tip_types <- sapply(tree_rooted$tip.label, function(x) {
    if(grepl("cop", x)) return("Ancient")
    if(grepl("Ammon|Bighorn|Snow|Orientalis|Urial|Mouflon", x)) return("Wild")
    return("Domestic")
})

tip_colors <- cols[tip_types]

pdf("tree_phylogeny_colored.pdf", width=14, height=16)
plot(tree_rooted, 
     tip.color = tip_colors,
     cex = 1.2,
     font = 2,
     edge.width = 2,
     main = "Phylogénie: Anciens, Domestiques et Sauvages",
     cex.main = 1.8)
add.scale.bar(cex = 1.3, lwd = 2)
legend("topleft", 
       legend = c("Ancien", "Domestique", "Sauvage"),
       text.col = c("red", "steelblue", "darkgreen"),
       bty = "n", cex = 1.5, pch = 19, col = c("red", "steelblue", "darkgreen"))
dev.off()

# Arbre circulaire (ggtree)
library(ggtree)

tip_data <- data.frame(
    label = tree_rooted$tip.label,
    Type = tip_types,
    stringsAsFactors = FALSE
)

p_tree <- ggtree(tree_rooted, layout="circular") %<+% tip_data +
    geom_tiplab(aes(color=Type), size=3.5, fontface="bold", offset=0.002) +
    geom_tippoint(aes(color=Type), size=3) +
    scale_color_manual(values=cols) +
    theme(legend.position="right", legend.text=element_text(size=12, face="bold")) +
    labs(title="Phylogénie circulaire")

ggsave("tree_circular.pdf", p_tree, width=14, height=14)

# ===== 4. ADMIXTURE =====
cat("\nADMIXTURE plots...\n")

cv_errors <- data.frame()
for(k in 2:10) {
    log_file <- paste0("admixture_K", k, ".log")
    if(file.exists(log_file)) {
        lines <- readLines(log_file)
        cv_line <- grep("CV error", lines, value=TRUE)
        if(length(cv_line) > 0) {
            cv_val <- as.numeric(sub(".*: ", "", cv_line))
            cv_errors <- rbind(cv_errors, data.frame(K=k, CV=cv_val))
        }
    }
}

if(nrow(cv_errors) > 0) {
    p_cv <- ggplot(cv_errors, aes(K, CV)) +
        geom_line(size=1.2) +
        geom_point(size=4) +
        theme_bw(base_size=14) +
        labs(title="ADMIXTURE Cross-Validation", x="K", y="CV error")
    
    ggsave("admixture_CV.pdf", p_cv, width=10, height=6)
    
    best_k <- cv_errors$K[which.min(cv_errors$CV)]
    cat(sprintf("Meilleur K: %d\n", best_k))
    
    # Plot ADMIXTURE
    q_file <- paste0("all_samples.", best_k, ".Q")
    if(file.exists(q_file)) {
        q <- read.table(q_file)
        fam <- read.table("all_samples.fam")
        
        q$Sample <- fam$V2
        q$Type <- annot_df[q$Sample, "Type"]
        q <- q[order(match(q$Type, c("Ancient", "Domestic", "Wild"))), ]
        q$Order <- 1:nrow(q)
        
        q_melt <- melt(q, id.vars=c("Sample", "Type", "Order"))
        
        p_admix <- ggplot(q_melt, aes(x=Order, y=value, fill=variable)) +
            geom_bar(stat="identity", width=1) +
            facet_grid(~Type, scales="free_x", space="free_x") +
            theme_minimal(base_size=14) +
            labs(title=paste0("ADMIXTURE K=", best_k), x="", y="Ancestry") +
            theme(axis.text.x=element_blank(), legend.position="none",
                  strip.text=element_text(face="bold", size=12))
        
        ggsave(paste0("admixture_K", best_k, ".pdf"), p_admix, width=16, height=6)
    }
}

cat("\n✓ Toutes les visualisations générées\n")
EOF

echo ""
echo "=========================================="
echo "✓ PIPELINE TERMINÉ AVEC SUCCÈS"
echo "=========================================="
echo ""
echo "📊 RÉSULTATS:"
echo ""
echo "Visualisations:"
echo "  • PCA: ${OUT_DIR}/05_analysis/pca_PC1_PC2.pdf"
echo "  • Heatmap: ${OUT_DIR}/05_analysis/distance_heatmap_full.pdf"
echo "  • Arbre rectangulaire: ${OUT_DIR}/05_analysis/tree_phylogeny_colored.pdf"
echo "  • Arbre circulaire: ${OUT_DIR}/05_analysis/tree_circular.pdf"
echo "  • ADMIXTURE: ${OUT_DIR}/05_analysis/admixture_K*.pdf"
echo ""
echo "Données:"
echo "  • VCF final: ${OUT_DIR}/04_merged/all_samples_filtered.vcf.gz"
echo "  • PLINK: ${OUT_DIR}/05_analysis/all_samples.{bed,bim,fam}"
echo "  • Arbres Newick: ${OUT_DIR}/05_analysis/tree_NJ*.nwk"
echo ""
echo "Échantillons analysés:"
echo "  • Anciens: cop408, cop410, cop412, cop414"
echo "  • Domestiques: Awassi, Charollais, East Friesian, Hu Sheep,"
echo "                 Kazak, Kermani, Polled Dorset, Rambouillet"
echo "  • Sauvages: Mouflon, Argali (Ammon polii), Bighorn (canadensis),"
echo "              Snow Sheep (nivicola), Orientalis, Urial (vignei)"
echo ""
echo "=========================================="
