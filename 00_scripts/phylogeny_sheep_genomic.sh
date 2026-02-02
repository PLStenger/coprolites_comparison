#!/bin/bash

#SBATCH --job-name=ancient_modern_imputation
#SBATCH --ntasks=8
#SBATCH -p smp
#SBATCH --mem=200G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/phylogeny_imputation.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/phylogeny_imputation.out"

# ========== INITIALISATION ==========
echo "Initialisation..."
module load conda/4.12.0
source ~/.bashrc
conda activate genomics_analysis

# Charger R si nécessaire
which Rscript > /dev/null 2>&1 || module load R/4.3.0

# Vérifier BEAGLE
if [[ ! -f "/home/plstenge/softwares/beagle.jar" ]]; then
    echo "Téléchargement de BEAGLE..."
    mkdir -p /home/plstenge/softwares
    cd /home/plstenge/softwares
    wget https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar
    ln -sf beagle.22Jul22.46e.jar beagle.jar
    cd -
fi

BEAGLE="/home/plstenge/softwares/beagle.jar"

# ========== CONFIG ==========
BAM_DIR="/home/plstenge/coprolites_comparison/12_mapdamage/unmerged_reads/Ovis_aries/paired_end"
REF="/home/plstenge/genomes/Ovis_aries/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa"
OUT_DIR="/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis"
MODERN_VCF="/home/plstenge/coprolites_comparison/15_modern_sheep_vcf/ISGC_SNP50_Breedv2.vcf.gz"

SAMPLES=(cop408 cop410 cop412 cop414)
THREADS=8

mkdir -p $OUT_DIR/{01_snps,02_merge,03_imputation,04_analysis,05_logs}

echo "=========================================="
echo "PIPELINE IMPUTATION + PHYLOGÉNIE"
echo "=========================================="
echo "Stratégie: Imputer les génotypes anciens"
echo "          sur les positions de l'array moderne"
echo "=========================================="

# ========== ÉTAPE 1 : SNP calling ADN ancien ==========
echo ""
echo "ÉTAPE 1: SNP calling échantillons anciens"
echo "----------------------------------------"

for SAMPLE in "${SAMPLES[@]}"; do
    BAM_FILE="${BAM_DIR}/${SAMPLE}_unmerged_Ovis_aries.sorted.bam"
    VCF_OUT="${OUT_DIR}/01_snps/${SAMPLE}.vcf.gz"
    
    if [[ -f "$VCF_OUT" ]]; then
        echo "  $SAMPLE: déjà fait"
        continue
    fi
    
    echo "  Processing $SAMPLE..."
    
    bcftools mpileup \
        -f $REF \
        -Q 20 -q 20 -A \
        --max-depth 1000000 \
        -a FORMAT/AD,FORMAT/DP,FORMAT/GP \
        $BAM_FILE | \
    bcftools call \
        -m -v \
        --ploidy 2 \
        -Oz -o $VCF_OUT
    
    bcftools index -f $VCF_OUT
    
    N=$(bcftools view -H $VCF_OUT | wc -l)
    echo "    → $N SNPs détectés"
done

# ========== ÉTAPE 2 : Préparer le VCF moderne ==========
echo ""
echo "ÉTAPE 2: Préparation VCF moderne (panel de référence)"
echo "----------------------------------------"

# Vérifier si déjà compressé avec bgzip
if [[ ! -f "${OUT_DIR}/02_merge/modern_ref.vcf.gz" ]]; then
    echo "  Recompression avec bgzip..."
    
    # Décompresser si gzip standard, puis recompresser avec bgzip
    if file $MODERN_VCF | grep -q "gzip compressed"; then
        gunzip -c $MODERN_VCF | bgzip -c > ${OUT_DIR}/02_merge/modern_ref.vcf.gz
    else
        bgzip -c $MODERN_VCF > ${OUT_DIR}/02_merge/modern_ref.vcf.gz
    fi
    
    bcftools index -f ${OUT_DIR}/02_merge/modern_ref.vcf.gz
else
    echo "  ✓ Panel moderne déjà préparé"
fi

MODERN_REF="${OUT_DIR}/02_merge/modern_ref.vcf.gz"
N_MODERN=$(bcftools view -H $MODERN_REF | wc -l)
echo "  → $N_MODERN positions dans l'array moderne"

# ========== ÉTAPE 3 : Merger échantillons anciens ==========
echo ""
echo "ÉTAPE 3: Merge des échantillons anciens"
echo "----------------------------------------"

VCF_LIST=""
for SAMPLE in "${SAMPLES[@]}"; do
    VCF_LIST="$VCF_LIST ${OUT_DIR}/01_snps/${SAMPLE}.vcf.gz"
done

bcftools merge \
    --force-samples \
    $VCF_LIST \
    -Oz -o ${OUT_DIR}/02_merge/ancient_merged.vcf.gz

bcftools index -f ${OUT_DIR}/02_merge/ancient_merged.vcf.gz

N_ANCIENT=$(bcftools view -H ${OUT_DIR}/02_merge/ancient_merged.vcf.gz | wc -l)
echo "  → $N_ANCIENT SNPs anciens"

# Filtrer
bcftools view \
    -i 'QUAL>20 && FORMAT/DP>2' \
    -v snps \
    ${OUT_DIR}/02_merge/ancient_merged.vcf.gz \
    -Oz -o ${OUT_DIR}/02_merge/ancient_filtered.vcf.gz

bcftools index -f ${OUT_DIR}/02_merge/ancient_filtered.vcf.gz

N_FILT=$(bcftools view -H ${OUT_DIR}/02_merge/ancient_filtered.vcf.gz | wc -l)
echo "  ✓ $N_FILT SNPs anciens filtrés"

# ========== ÉTAPE 4 : Créer VCF des positions modernes pour les anciens ==========
echo ""
echo "ÉTAPE 4: Extraction des positions modernes dans les BAMs anciens"
echo "----------------------------------------"

echo "  Création du fichier de positions..."
bcftools query -f '%CHROM\t%POS\n' $MODERN_REF > ${OUT_DIR}/03_imputation/modern_positions.txt

N_POS=$(wc -l < ${OUT_DIR}/03_imputation/modern_positions.txt)
echo "  → $N_POS positions à extraire"

# Appeler génotypes sur TOUTES les positions de l'array (même sans couverture)
echo "  Calling sur positions modernes (mode ALL sites)..."

for SAMPLE in "${SAMPLES[@]}"; do
    BAM_FILE="${BAM_DIR}/${SAMPLE}_unmerged_Ovis_aries.sorted.bam"
    VCF_OUT="${OUT_DIR}/03_imputation/${SAMPLE}_at_modern_pos.vcf.gz"
    
    echo "    $SAMPLE..."
    
    # Utiliser -T pour cibler uniquement les positions de l'array
    bcftools mpileup \
        -f $REF \
        -T ${OUT_DIR}/03_imputation/modern_positions.txt \
        -Q 20 -q 20 -A \
        --max-depth 1000000 \
        -a FORMAT/AD,FORMAT/DP \
        $BAM_FILE | \
    bcftools call \
        -m \
        --ploidy 2 \
        -Oz -o $VCF_OUT
    
    bcftools index -f $VCF_OUT
    
    # Stats couverture
    N_COVERED=$(bcftools query -f '[%DP\n]' $VCF_OUT | awk '$1>0' | wc -l)
    echo "      → $N_COVERED positions avec couverture (sur $N_POS)"
done

# Merger les anciens sur positions modernes
echo "  Merge des anciens..."

VCF_LIST_MOD=""
for SAMPLE in "${SAMPLES[@]}"; do
    VCF_LIST_MOD="$VCF_LIST_MOD ${OUT_DIR}/03_imputation/${SAMPLE}_at_modern_pos.vcf.gz"
done

bcftools merge \
    --force-samples \
    --missing-to-ref \
    $VCF_LIST_MOD \
    -Oz -o ${OUT_DIR}/03_imputation/ancient_at_modern_pos.vcf.gz

bcftools index -f ${OUT_DIR}/03_imputation/ancient_at_modern_pos.vcf.gz

# ========== ÉTAPE 5 : IMPUTATION avec BEAGLE ==========
echo ""
echo "ÉTAPE 5: Imputation avec BEAGLE"
echo "----------------------------------------"

echo "  Lancement de BEAGLE (peut prendre 10-30 min)..."
echo "  Panel de référence: moutons modernes"
echo "  Cibles: échantillons anciens"

java -Xmx180g -jar $BEAGLE \
    gt=${OUT_DIR}/03_imputation/ancient_at_modern_pos.vcf.gz \
    ref=$MODERN_REF \
    out=${OUT_DIR}/03_imputation/ancient_imputed \
    nthreads=$THREADS \
    impute=true \
    gp=true \
    ne=10000

if [[ $? -ne 0 ]]; then
    echo "  ✗ BEAGLE a échoué"
    echo ""
    echo "  Solutions de repli:"
    echo "  1. Vérifier compatibilité des chromosomes entre ancien/moderne"
    echo "  2. Réduire le nombre de positions"
    echo "  3. Utiliser uniquement les positions avec couverture"
    exit 1
fi

# Renommer sortie
mv ${OUT_DIR}/03_imputation/ancient_imputed.vcf.gz ${OUT_DIR}/03_imputation/ancient_imputed_raw.vcf.gz
bcftools index -f ${OUT_DIR}/03_imputation/ancient_imputed_raw.vcf.gz

echo "  ✓ Imputation terminée"

# Filtrer génotypes imputés de faible qualité (GP < 0.8)
echo "  Filtrage des génotypes de faible qualité..."

bcftools +setGT ${OUT_DIR}/03_imputation/ancient_imputed_raw.vcf.gz \
    -- -t q -i 'FMT/GP<0.8' -n . | \
bcftools view -Oz -o ${OUT_DIR}/03_imputation/ancient_imputed_filtered.vcf.gz

bcftools index -f ${OUT_DIR}/03_imputation/ancient_imputed_filtered.vcf.gz

# ========== ÉTAPE 6 : Merge ancien imputé + moderne ==========
echo ""
echo "ÉTAPE 6: Merge final ancien imputé + moderne"
echo "----------------------------------------"

bcftools merge \
    --force-samples \
    ${OUT_DIR}/03_imputation/ancient_imputed_filtered.vcf.gz \
    $MODERN_REF \
    -Oz -o ${OUT_DIR}/04_analysis/combined_imputed.vcf.gz

bcftools index -f ${OUT_DIR}/04_analysis/combined_imputed.vcf.gz

N_COMBINED=$(bcftools view -H ${OUT_DIR}/04_analysis/combined_imputed.vcf.gz | wc -l)
echo "  ✓ $N_COMBINED SNPs dans le dataset combiné"

# ========== ÉTAPE 7 : Filtres finaux pour phylogénie ==========
echo ""
echo "ÉTAPE 7: Filtrage final"
echo "----------------------------------------"

# Retirer SNPs avec trop de missing data
bcftools view \
    -i 'F_MISSING<0.2' \
    ${OUT_DIR}/04_analysis/combined_imputed.vcf.gz \
    -Oz -o ${OUT_DIR}/04_analysis/combined_clean.vcf.gz

bcftools index -f ${OUT_DIR}/04_analysis/combined_clean.vcf.gz

N_CLEAN=$(bcftools view -H ${OUT_DIR}/04_analysis/combined_clean.vcf.gz | wc -l)
echo "  ✓ $N_CLEAN SNPs après filtrage (missing<20%)"

# Conversion PLINK
echo "  Conversion PLINK..."
plink --vcf ${OUT_DIR}/04_analysis/combined_clean.vcf.gz \
      --sheep \
      --allow-extra-chr \
      --make-bed \
      --set-missing-var-ids @:# \
      --mind 0.5 \
      --geno 0.2 \
      --maf 0.01 \
      --out ${OUT_DIR}/04_analysis/combined

N_FINAL=$(wc -l < ${OUT_DIR}/04_analysis/combined.bim)
echo "  → $N_FINAL SNPs dans dataset final"

# ========== ÉTAPE 8 : PCA ==========
echo ""
echo "ÉTAPE 8: PCA"
echo "----------------------------------------"

plink --bfile ${OUT_DIR}/04_analysis/combined \
      --pca 10 \
      --allow-extra-chr \
      --out ${OUT_DIR}/04_analysis/pca

Rscript <<'EOF'
library(ggplot2)
library(ggrepel)

# Lire PCA
pca <- read.table("/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis/04_analysis/pca.eigenvec", 
                  header=FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:10))

# Identifier anciens vs modernes
pca$Type <- ifelse(grepl("cop", pca$IID), "Ancient", "Modern")
pca$Label <- ifelse(pca$Type == "Ancient", as.character(pca$IID), "")

# Plot PC1 vs PC2
p1 <- ggplot(pca, aes(PC1, PC2, color=Type)) +
  geom_point(size=3, alpha=0.7) +
  geom_text_repel(aes(label=Label), size=4, fontface="bold", 
                  max.overlaps=20, color="red") +
  scale_color_manual(values=c("Ancient"="red", "Modern"="steelblue")) +
  theme_bw(base_size=14) +
  labs(title="PCA: Moutons anciens (imputés) vs modernes",
       subtitle="Basé sur SNP50 array après imputation",
       x=paste0("PC1 (", round(100*0.15, 1), "%)"),  # Ajuster avec variance expliquée
       y=paste0("PC2 (", round(100*0.10, 1), "%)")) +
  theme(legend.position="top")

ggsave("/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis/04_analysis/pca_PC1_PC2.pdf",
       p1, width=12, height=10)
ggsave("/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis/04_analysis/pca_PC1_PC2.png",
       p1, width=12, height=10, dpi=300)

# Plot PC2 vs PC3
p2 <- ggplot(pca, aes(PC2, PC3, color=Type)) +
  geom_point(size=3, alpha=0.7) +
  geom_text_repel(aes(label=Label), size=4, fontface="bold",
                  max.overlaps=20, color="red") +
  scale_color_manual(values=c("Ancient"="red", "Modern"="steelblue")) +
  theme_bw(base_size=14) +
  labs(title="PCA: Moutons anciens vs modernes (PC2-PC3)",
       x="PC2", y="PC3") +
  theme(legend.position="top")

ggsave("/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis/04_analysis/pca_PC2_PC3.pdf",
       p2, width=12, height=10)

cat("✓ PCA plots sauvegardés\n")
EOF

echo "  ✓ PCA terminée"

# ========== ÉTAPE 9 : ADMIXTURE ==========
echo ""
echo "ÉTAPE 9: ADMIXTURE"
echo "----------------------------------------"

which admixture > /dev/null 2>&1
if [[ $? -ne 0 ]]; then
    echo "  ⚠️  ADMIXTURE non disponible, installation..."
    conda install -c bioconda admixture -y
fi

cd ${OUT_DIR}/04_analysis

# Tester K de 2 à 8
for K in {2..8}; do
    echo "  Running K=$K..."
    admixture --cv combined.bed $K -j$THREADS > admixture_K${K}.log 2>&1
done

echo ""
echo "  Cross-validation errors:"
grep "CV error" admixture_K*.log | sort -k4 -n | head -5

# Plot ADMIXTURE
Rscript <<'EOF'
library(ggplot2)
library(reshape2)

setwd("/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis/04_analysis")

# Trouver meilleur K
cv_errors <- data.frame()
for(k in 2:8) {
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

# Plot CV error
p_cv <- ggplot(cv_errors, aes(K, CV)) +
    geom_line(size=1) +
    geom_point(size=3) +
    theme_bw(base_size=14) +
    labs(title="ADMIXTURE Cross-Validation",
         x="K (nombre de populations)",
         y="Cross-validation error") +
    scale_x_continuous(breaks=2:8)

ggsave("admixture_CV_plot.pdf", p_cv, width=8, height=6)

# Identifier meilleur K
best_k <- cv_errors$K[which.min(cv_errors$CV)]
cat(paste0("\n✓ Meilleur K: ", best_k, "\n"))

# Plot ADMIXTURE pour meilleur K
q_file <- paste0("combined.", best_k, ".Q")
if(file.exists(q_file)) {
    q <- read.table(q_file)
    fam <- read.table("combined.fam")
    
    q$Sample <- fam$V2
    q$Type <- ifelse(grepl("cop", q$Sample), "Ancient", "Modern")
    q$Order <- 1:nrow(q)
    
    # Trier: anciens en premier
    q <- q[order(q$Type, decreasing=TRUE), ]
    q$Order <- 1:nrow(q)
    
    # Reshape pour ggplot
    q_melt <- melt(q, id.vars=c("Sample", "Type", "Order"),
                   variable.name="Cluster", value.name="Proportion")
    
    p_admix <- ggplot(q_melt, aes(x=Order, y=Proportion, fill=Cluster)) +
        geom_bar(stat="identity", width=1) +
        theme_minimal(base_size=12) +
        labs(title=paste0("ADMIXTURE K=", best_k),
             x="Échantillons", y="Ancestry proportion") +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              panel.grid=element_blank()) +
        facet_grid(~Type, scales="free_x", space="free_x")
    
    ggsave(paste0("admixture_K", best_k, "_plot.pdf"), p_admix, width=14, height=6)
    cat("✓ ADMIXTURE plot sauvegardé\n")
}
EOF

cd - > /dev/null

# ========== ÉTAPE 10 : Distance génétique et arbre ==========
echo ""
echo "ÉTAPE 10: Distance génétique et phylogénie"
echo "----------------------------------------"

plink --bfile ${OUT_DIR}/04_analysis/combined \
      --distance square \
      --allow-extra-chr \
      --out ${OUT_DIR}/04_analysis/distances

Rscript <<'EOF'
library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)

# Lire distance
dist_file <- "/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis/04_analysis/distances.dist"
dist_id <- "/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis/04_analysis/distances.dist.id"

ids <- read.table(dist_id, header=FALSE)$V1
dist_mat <- as.matrix(read.table(dist_file, header=FALSE))
rownames(dist_mat) <- ids
colnames(dist_mat) <- ids

dist_obj <- as.dist(dist_mat)

# Arbre NJ
tree_nj <- nj(dist_obj)

# Identifier anciens
ancient_samples <- grep("cop", tree_nj$tip.label, value=TRUE)
tree_nj$tip.label.type <- ifelse(tree_nj$tip.label %in% ancient_samples, "Ancient", "Modern")

# Sauvegarder
write.tree(tree_nj, "/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis/04_analysis/tree_NJ.nwk")

# Plot avec ggtree
pdf("/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis/04_analysis/tree_NJ_colored.pdf",
    width=12, height=14)

tip_colors <- ifelse(grepl("cop", tree_nj$tip.label), "red", "steelblue")
tip_sizes <- ifelse(grepl("cop", tree_nj$tip.label), 3, 1.5)

ggtree(tree_nj, layout="rectangular") +
    geom_tiplab(aes(label=label), size=2.5, 
                color=tip_colors,
                fontface=ifelse(grepl("cop", tree_nj$tip.label), "bold", "plain")) +
    geom_tippoint(color=tip_colors, size=tip_sizes) +
    ggtitle("Arbre phylogénétique: Moutons anciens vs modernes") +
    theme_tree2()

dev.off()

cat("✓ Arbre phylogénétique sauvegardé\n")
EOF

# ========== STATISTIQUES FINALES ==========
echo ""
echo "ÉTAPE 11: Statistiques d'imputation"
echo "----------------------------------------"

Rscript <<'EOF'
# Évaluer qualité de l'imputation
vcf_imp <- "/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis/03_imputation/ancient_imputed_raw.vcf.gz"

# Extraire GP scores
system(paste0("bcftools query -f '[%SAMPLE\\t%GP\\n]' ", vcf_imp, " -s cop408,cop410,cop412,cop414 > /tmp/gp_scores.txt"))

gp <- read.table("/tmp/gp_scores.txt", header=FALSE, sep="\t")
colnames(gp) <- c("Sample", "GP")

# Parser GP (format: prob(0/0),prob(0/1),prob(1/1))
gp$GP_max <- sapply(strsplit(as.character(gp$GP), ","), function(x) {
    probs <- as.numeric(x)
    if(all(is.na(probs))) return(NA)
    max(probs, na.rm=TRUE)
})

# Résumé par échantillon
cat("\nQualité de l'imputation (probabilité génotype max):\n")
for(s in c("cop408", "cop410", "cop412", "cop414")) {
    gp_sample <- gp$GP_max[gp$Sample == s]
    gp_sample <- gp_sample[!is.na(gp_sample)]
    
    cat(sprintf("  %s:\n", s))
    cat(sprintf("    - Médiane GP: %.3f\n", median(gp_sample)))
    cat(sprintf("    - GP > 0.9: %.1f%%\n", 100*mean(gp_sample > 0.9)))
    cat(sprintf("    - GP > 0.8: %.1f%%\n", 100*mean(gp_sample > 0.8)))
    cat(sprintf("    - GP < 0.5: %.1f%%\n", 100*mean(gp_sample < 0.5)))
}

# Plot distribution
library(ggplot2)

gp_clean <- gp[!is.na(gp$GP_max), ]

p <- ggplot(gp_clean, aes(x=GP_max, fill=Sample)) +
    geom_histogram(bins=50, alpha=0.7) +
    facet_wrap(~Sample, ncol=1) +
    theme_bw(base_size=14) +
    labs(title="Distribution des probabilités de génotypes imputés",
         x="Probabilité max du génotype",
         y="Nombre de SNPs") +
    geom_vline(xintercept=0.8, linetype="dashed", color="red") +
    theme(legend.position="none")

ggsave("/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis/05_logs/imputation_quality.pdf",
       p, width=10, height=10)

cat("\n✓ Statistiques sauvegardées\n")
EOF

# ========== RÉSUMÉ FINAL ==========
echo ""
echo "=========================================="
echo "✓ PIPELINE TERMINÉ AVEC SUCCÈS"
echo "=========================================="
echo ""
echo "Résultats principaux:"
echo "  • PCA: ${OUT_DIR}/04_analysis/pca_PC1_PC2.pdf"
echo "  • ADMIXTURE: ${OUT_DIR}/04_analysis/admixture_K*_plot.pdf"
echo "  • Arbre: ${OUT_DIR}/04_analysis/tree_NJ_colored.pdf"
echo "  • Newick: ${OUT_DIR}/04_analysis/tree_NJ.nwk"
echo "  • Qualité imputation: ${OUT_DIR}/05_logs/imputation_quality.pdf"
echo ""
echo "Fichiers VCF:"
echo "  • Anciens bruts: ${OUT_DIR}/02_merge/ancient_filtered.vcf.gz"
echo "  • Anciens imputés: ${OUT_DIR}/03_imputation/ancient_imputed_filtered.vcf.gz"
echo "  • Combiné final: ${OUT_DIR}/04_analysis/combined_clean.vcf.gz"
echo ""
echo "Statistiques:"
echo "  • SNPs anciens détectés: $N_FILT"
echo "  • Positions array moderne: $N_MODERN"
echo "  • SNPs finaux après imputation: $N_FINAL"
echo ""
echo "=========================================="
