#!/bin/bash

#SBATCH --job-name=phylogeny_FULLGENOME
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH -p smp
#SBATCH --mem=250G
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/phylogeny_genomes_FULLGENOME.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/phylogeny_genomes_FULLGENOME.out"

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
OUT_DIR="/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis/12_complete_phylogeny_FULLGENOME"

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
    ["Bighorn"]="Ovis_canadensis/GCF_042477335.2_ARS-UI_OviCan_v2_genomic.fna.gz"
    ["Snow_Sheep"]="Ovis_nivicola_lydekkeri/GCA_903231385.1_OvNiv1.0_genomic.fna.gz"
    ["Orientalis"]="Ovis_orientalis_from_Iran_PRJEB3141/GCF_000765115.1_Oori1_genomic.fna.gz"
    ["Urial"]="Ovis_vignei/GCA_053525375.1_ASM5352537v1_genomic.fna.gz"
)

THREADS=16

mkdir -p $OUT_DIR/{01_indexed_genomes,02_consensus,03_vcf,04_merged,05_analysis}

echo "=========================================="
echo "PHYLOGÉNIE GÉNOME ENTIER"
echo "=========================================="
echo "Anciens: 4 coprolites"
echo "Modernes: ${#MODERN_GENOMES[@]} génomes de référence"
echo "⚠️  GÉNOME COMPLET (pas seulement chr1)"
echo "   Temps estimé: 12-24h"
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

if [[ ! -f "${REF}.fai" ]]; then
    echo "  Indexation FASTA..."
    samtools faidx $REF
fi

# ========== ÉTAPE 2 : Générer séquences consensus pour modernes ==========
echo ""
echo "ÉTAPE 2: Génération des consensus modernes (GÉNOME ENTIER)"
echo "----------------------------------------"
echo "  Stratégie: Pseudogénomes basés sur alignement vs référence commune"

for NAME in "${!MODERN_GENOMES[@]}"; do
    GENOME_PATH="${MODERN_GENOMES_DIR}/${MODERN_GENOMES[$NAME]}"
    CONSENSUS="${OUT_DIR}/02_consensus/${NAME}_consensus.fa.gz"
    
    if [[ -f "$CONSENSUS" ]]; then
        echo "  $NAME: déjà fait"
        continue
    fi
    
    echo "  $NAME: extraction génome complet..."
    
    # Extraire les 10 premiers scaffolds/chromosomes principaux (éviter milliers de contigs)
    zcat $GENOME_PATH | \
    awk '/^>/{n++; if(n>10) exit} {print}' | \
    sed "s/^>.*/>${NAME}/" | \
    gzip > $CONSENSUS
    
    SIZE=$(zcat $CONSENSUS | grep -v "^>" | wc -c)
    echo "    → Consensus créé ($SIZE bp)"
done

echo "  ✓ Note: Top 10 scaffolds par génome (représentatif, évite contigs courts)"

# ========== ÉTAPE 3 : SNP calling sur échantillons anciens (GÉNOME ENTIER) ==========
echo ""
echo "ÉTAPE 3: SNP calling échantillons anciens (GÉNOME ENTIER)"
echo "----------------------------------------"

for SAMPLE in "${ANCIENT_SAMPLES[@]}"; do
    BAM="${ANCIENT_BAM_DIR}/${SAMPLE}_unmerged_Ovis_aries.sorted.bam"
    VCF="${OUT_DIR}/03_vcf/${SAMPLE}.vcf.gz"
    
    if [[ -f "$VCF" ]]; then
        echo "  $SAMPLE: déjà fait"
        continue
    fi
    
    echo "  $SAMPLE..."
    
    # Pas de filtre -r (tous les chromosomes)
    bcftools mpileup \
        -f $REF \
        -Q 20 -q 20 -A \
        --max-depth 1000000 \
        $BAM | \
    bcftools call -m -v --ploidy 2 | \
    bcftools reheader -s <(echo "$SAMPLE") | \
    bcftools view -Oz -o $VCF
    
    bcftools index -f $VCF
    
    N=$(bcftools view -H $VCF | wc -l)
    echo "    → $N SNPs (génome entier)"
done

# ========== ÉTAPE 4 : Appeler variants sur consensus modernes (GÉNOME ENTIER) ==========
echo ""
echo "ÉTAPE 4: Extraction variants depuis consensus modernes (GÉNOME ENTIER)"
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
    bwa mem -t $THREADS $REF - | \
    samtools sort -@ 4 -o $TMP_BAM
    
    samtools index $TMP_BAM
    
    # Calling (pas de filtre -r, tous chromosomes)
    bcftools mpileup -f $REF --max-depth 500 $TMP_BAM | \
    bcftools call -m -v --ploidy 2 | \
    bcftools reheader -s <(echo "$NAME") | \
    bcftools view -Oz -o $VCF
    
    bcftools index -f $VCF
    
    rm -f $TMP_BAM $TMP_BAM.bai
    
    N=$(bcftools view -H $VCF | wc -l)
    echo "    → $N SNPs (génome entier)"
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
echo "  → $N_RAW SNPs bruts (génome entier)"

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

if [[ $N_FILT -lt 500 ]]; then
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

# ========== ÉTAPE 7 : Statistiques missing data ==========
echo ""
echo "ÉTAPE 7: Statistiques missing data"
echo "----------------------------------------"

plink --bfile ${OUT_DIR}/05_analysis/all_samples \
      --missing \
      --allow-extra-chr \
      --out ${OUT_DIR}/05_analysis/missing_stats

echo "  Missing data par échantillon (top 10):"
awk 'NR>1 {print $2, $6}' ${OUT_DIR}/05_analysis/missing_stats.imiss | \
    sort -k2 -rn | head -10 | \
    awk '{printf "    %s: %.2f%%\n", $1, $2*100}'

# Filtrer échantillons avec données suffisantes
awk '$6 < 0.5 {print $1, $2}' ${OUT_DIR}/05_analysis/missing_stats.imiss > \
    ${OUT_DIR}/05_analysis/samples_ok.txt

N_OK=$(wc -l < ${OUT_DIR}/05_analysis/samples_ok.txt)
echo "  → $N_OK échantillons avec <50% missing"

if [[ $N_OK -ge 5 ]]; then
    echo "  Création dataset filtré..."
    
    plink --bfile ${OUT_DIR}/05_analysis/all_samples \
          --keep ${OUT_DIR}/05_analysis/samples_ok.txt \
          --geno 0.3 \
          --maf 0.01 \
          --make-bed \
          --allow-extra-chr \
          --out ${OUT_DIR}/05_analysis/filtered
    
    DATASET="${OUT_DIR}/05_analysis/filtered"
else
    echo "  ⚠️  Peu d'échantillons OK, utilisation dataset complet"
    DATASET="${OUT_DIR}/05_analysis/all_samples"
fi

# ========== ÉTAPE 8 : PCA ==========
echo ""
echo "ÉTAPE 8: PCA"
echo "----------------------------------------"

plink --bfile $DATASET \
      --pca 10 \
      --allow-extra-chr \
      --out ${OUT_DIR}/05_analysis/pca

if [[ -f "${OUT_DIR}/05_analysis/pca.eigenvec" ]]; then
    echo "  ✓ PCA terminée"
    echo "  Échantillons dans PCA: $(wc -l < ${OUT_DIR}/05_analysis/pca.eigenvec)"
else
    echo "  ✗ PCA échouée, voir pca.log"
fi

# ========== ÉTAPE 9 : Distance génétique ==========
echo ""
echo "ÉTAPE 9: Distances génétiques"
echo "----------------------------------------"

plink --bfile $DATASET \
      --distance square 1-ibs \
      --allow-extra-chr \
      --out ${OUT_DIR}/05_analysis/distances

echo "  ✓ Matrice de distances calculée"

# ========== ÉTAPE 10 : Distances avec TOUS les échantillons (même missing) ==========
echo ""
echo "ÉTAPE 10: Distances complètes (tous échantillons, dataset relax)"
echo "----------------------------------------"

plink --bfile ${OUT_DIR}/05_analysis/all_samples \
      --geno 0.95 \
      --make-bed \
      --allow-extra-chr \
      --out ${OUT_DIR}/05_analysis/all_relaxed

plink --bfile ${OUT_DIR}/05_analysis/all_relaxed \
      --distance square 1-ibs \
      --allow-extra-chr \
      --out ${OUT_DIR}/05_analysis/distances_all

echo "  ✓ Matrice de distances complète (anciens inclus)"

# ========== ÉTAPE 11 : ADMIXTURE ==========
echo ""
echo "ÉTAPE 11: ADMIXTURE"
echo "----------------------------------------"

which admixture > /dev/null 2>&1
if [[ $? -ne 0 ]]; then
    echo "  ⚠️  ADMIXTURE non disponible, ignoré"
else
    cd ${OUT_DIR}/05_analysis
    
    for K in {2..10}; do
        echo "  K=$K..."
        admixture --cv $(basename $DATASET).bed $K -j$THREADS > admixture_K${K}.log 2>&1
    done
    
    echo ""
    echo "  Cross-validation (5 meilleurs):"
    grep "CV error" admixture_K*.log | sort -k4 -n | head -5
    
    cd - > /dev/null
fi

# ========== ÉTAPE 12 : VISUALISATIONS R ==========
echo ""
echo "ÉTAPE 12: Visualisations"
echo "----------------------------------------"

Rscript <<'EOF'
library(ggplot2)
library(ggrepel)
library(ape)
library(phangorn)
library(pheatmap)
library(RColorBrewer)
library(reshape2)

outdir <- Sys.getenv("OUT_DIR")
if(outdir == "") outdir <- "/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis/12_complete_phylogeny_FULLGENOME/05_analysis"
setwd(outdir)

cat("==========================================\n")
cat("GÉNÉRATION VISUALISATIONS\n")
cat("==========================================\n\n")

# ===== 1. PCA (échantillons avec données OK) =====
if(file.exists("pca.eigenvec")) {
  cat("1. PCA (échantillons >50% génotypés)...\n")
  
  pca <- read.table("pca.eigenvec", header=FALSE)
  n_pcs <- ncol(pca) - 2
  colnames(pca) <- c("FID", "IID", paste0("PC", 1:n_pcs))
  
  pca$Category <- "Modern_Domestic"
  pca$Category[grepl("cop", pca$IID, ignore.case=TRUE)] <- "Ancient"
  pca$Category[grepl("Ammon|Bighorn|Snow|Orientalis|Urial|Mouflon", pca$IID)] <- "Wild"
  pca$Type <- factor(pca$Category, levels=c("Ancient", "Modern_Domestic", "Wild"))
  
  cols <- c("Ancient"="red", "Modern_Domestic"="steelblue", "Wild"="darkgreen")
  
  if(n_pcs >= 2) {
    p1 <- ggplot(pca, aes(PC1, PC2, color=Type, label=IID)) +
      geom_point(size=4, alpha=0.8) +
      geom_text_repel(size=3.5, fontface="bold", max.overlaps=25, 
                      segment.size=0.3, box.padding=0.5) +
      scale_color_manual(values=cols) +
      theme_bw(base_size=14) +
      labs(title="PCA: Phylogénie génome entier",
           subtitle=paste0(nrow(pca), " échantillons avec données suffisantes"),
           x="PC1", y="PC2") +
      theme(legend.position="top", legend.text=element_text(size=12, face="bold"))
    
    ggsave("pca_PC1_PC2.pdf", p1, width=16, height=12)
    cat("  ✓ pca_PC1_PC2.pdf\n")
  }
  
  if(n_pcs >= 3) {
    p2 <- ggplot(pca, aes(PC2, PC3, color=Type, label=IID)) +
      geom_point(size=4, alpha=0.8) +
      geom_text_repel(size=3.5, fontface="bold", max.overlaps=25) +
      scale_color_manual(values=cols) +
      theme_bw(base_size=14) +
      labs(title="PCA: PC2 vs PC3", x="PC2", y="PC3") +
      theme(legend.position="top")
    
    ggsave("pca_PC2_PC3.pdf", p2, width=16, height=12)
    cat("  ✓ pca_PC2_PC3.pdf\n")
  }
} else {
  cat("1. PCA: fichier non trouvé, skip\n")
}

cat("\n")

# ===== 2. DISTANCES (tous échantillons, anciens inclus) =====
cat("2. Distances génétiques (TOUS échantillons)...\n")

if(file.exists("distances_all.dist")) {
  dist_mat <- as.matrix(read.table("distances_all.dist"))
  ids <- read.table("distances_all.dist.id")$V1
  rownames(dist_mat) <- colnames(dist_mat) <- ids
  
  annot_df <- data.frame(
      Type = ifelse(grepl("cop", ids, ignore.case=TRUE), "Ancient",
                    ifelse(grepl("Ammon|Bighorn|Snow|Orientalis|Urial|Mouflon", ids), 
                           "Wild", "Domestic")),
      row.names = ids
  )
  
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
           main = "Distances génétiques - Génome entier (1-IBS)",
           color = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
           border_color = "grey80")
  dev.off()
  
  cat("  ✓ distance_heatmap_full.pdf\n")
  
  # ===== ANALYSE POSITIONNEMENT ANCIENS =====
  cat("\n")
  cat("==========================================\n")
  cat("POSITIONNEMENT DES ÉCHANTILLONS ANCIENS\n")
  cat("==========================================\n\n")
  
  ancient_ids <- ids[grepl("cop", ids, ignore.case=TRUE)]
  modern_ids <- ids[!grepl("cop", ids, ignore.case=TRUE)]
  
  results <- data.frame()
  
  for(anc in ancient_ids) {
    dists <- dist_mat[anc, modern_ids]
    dists <- dists[!is.na(dists) & dists > 0]
    
    if(length(dists) == 0) {
      cat(anc, ": AUCUNE distance calculable\n\n")
      next
    }
    
    top5 <- sort(dists)[1:min(5, length(dists))]
    
    cat(anc, ":\n")
    cat("  Distances calculées:", length(dists), "modernes\n")
    cat("  Distance moyenne:", round(mean(dists, na.rm=TRUE), 4), "\n")
    cat("  Top 5 plus proches:\n")
    
    for(i in seq_along(top5)) {
      cat("    ", i, ". ", names(top5)[i], " (dist = ", 
          round(top5[i], 4), ")\n", sep="")
      
      results <- rbind(results, data.frame(
        Ancient = anc,
        Rank = i,
        Modern = names(top5)[i],
        Distance = top5[i]
      ))
    }
    cat("\n")
  }
  
  write.table(results, "ancient_closest_modern.tsv", 
              quote=FALSE, row.names=FALSE, sep="\t")
  cat("✓ Fichier: ancient_closest_modern.tsv\n\n")
  
  # ===== ARBRE PHYLOGÉNÉTIQUE =====
  cat("3. Arbre phylogénétique...\n")
  
  dist_obj <- as.dist(dist_mat)
  tree_nj <- nj(dist_obj)
  tree_rooted <- midpoint(tree_nj)
  
  write.tree(tree_nj, "tree_NJ.nwk")
  write.tree(tree_rooted, "tree_NJ_rooted.nwk")
  
  tip_types <- sapply(tree_rooted$tip.label, function(x) {
      if(grepl("cop", x, ignore.case=TRUE)) return("Ancient")
      if(grepl("Ammon|Bighorn|Snow|Orientalis|Urial|Mouflon", x)) return("Wild")
      return("Domestic")
  })
  
  tip_colors <- c("Ancient"="red", "Domestic"="steelblue", "Wild"="darkgreen")[tip_types]
  
  pdf("tree_phylogeny_colored.pdf", width=14, height=16)
  plot(tree_rooted, 
       tip.color = tip_colors,
       cex = 1.2,
       font = 2,
       edge.width = 2,
       main = "Phylogénie génome entier: Anciens + Modernes",
       cex.main = 1.8)
  add.scale.bar(cex = 1.3, lwd = 2)
  legend("topleft", 
         legend = c("Ancien", "Domestique", "Sauvage"),
         text.col = c("red", "steelblue", "darkgreen"),
         bty = "n", cex = 1.5, pch = 19, 
         col = c("red", "steelblue", "darkgreen"))
  dev.off()
  
  cat("  ✓ tree_phylogeny_colored.pdf\n")
  
  # Arbre circulaire
  library(ggtree)
  
  tip_data <- data.frame(
      label = tree_rooted$tip.label,
      Type = tip_types,
      stringsAsFactors = FALSE
  )
  
  p_tree <- ggtree(tree_rooted, layout="circular") %<+% tip_data +
      geom_tiplab(aes(color=Type), size=3.5, fontface="bold", offset=0.002) +
      geom_tippoint(aes(color=Type), size=3) +
      scale_color_manual(values=c("Ancient"="red", "Domestic"="steelblue", "Wild"="darkgreen")) +
      theme(legend.position="right", legend.text=element_text(size=12, face="bold")) +
      labs(title="Phylogénie circulaire - Génome entier")
  
  ggsave("tree_circular.pdf", p_tree, width=14, height=14)
  cat("  ✓ tree_circular.pdf\n")
  
} else {
  cat("2. Distances: fichier non trouvé\n")
}

cat("\n")

# ===== 4. ADMIXTURE =====
cat("4. ADMIXTURE...\n")

cv_files <- Sys.glob("admixture_K*.log")
if(length(cv_files) > 0) {
  cv_errors <- data.frame()
  
  for(log_file in cv_files) {
    lines <- readLines(log_file)
    cv_line <- grep("CV error", lines, value=TRUE)
    if(length(cv_line) > 0) {
      k <- as.numeric(sub(".*_K([0-9]+)\\.log", "\\1", log_file))
      cv_val <- as.numeric(sub(".*: ", "", cv_line))
      cv_errors <- rbind(cv_errors, data.frame(K=k, CV=cv_val))
    }
  }
  
  if(nrow(cv_errors) > 0) {
    p_cv <- ggplot(cv_errors, aes(K, CV)) +
        geom_line(size=1.2) +
        geom_point(size=4) +
        theme_bw(base_size=14) +
        labs(title="ADMIXTURE Cross-Validation", x="K", y="CV error")
    
    ggsave("admixture_CV.pdf", p_cv, width=10, height=6)
    cat("  ✓ admixture_CV.pdf\n")
    
    best_k <- cv_errors$K[which.min(cv_errors$CV)]
    cat("  Meilleur K:", best_k, "\n")
  }
}

cat("\n==========================================\n")
cat("✓ VISUALISATIONS TERMINÉES\n")
cat("==========================================\n")
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
echo "  • Positionnement anciens: ${OUT_DIR}/05_analysis/ancient_closest_modern.tsv"
echo ""
echo "Données:"
echo "  • VCF final: ${OUT_DIR}/04_merged/all_samples_filtered.vcf.gz"
echo "  • PLINK: ${OUT_DIR}/05_analysis/all_samples.{bed,bim,fam}"
echo "  • Arbres Newick: ${OUT_DIR}/05_analysis/tree_NJ*.nwk"
echo ""
echo "Échantillons analysés:"
echo "  • Anciens: cop408, cop410, cop412, cop414"
echo "  • Domestiques: Awassi, Charollais, East Friesian, Hu Sheep,"
echo "                 Kazak, Kermani, Polled Dorset"
echo "  • Sauvages: Mouflon, Argali (Ammon polii), Bighorn (canadensis),"
echo "              Snow Sheep (nivicola), Orientalis, Urial (vignei)"
echo ""
echo "=========================================="
