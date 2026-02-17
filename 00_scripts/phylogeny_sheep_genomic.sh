#!/bin/bash

#SBATCH --job-name=phylogeny_chr1_ancient_modern
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -p smp
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/phylogeny_chr1_genomes.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/phylogeny_chr1_genomes.out"

# ========== INITIALISATION ==========
echo "Initialisation..."

module load conda/4.12.0
source ~/.bashrc
conda activate genomics_analysis
which Rscript > /dev/null 2>&1 || module load R/4.3.0

for tool in bwa samtools bcftools plink; do
    which $tool > /dev/null 2>&1 || { echo "⚠️  $tool manquant"; exit 1; }
done

# ========== CONFIGURATION ==========
ANCIENT_UNMERGED_DIR="/home/plstenge/coprolites_comparison/12_mapdamage/unmerged_reads/Ovis_aries/paired_end"
ANCIENT_MERGED_DIR="/home/plstenge/coprolites_comparison/12_mapdamage/merged_reads/Ovis_aries"
MODERN_GENOMES_DIR="/home/plstenge/coprolites_comparison/16_genomes"

REF="/home/plstenge/genomes/Ovis_aries/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa"
REF_CHR1="/home/plstenge/genomes/Ovis_aries/Ovis_aries.ARS-UI_Ramb_v3.0.chr1.fa"

OUT_DIR="/home/plstenge/coprolites_comparison/14_nuclear_genome_analysis/13_phylogeny_chr1_merged_unmerged"

ANCIENT_SAMPLES=(cop408 cop410 cop412 cop414)

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

mkdir -p $OUT_DIR/{01_indexed_genomes,02_consensus,03_vcf,04_merged,05_analysis,06_bam_merged}

echo "=========================================="
echo "PHYLOGÉNIE CHR1 (ANCIENTS MERGED+UNMERGED + MODERNES)"
echo "=========================================="
echo "Anciens: ${#ANCIENT_SAMPLES[@]} coprolites"
echo "Modernes: ${#MODERN_GENOMES[@]} génomes de référence"
echo "=========================================="

# ========== ÉTAPE 0 : Merge UNMERGED + MERGED par coprolite ==========
echo ""
echo "ÉTAPE 0: Merge BAM unmerged + merged par coprolite"
echo "----------------------------------------"

for SAMPLE in "${ANCIENT_SAMPLES[@]}"; do
    UNMERGED_BAM="${ANCIENT_UNMERGED_DIR}/${SAMPLE}_unmerged_Ovis_aries.sorted.bam"
    MERGED_BAM="${ANCIENT_MERGED_DIR}/${SAMPLE}_merged_Ovis_aries_merged.sorted.bam"
    OUT_BAM="${OUT_DIR}/06_bam_merged/${SAMPLE}_merged_unmerged_Ovis_aries.sorted.bam"

    mkdir -p "${OUT_DIR}/06_bam_merged"

    if [[ -f "$OUT_BAM" ]]; then
        echo "  $SAMPLE: BAM merged déjà présent"
    else
        echo "  $SAMPLE: merge unmerged + merged..."

        if [[ ! -f "$UNMERGED_BAM" ]]; then
            echo "    ⚠️  UNMERGED manquant: $UNMERGED_BAM"
        fi
        if [[ ! -f "$MERGED_BAM" ]]; then
            echo "    ⚠️  MERGED manquant: $MERGED_BAM"
        fi

        samtools merge -@ $THREADS -f "$OUT_BAM" \
            $( [[ -f "$UNMERGED_BAM" ]] && echo "$UNMERGED_BAM" ) \
            $( [[ -f "$MERGED_BAM" ]] && echo "$MERGED_BAM" )

        samtools index "$OUT_BAM"
    fi

    N_READS=$(samtools view -c "$OUT_BAM")
    echo "    → $N_READS reads dans BAM merged"
done

# ========== ÉTAPE 1 : Préparation REF chr1 ==========
echo ""
echo "ÉTAPE 1: Préparation référence chr1"
echo "----------------------------------------"

if [[ ! -f "${REF}.bwt" ]]; then
    echo "  Indexation BWA sur REF complet..."
    bwa index "$REF"
fi

if [[ ! -f "${REF}.fai" ]]; then
    echo "  Indexation FASTA sur REF complet..."
    samtools faidx "$REF"
fi

if [[ ! -f "$REF_CHR1" ]]; then
    echo "  Extraction du chr1 de la référence..."
    samtools faidx "$REF" 1 > "$REF_CHR1"
fi

if [[ ! -f "${REF_CHR1}.bwt" ]]; then
    echo "  Indexation BWA sur REF chr1..."
    bwa index "$REF_CHR1"
fi

if [[ ! -f "${REF_CHR1}.fai" ]]; then
    echo "  Indexation FASTA sur REF chr1..."
    samtools faidx "$REF_CHR1"
fi

# ========== ÉTAPE 2 : Consensus modernes CHR1 ==========
echo ""
echo "ÉTAPE 2: Génération des consensus modernes (CHR1)"
echo "----------------------------------------"
echo "  Stratégie: extraction du contig correspondant à chr1"

for NAME in "${!MODERN_GENOMES[@]}"; do
    GENOME_PATH="${MODERN_GENOMES_DIR}/${MODERN_GENOMES[$NAME]}"
    CONSENSUS="${OUT_DIR}/02_consensus/${NAME}_chr1_consensus.fa.gz"

    if [[ -f "$CONSENSUS" ]]; then
        echo "  $NAME: consensus déjà présent"
        continue
    fi

    echo "  $NAME: extraction contig chr1..."

    zcat "$GENOME_PATH" | awk '
        BEGIN { keep=0; n=0 }
        /^>/ {
            if ($0 ~ /^>chr1(\s|$)/ || $0 ~ /^>1(\s|$)/) {
                print $0; keep=1; next
            }
            if (keep==1) { exit }
            headers[++n]=$0
            next
        }
        {
            if (keep==1) { print; next }
            seq[n]=seq[n] $0
        }
        END {
            if (keep==0 && n>0) {
                print headers[1]
                print seq[1]
            }
        }
    ' | gzip > "$CONSENSUS"

    echo "    → Consensus CHR1 créé: $CONSENSUS"
done

# ========== ÉTAPE 3 : SNP calling échantillons anciens (FORCÉ) ==========
echo ""
echo "ÉTAPE 3: SNP calling échantillons anciens (CHR1, merged+unmerged)"
echo "----------------------------------------"
echo "  ⚠️  Recalcul forcé : suppression des VCF anciens"

mkdir -p "${OUT_DIR}/03_vcf"

for SAMPLE in "${ANCIENT_SAMPLES[@]}"; do
    VCF="${OUT_DIR}/03_vcf/${SAMPLE}.vcf.gz"
    VCF_INDEX="${VCF}.csi"
    if [[ -f "$VCF" ]]; then
        echo "  Suppression ancien VCF: $VCF"
        rm -f "$VCF" "$VCF_INDEX"
    fi
done

for SAMPLE in "${ANCIENT_SAMPLES[@]}"; do
    BAM="${OUT_DIR}/06_bam_merged/${SAMPLE}_merged_unmerged_Ovis_aries.sorted.bam"
    VCF="${OUT_DIR}/03_vcf/${SAMPLE}.vcf.gz"

    if [[ ! -f "$BAM" ]]; then
        echo "  ⚠️  BAM merged manquant pour $SAMPLE, skip"
        continue
    fi

    echo "  $SAMPLE..."

    bcftools mpileup \
        -f "$REF_CHR1" \
        -Q 20 -q 20 -A \
        --max-depth 200000 \
        -r 1 \
        "$BAM" | \
    bcftools call -m -v --ploidy 2 | \
    bcftools reheader -s <(echo "$SAMPLE") | \
    bcftools view -Oz -o "$VCF"

    bcftools index -f "$VCF"

    N=$(bcftools view -H "$VCF" | wc -l)
    echo "    → $N SNPs (chromosome 1)"
done

# ========== ÉTAPE 4 : SNPs sur consensus modernes (CHR1) ==========
echo ""
echo "ÉTAPE 4: Extraction variants depuis consensus modernes (CHR1)"
echo "----------------------------------------"

for NAME in "${!MODERN_GENOMES[@]}"; do
    CONSENSUS="${OUT_DIR}/02_consensus/${NAME}_chr1_consensus.fa.gz"
    VCF="${OUT_DIR}/03_vcf/${NAME}.vcf.gz"

    if [[ -f "$VCF" ]]; then
        echo "  $NAME: VCF déjà présent, skip"
        continue
    fi

    if [[ ! -f "$CONSENSUS" ]]; then
        echo "  ⚠️  Consensus manquant pour $NAME, skip"
        continue
    fi

    echo "  $NAME: mapping + SNP calling (CHR1)..."

    TMP_BAM="${OUT_DIR}/03_vcf/${NAME}_tmp.bam"
    TMP_LOG="${OUT_DIR}/03_vcf/${NAME}_bwa.log"

    zcat "$CONSENSUS" | \
    bwa mem -t $THREADS "$REF_CHR1" - 2> "$TMP_LOG" | \
    samtools sort -@ $THREADS -o "$TMP_BAM"

    samtools index "$TMP_BAM"

    bcftools mpileup -f "$REF_CHR1" -r 1 --max-depth 500 "$TMP_BAM" | \
    bcftools call -m -v --ploidy 2 | \
    bcftools reheader -s <(echo "$NAME") | \
    bcftools view -Oz -o "$VCF"

    bcftools index -f "$VCF"

    rm -f "$TMP_BAM" "$TMP_BAM.bai"

    N=$(bcftools view -H "$VCF" | wc -l)
    echo "    → $N SNPs (CHR1)"
done

# ========== ÉTAPE 5 : Merge VCFs ==========
echo ""
echo "ÉTAPE 5: Merge de tous les échantillons"
echo "----------------------------------------"

ALL_VCFS=""

for SAMPLE in "${ANCIENT_SAMPLES[@]}"; do
    VCF="${OUT_DIR}/03_vcf/${SAMPLE}.vcf.gz"
    if [[ -f "$VCF" ]]; then
        ALL_VCFS="$ALL_VCFS $VCF"
    else
        echo "  ⚠️  VCF manquant pour $SAMPLE, non inclus"
    fi
done

for NAME in "${!MODERN_GENOMES[@]}"; do
    VCF="${OUT_DIR}/03_vcf/${NAME}.vcf.gz"
    if [[ -f "$VCF" ]]; then
        ALL_VCFS="$ALL_VCFS $VCF"
    else
        echo "  ⚠️  VCF moderne manquant pour $NAME, non inclus"
    fi
done

N_VCFS=$(echo $ALL_VCFS | wc -w)
echo "  Merge de $N_VCFS VCFs..."

if [[ $N_VCFS -lt 2 ]]; then
    echo "  ⚠️  Moins de 2 VCFs, arrêt."
    exit 1
fi

mkdir -p "${OUT_DIR}/04_merged"

bcftools merge $ALL_VCFS -Oz -o ${OUT_DIR}/04_merged/all_samples_raw.vcf.gz
bcftools index -f ${OUT_DIR}/04_merged/all_samples_raw.vcf.gz

N_RAW=$(bcftools view -H ${OUT_DIR}/04_merged/all_samples_raw.vcf.gz | wc -l)
echo "  → $N_RAW SNPs bruts (chr1)"

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
    echo "  ⚠️  Peu de SNPs, relaxation des filtres (F_MISSING<0.7)..."

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

mkdir -p ${OUT_DIR}/05_analysis

plink --vcf ${OUT_DIR}/04_merged/all_samples_filtered.vcf.gz \
      --sheep \
      --allow-extra-chr \
      --make-bed \
      --set-missing-var-ids @:# \
      --out ${OUT_DIR}/05_analysis/all_samples

N_SNPS_PLINK=$(wc -l < ${OUT_DIR}/05_analysis/all_samples.bim)
echo "  ✓ $N_SNPS_PLINK SNPs dans PLINK"

# ========== ÉTAPE 7 : PCA + distances + ADMIXTURE + R ==========
echo ""
echo "ÉTAPE 7: Analyses PLINK"
echo "----------------------------------------"

cd ${OUT_DIR}/05_analysis

plink --bfile all_samples \
      --pca 10 \
      --allow-extra-chr \
      --out pca

plink --bfile all_samples \
      --distance square 1-ibs \
      --allow-extra-chr \
      --out distances

if command -v admixture >/dev/null 2>&1; then
    for K in {2..10}; do
        echo "  ADMIXTURE K=$K..."
        admixture --cv all_samples.bed $K -j4 > admixture_K${K}.log 2>&1
    done
fi

Rscript <<'EOF'
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(ape)
library(phangorn)
library(RColorBrewer)
library(reshape2)
library(ggtree)

setwd(".")

if(file.exists("pca.eigenvec")){
  pca <- read.table("pca.eigenvec", header=FALSE)
  colnames(pca) <- c("FID", "IID", paste0("PC", 1:(ncol(pca)-2)))

  pca$Category <- "Modern_Domestic"
  pca$Category[grepl("cop", pca$IID, ignore.case=TRUE)] <- "Ancient"
  pca$Category[grepl("Ammon|Bighorn|Snow|Orientalis|Urial|Mouflon", pca$IID)] <- "Wild"
  pca$Type <- factor(pca$Category, levels=c("Ancient","Modern_Domestic","Wild"))

  cols <- c("Ancient"="red","Modern_Domestic"="steelblue","Wild"="darkgreen")

  if(ncol(pca) >= 4){
    p1 <- ggplot(pca, aes(PC1, PC2, color=Type, label=IID)) +
      geom_point(size=4, alpha=0.8) +
      geom_text_repel(size=3.5, max.overlaps=30) +
      scale_color_manual(values=cols) +
      theme_bw(base_size=14) +
      labs(title="PCA chr1: anciens + modernes",
           x="PC1", y="PC2") +
      theme(legend.position="top")
    ggsave("pca_PC1_PC2.pdf", p1, width=10, height=8)
  }
}

if(file.exists("distances.dist")){
  dist_mat <- as.matrix(read.table("distances.dist"))
  ids <- read.table("distances.dist.id")$V1
  rownames(dist_mat) <- colnames(dist_mat) <- ids

  annot_df <- data.frame(
    Type = ifelse(grepl("cop", ids, ignore.case=TRUE), "Ancient",
           ifelse(grepl("Ammon|Bighorn|Snow|Orientalis|Urial|Mouflon", ids),
                  "Wild","Domestic")),
    row.names = ids
  )
  annot_colors <- list(Type=c(Ancient="red",Domestic="steelblue",Wild="darkgreen"))

  pdf("distance_heatmap_full.pdf", width=11, height=10)
  pheatmap(dist_mat,
           annotation_row=annot_df,
           annotation_col=annot_df,
           annotation_colors=annot_colors,
           main="Distances génétiques chr1 (1-IBS)")
  dev.off()

  dist_obj <- as.dist(dist_mat)
  tree_nj <- nj(dist_obj)
  tree_rooted <- midpoint(tree_nj)
  write.tree(tree_rooted, "tree_NJ_rooted.nwk")

  tip_types <- sapply(tree_rooted$tip.label, function(x){
    if(grepl("cop", x, ignore.case=TRUE)) return("Ancient")
    if(grepl("Ammon|Bighorn|Snow|Orientalis|Urial|Mouflon", x)) return("Wild")
    return("Domestic")
  })
  cols <- c(Ancient="red",Domestic="steelblue",Wild="darkgreen")

  pdf("tree_phylogeny_colored.pdf", width=12, height=10)
  plot(tree_rooted, tip.color=cols[tip_types], cex=1.0,
       main="Phylogénie chr1: anciens + modernes")
  legend("topleft", legend=names(cols), col=cols, pch=19, bty="n")
  dev.off()

  p_tree <- ggtree(tree_rooted, layout="circular") %<+%
    data.frame(label=tree_rooted$tip.label,Type=tip_types) +
    geom_tiplab(aes(color=Type), size=3) +
    scale_color_manual(values=cols)
  ggsave("tree_circular.pdf", p_tree, width=10, height=10)
}
EOF

echo ""
echo "=========================================="
echo "✓ PIPELINE CHR1 TERMINÉ"
echo "=========================================="
echo "Résultats dans: ${OUT_DIR}/05_analysis"
