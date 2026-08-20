#!/usr/bin/env bash
#SBATCH --job-name=15_run_competitive_mapping_lca
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=64
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=/home/plstenge/coprolites_comparison/00_scripts/15_run_competitive_mapping_lca.err
#SBATCH --output=/home/plstenge/coprolites_comparison/00_scripts/15_run_competitive_mapping_lca.out

module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics

########################################
# CONFIGURATION
########################################

BASE=/home/plstenge/coprolites_comparison
GENOMES=/home/plstenge/genomes
OUTDIR=${BASE}/14_competitive_mapping_lca
THREADS=16

SAMPLES=(cop408 cop410 cop412 cop414 cop417)

mkdir -p ${OUTDIR}/{ref,mapped,dedup,filtered,lca,coverage,logs}

########################################
# ETAPE 1 : construire la reference competitive
########################################
# On liste ici les fichiers genomiques a utiliser pour chaque espece.
# IMPORTANT : Ovis_aries dispose d'un mitogenome dedie -> on l'ajoute aussi,
# ca aide beaucoup pour l'authentification aDNA (haute couverture mito).

REF_LIST=(
  "Alces_alces:${GENOMES}/Alces_alces/GCA_059051365.1_mAlcAlc2_p1.1_genomic.fna"
  "Capra_hircus:${GENOMES}/Capra_hircus/Capra_hircus.ARS1.dna.toplevel.fixed.fa"
  "Capra_ibex:${GENOMES}/Capra_ibex/GCA_054642885.1_CapIbe1.0_genomic.fna"
  "Capra_sibirica:${GENOMES}/Capra_sibirica/GCA_003182615.2_ASM318261v2_genomic.fna"
  "Capra_aegagrus:${GENOMES}/Capra_aegagrus/GCA_000978405.1_CapAeg_1.0_genomic.fna"
  "Ovis_aries:${GENOMES}/Ovis_aries/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fixed.fa"
  "Rangifer_tarandus:${GENOMES}/Rangifer_tarandus/GCA_949782905.1_mRanTar1.h1.1_genomic.fna"
  "Mus_musculus:${GENOMES}/Mus_musculus/Mus_musculus.GRCm39.dna.toplevel.fa"
  "Bos_taurus:${GENOMES}/Bos_taurus/GCF_002263795.3_ARS-UCD2.0_genomic.fna"
)

# Optionnel: mitogenome Ovis_aries dedie, si present
if [ -d "${GENOMES}/Ovis_aries/Ovis_aries_mitochondrion" ]; then
  MITO_FILE=$(find "${GENOMES}/Ovis_aries/Ovis_aries_mitochondrion" -name "*.fa*" | head -n1)
  if [ -n "${MITO_FILE:-}" ]; then
    REF_LIST+=("Ovis_aries_mito:${MITO_FILE}")
  fi
fi

COMPETITIVE_REF=${OUTDIR}/ref/competitive_reference.fasta
ACC2TAX=${OUTDIR}/ref/acc2taxid.custom.tsv
TAXID_LIST=${OUTDIR}/ref/taxid_species.tsv

# Table taxid NCBI officiels pour chaque espece cible (fixes, verifies sur NCBI Taxonomy)
declare -A TAXIDS=(
  [Alces_alces]=9852
  [Capra_hircus]=9925
  [Capra_ibex]=90542
  [Capra_sibirica]=97835
  [Capra_aegagrus]=9838
  [Ovis_aries]=9940
  [Ovis_aries_mito]=9940
  [Rangifer_tarandus]=9870
  [Mus_musculus]=10090
  [Bos_taurus]=9913
)

echo -e "species\ttaxid" > ${TAXID_LIST}
for k in "${!TAXIDS[@]}"; do echo -e "${k}\t${TAXIDS[$k]}"; done >> ${TAXID_LIST}

echo "[1/9] Construction de la reference competitive avec headers prefixes par espece..."
> ${COMPETITIVE_REF}
> ${ACC2TAX}

for entry in "${REF_LIST[@]}"; do
  species="${entry%%:*}"
  fasta="${entry#*:}"
  taxid="${TAXIDS[$species]}"

  if [ ! -f "$fasta" ]; then
    echo "  ATTENTION: fichier manquant pour ${species}: ${fasta}" >&2
    continue
  fi

  echo "  -> ajout de ${species} (taxid ${taxid}) depuis $(basename ${fasta})"

  # Renommer chaque header pour inclure l'espece, garder l'accession d'origine
  awk -v sp="${species}" '
    /^>/ {
      acc=substr($1,2)
      print ">"sp"__"acc" "substr($0, index($0,$2))
      next
    }
    { print }
  ' "$fasta" >> ${COMPETITIVE_REF}

  # Construire acc2taxid : chaque header devient sp__acc -> taxid
  grep "^>" "$fasta" | awk -v sp="${species}" -v tid="${taxid}" '
    { acc=substr($1,2); print sp"__"acc"\t"tid }
  ' >> ${ACC2TAX}
done

echo "  Reference competitive: ${COMPETITIVE_REF}"
echo "  Table acc2taxid: ${ACC2TAX}"

########################################
# ETAPE 2 : indexer la reference avec bwa + samtools faidx
########################################
echo "[2/9] Indexation bwa de la reference competitive (peut prendre du temps)..."
if [ ! -f "${COMPETITIVE_REF}.bwt" ]; then
  bwa index ${COMPETITIVE_REF} 2> ${OUTDIR}/logs/bwa_index.log
fi
samtools faidx ${COMPETITIVE_REF}

########################################
# ETAPE 3 a 7 : par echantillon : bwa aln -> sam -> bam -> tri -> dedup -> filtre MAPQ
########################################
for sample in "${SAMPLES[@]}"; do
  FASTQ=${BASE}/06_fastp/clean_${sample}_grouped_dedup_fastp_merged.fastq.gz

  if [ ! -f "$FASTQ" ]; then
    echo "ATTENTION: fastq manquant pour ${sample}: ${FASTQ}" >&2
    continue
  fi

  echo "[3/9] Mapping bwa aln (parametres aDNA) pour ${sample}..."
  SAI=${OUTDIR}/mapped/${sample}.sai
  bwa aln -t ${THREADS} -n 0.01 -o 2 -l 1024 \
    ${COMPETITIVE_REF} ${FASTQ} > ${SAI} 2> ${OUTDIR}/logs/${sample}_bwa_aln.log

  echo "[4/9] bwa samse -> BAM trie pour ${sample}..."
  BAM_SORTED=${OUTDIR}/mapped/${sample}.sorted.bam
  bwa samse -r "@RG\tID:${sample}\tSM:${sample}\tPL:illumina" \
    ${COMPETITIVE_REF} ${SAI} ${FASTQ} 2> ${OUTDIR}/logs/${sample}_samse.log \
    | samtools view -b -F 4 -@ ${THREADS} - \
    | samtools sort -@ ${THREADS} -o ${BAM_SORTED} -
  samtools index ${BAM_SORTED}

  echo "[5/9] Deduplication (samtools markdup, single-end) pour ${sample}..."
  BAM_FIXMATE=${OUTDIR}/dedup/${sample}.fixmate.bam
  BAM_DEDUP=${OUTDIR}/dedup/${sample}.dedup.bam
  samtools sort -n -@ ${THREADS} -o ${OUTDIR}/dedup/${sample}.nsorted.bam ${BAM_SORTED}
  samtools fixmate -m ${OUTDIR}/dedup/${sample}.nsorted.bam ${BAM_FIXMATE}
  samtools sort -@ ${THREADS} -o ${OUTDIR}/dedup/${sample}.possorted.bam ${BAM_FIXMATE}
  samtools markdup -r -@ ${THREADS} ${OUTDIR}/dedup/${sample}.possorted.bam ${BAM_DEDUP}
  samtools index ${BAM_DEDUP}
  rm -f ${OUTDIR}/dedup/${sample}.nsorted.bam ${OUTDIR}/dedup/${sample}.possorted.bam

  echo "[6/9] Filtrage MAPQ (>=25) et edit distance pour ${sample}..."
  BAM_FILTERED=${OUTDIR}/filtered/${sample}.filtered.bam
  samtools view -h -q 25 -@ ${THREADS} ${BAM_DEDUP} \
    | awk '/^@/ || ($0 ~ /NM:i:[0-3]([^0-9]|$)/)' \
    | samtools view -b -@ ${THREADS} -o ${BAM_FILTERED} -
  samtools index ${BAM_FILTERED}

  echo "[7/9] Statistiques rapides pour ${sample}..."
  samtools flagstat ${BAM_FILTERED} > ${OUTDIR}/logs/${sample}_flagstat.txt
  samtools idxstats ${BAM_FILTERED} > ${OUTDIR}/coverage/${sample}_idxstats.tsv
done

########################################
# ETAPE 8 : ngsLCA / metaDMG-cpp lca
########################################
echo "[8/9] Assignation LCA avec metaDMG-cpp lca..."

for sample in "${SAMPLES[@]}"; do
  BAM_FILTERED=${OUTDIR}/filtered/${sample}.filtered.bam
  [ -f "$BAM_FILTERED" ] || continue

  metaDMG-cpp lca \
    -bam ${BAM_FILTERED} \
    -names ~/ncbi_tax_dmp/names.dmp \
    -nodes ~/ncbi_tax_dmp/nodes.dmp \
    -acc2tax ${ACC2TAX} \
    -edit_dist_min 0 -edit_dist_max 5 \
    -simscorelow 0.95 -simscorehigh 1.0 \
    -out ${OUTDIR}/lca/${sample} \
    > ${OUTDIR}/logs/${sample}_lca.log 2>&1

  echo "  -> resultats LCA pour ${sample} : ${OUTDIR}/lca/${sample}.lca"
done

########################################
# ETAPE 9 : verification distribution des reads le long du genome
# (evite un signal artefactuel du a une seule region ultra-conservee)
########################################
echo "[9/9] Calcul de la couverture par fenetre pour verifier la distribution des reads..."

for sample in "${SAMPLES[@]}"; do
  BAM_FILTERED=${OUTDIR}/filtered/${sample}.filtered.bam
  [ -f "$BAM_FILTERED" ] || continue

  BEDCOV=${OUTDIR}/coverage/${sample}_windows_1kb.bed
  awk 'BEGIN{OFS="\t"} {print $1, $2}' ${COMPETITIVE_REF}.fai > ${OUTDIR}/coverage/genome.sizes

  samtools bedcov <(awk 'BEGIN{OFS="\t"} {for(i=0;i<$2;i+=1000) print $1,i,(i+1000>$2?$2:i+1000)}' ${OUTDIR}/coverage/genome.sizes) \
    ${BAM_FILTERED} > ${BEDCOV} 2> ${OUTDIR}/logs/${sample}_bedcov.log || true

  echo "  -> couverture par fenetres 1kb pour ${sample} : ${BEDCOV}"
done

echo ""
echo "=== PIPELINE TERMINE ==="
echo "Reference competitive : ${COMPETITIVE_REF}"
echo "Table acc2taxid       : ${ACC2TAX}"
echo "BAMs filtres          : ${OUTDIR}/filtered/"
echo "Resultats LCA         : ${OUTDIR}/lca/"
echo "Couverture par fenetre: ${OUTDIR}/coverage/"
