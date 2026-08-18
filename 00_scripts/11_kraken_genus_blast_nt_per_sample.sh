#!/usr/bin/env bash
#SBATCH --job-name=11_genus_blast_nt
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=36
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=/home/plstenge/coprolites_comparison/00_scripts/11_kraken_genus_blast_nt_per_sample.err
#SBATCH --output=/home/plstenge/coprolites_comparison/00_scripts/11_kraken_genus_blast_nt_per_sample.out

# ==============================================================================
# Extraction Kraken2 + validation BLASTn par ECHANTILLON et par GENRE.
#
# Resultats attendus : 3 genres x N echantillons Cop*
#   Capra    (TaxID 9924, descendants inclus)
#   Alces    (TaxID 9851, descendants inclus)
#   Rangifer (TaxID 9870, descendants inclus)
#
# Chaque combinaison Cop*/genre produit :
#   - FASTQ Kraken concatene
#   - FASTA BLAST (conserve, .fasta et .fasta.gz)
#   - BLAST brut (conserve, .tsv et .tsv.gz)
#   - meilleurs hits BLAST au rang species
#   - table de proportions des especes
#   - diagramme circulaire PDF individuel
#
# Une figure globale, avec un panneau par Cop*/genre, est egalement produite.
# ==============================================================================

module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics

WORKDIR="/home/plstenge/coprolites_comparison"
FASTP_DIR="${WORKDIR}/06_fastp"
FOLLOWING_DIR="${WORKDIR}/09_kraken_following"
K35_DIR="${FOLLOWING_DIR}/k35"
K29_DIR="${FOLLOWING_DIR}/k29"
K25_DIR="${FOLLOWING_DIR}/k25"
KRAKENTOOLS_DIR="${WORKDIR}/08_krakentools/KrakenTools"

# nt local fourni. L'index BLAST est construit une seule fois si absent.
NT_FASTA="/storage/biodatabanks/ncbi/NT/current/fasta/All/all.fasta"
BLAST_DB_PREFIX="/home/plstenge/blastdb/nt"

THREADS="${SLURM_CPUS_PER_TASK:-36}"
SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")

OUT_ROOT="${WORKDIR}/11_kraken_genus_blast_nt"
TMP_ROOT="${OUT_ROOT}/tmp_cascade_fastq"
SUMMARY_DIR="${OUT_ROOT}/summary"
mkdir -p "${OUT_ROOT}" "${TMP_ROOT}" "${SUMMARY_DIR}" "$(dirname "${BLAST_DB_PREFIX}")"

LOGFILE="${WORKDIR}/00_scripts/11_kraken_genus_blast_nt_per_sample_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOGFILE}") 2>&1

# TaxIDs de genres : --include-children recupere le genre et toutes ses especes.
declare -A GENERA=(
    ["Capra"]="9924"
    ["Alces"]="9851"
    ["Rangifer"]="9870"
)

# Seuils de conservation des hits BLAST.
EVALUE="1e-5"
WORD_SIZE="7"
MIN_PIDENT="90"
MIN_ALN_LEN="30"
MAX_TARGET_SEQS="25"

for CMD in blastn makeblastdb taxonkit seqkit Rscript; do
    command -v "${CMD}" >/dev/null || {
        echo "ERREUR: ${CMD} est introuvable dans l'environnement metagenomics." >&2
        exit 1
    }
done

[[ -s "${NT_FASTA}" ]] || {
    echo "ERREUR: FASTA nt absent/vide : ${NT_FASTA}" >&2
    exit 1
}

count_fastq_reads() {
    [[ -s "$1" ]] && awk 'END { print int(NR / 4) }' "$1" || echo 0
}

extract_unclassified() {
    local kraken="$1" fastq="$2" output_fastq="$3"
    python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
        -k "${kraken}" -s "${fastq}" -t 0 -o "${output_fastq}" --fastq-output
}

extract_genus_and_children() {
    local kraken="$1" report="$2" fastq="$3" taxid="$4" output_fastq="$5"
    python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
        -k "${kraken}" -r "${report}" -s "${fastq}" -t "${taxid}" \
        --include-children -o "${output_fastq}" --fastq-output
}

# BLAST DB v5 (.ndb) ou v4 (.nin). La construction peut etre tres longue et
# necessite beaucoup d'espace; elle est ignoree si la base existe deja.
if [[ ! -e "${BLAST_DB_PREFIX}.ndb" && ! -e "${BLAST_DB_PREFIX}.nin" ]]; then
    echo "Construction de la base BLAST nt : ${BLAST_DB_PREFIX}"
    makeblastdb \
        -in "${NT_FASTA}" \
        -dbtype nucl \
        -parse_seqids \
        -blastdb_version 5 \
        -title "NCBI_nt_current" \
        -out "${BLAST_DB_PREFIX}"
fi

# ------------------------------------------------------------------------------
# Etape 1 : extraction des reads Kraken pour chaque combinaison Cop*/genre.
# ------------------------------------------------------------------------------
printf 'sample\tgenus\treads_k35\treads_k29\treads_k25\treads_total\n' \
    > "${SUMMARY_DIR}/kraken_extraction_counts.tsv"

for SAMPLE in "${SAMPLES[@]}"; do
    ORIGINAL_FASTQ="${FASTP_DIR}/clean_${SAMPLE}_grouped_dedup_fastp_merged.fastq.gz"
    K35_KRAKEN="${K35_DIR}/${SAMPLE}_merged.kraken"
    K35_REPORT="${K35_DIR}/${SAMPLE}_merged.report"
    K29_KRAKEN="${K29_DIR}/${SAMPLE}_merged_unassigned_k35.kraken"
    K29_REPORT="${K29_DIR}/${SAMPLE}_merged_unassigned_k35.report"
    K25_KRAKEN="${K25_DIR}/${SAMPLE}_merged_unassigned_k29.kraken"
    K25_REPORT="${K25_DIR}/${SAMPLE}_merged_unassigned_k29.report"

    if [[ ! -s "${ORIGINAL_FASTQ}" || ! -s "${K35_KRAKEN}" || ! -s "${K35_REPORT}" ]]; then
        echo "ATTENTION: FASTQ/resultats k35 absents pour ${SAMPLE}; echantillon ignore."
        continue
    fi

    SAMPLE_TMP="${TMP_ROOT}/${SAMPLE}"
    mkdir -p "${SAMPLE_TMP}"
    UN_K35="${SAMPLE_TMP}/${SAMPLE}_unassigned_after_k35.fastq"
    UN_K29="${SAMPLE_TMP}/${SAMPLE}_unassigned_after_k29.fastq"

    # Reconstruction de la cascade : k29 = U de k35; k25 = U de k29.
    extract_unclassified "${K35_KRAKEN}" "${ORIGINAL_FASTQ}" "${UN_K35}"
    if [[ -s "${K29_KRAKEN}" && -s "${UN_K35}" ]]; then
        extract_unclassified "${K29_KRAKEN}" "${UN_K35}" "${UN_K29}"
    else
        : > "${UN_K29}"
    fi

    for GENUS in "${!GENERA[@]}"; do
        TAXID="${GENERA[${GENUS}]}"
        OUTDIR="${OUT_ROOT}/${SAMPLE}/${GENUS}"
        mkdir -p "${OUTDIR}"

        FQ_K35="${OUTDIR}/${SAMPLE}_${GENUS}_k35.fastq"
        FQ_K29="${OUTDIR}/${SAMPLE}_${GENUS}_k29.fastq"
        FQ_K25="${OUTDIR}/${SAMPLE}_${GENUS}_k25.fastq"
        FQ_ALL="${OUTDIR}/${SAMPLE}_${GENUS}_all_kraken.fastq"

        echo "Extraction : ${SAMPLE} / ${GENUS} (TaxID ${TAXID}, enfants inclus)"
        extract_genus_and_children "${K35_KRAKEN}" "${K35_REPORT}" "${ORIGINAL_FASTQ}" "${TAXID}" "${FQ_K35}"

        if [[ -s "${K29_KRAKEN}" && -s "${K29_REPORT}" && -s "${UN_K35}" ]]; then
            extract_genus_and_children "${K29_KRAKEN}" "${K29_REPORT}" "${UN_K35}" "${TAXID}" "${FQ_K29}"
        else
            : > "${FQ_K29}"
        fi

        if [[ -s "${K25_KRAKEN}" && -s "${K25_REPORT}" && -s "${UN_K29}" ]]; then
            extract_genus_and_children "${K25_KRAKEN}" "${K25_REPORT}" "${UN_K29}" "${TAXID}" "${FQ_K25}"
        else
            : > "${FQ_K25}"
        fi

        cat "${FQ_K35}" "${FQ_K29}" "${FQ_K25}" > "${FQ_ALL}"
        N35=$(count_fastq_reads "${FQ_K35}")
        N29=$(count_fastq_reads "${FQ_K29}")
        N25=$(count_fastq_reads "${FQ_K25}")
        NTOTAL=$(count_fastq_reads "${FQ_ALL}")
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${SAMPLE}" "${GENUS}" "${N35}" "${N29}" "${N25}" "${NTOTAL}" \
            >> "${SUMMARY_DIR}/kraken_extraction_counts.tsv"
    done

    rm -rf "${SAMPLE_TMP}"
done

# ------------------------------------------------------------------------------
# Etape 2 : BLAST et tableaux d'especes, independamment pour chaque Cop*/genre.
# ------------------------------------------------------------------------------
printf 'sample\tgenus\tspecies\treads\tpercent_within_sample_genus\n' \
    > "${SUMMARY_DIR}/blast_species_counts_all_results.tsv"

for SAMPLE in "${SAMPLES[@]}"; do
    for GENUS in "${!GENERA[@]}"; do
        OUTDIR="${OUT_ROOT}/${SAMPLE}/${GENUS}"
        FQ_ALL="${OUTDIR}/${SAMPLE}_${GENUS}_all_kraken.fastq"
        BLAST_INPUT="${OUTDIR}/${SAMPLE}_${GENUS}.blast_input.fasta"
        BLAST_OUT="${OUTDIR}/${SAMPLE}_${GENUS}.blastn.tsv"
        BEST_HITS="${OUTDIR}/${SAMPLE}_${GENUS}.best_hits.tsv"
        SPECIES_COUNTS="${OUTDIR}/${SAMPLE}_${GENUS}.blast_species_counts.tsv"

        if [[ ! -s "${FQ_ALL}" ]]; then
            echo "${SAMPLE}/${GENUS}: aucun read Kraken; pas de BLAST."
            printf 'species\treads\tpercent\n' > "${SPECIES_COUNTS}"
            continue
        fi

        # Des reads strictement identiques sont regroupes afin de reduire le
        # nombre de recherches dans nt. La multiplicite est encodee dans chaque
        # header FASTA (_count=N) puis reintroduite dans les comptages finaux.
        seqkit seq -s -w 0 "${FQ_ALL}" \
            | awk '\
                /^>/ { if (seq != "") n[seq]++; seq=""; next } \
                { seq = seq $0 } \
                END { if (seq != "") n[seq]++; for (s in n) print n[s] "\\t" s }' \
            | sort -k2,2 \
            | awk -v p="${SAMPLE}_${GENUS}" 'BEGIN { OFS="" } { print ">" p "_unique_" NR "_count=" $1 "\\n" $2 }' \
            > "${BLAST_INPUT}"
        gzip -c "${BLAST_INPUT}" > "${BLAST_INPUT}.gz"

        N_RAW=$(count_fastq_reads "${FQ_ALL}")
        N_UNIQUE=$(grep -c '^>' "${BLAST_INPUT}" || true)
        printf 'sample\tgenus\traw_kraken_reads\tunique_sequences_sent_to_blast\n%s\t%s\t%s\t%s\n' \
            "${SAMPLE}" "${GENUS}" "${N_RAW}" "${N_UNIQUE}" \
            > "${OUTDIR}/${SAMPLE}_${GENUS}.blast_input_summary.tsv"

        echo "BLASTN : ${SAMPLE} / ${GENUS} : ${N_UNIQUE} sequences uniques (${N_RAW} reads)."
        blastn \
            -task blastn-short \
            -query "${BLAST_INPUT}" \
            -db "${BLAST_DB_PREFIX}" \
            -num_threads "${THREADS}" \
            -evalue "${EVALUE}" \
            -word_size "${WORD_SIZE}" \
            -dust no \
            -soft_masking false \
            -max_target_seqs "${MAX_TARGET_SEQS}" \
            -outfmt '6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms' \
            -out "${BLAST_OUT}"
        gzip -c "${BLAST_OUT}" > "${BLAST_OUT}.gz"

        # Filtre, un meilleur hit par sequence, puis recuperation de lineage et
        # du nom au rang espece avec TaxonKit.
        {
            printf 'qseqid\tread_count\tsaccver\tpident\tlength\tevalue\tbitscore\tstaxids\tsscinames\tlineage\tranks\tspecies\n'
            if [[ -s "${BLAST_OUT}" ]]; then
                awk -F'\t' -v minpid="${MIN_PIDENT}" -v minlen="${MIN_ALN_LEN}" '\
                    ($3 >= minpid && $4 >= minlen) { print $1 "\\t" $2 "\\t" $3 "\\t" $4 "\\t" $11 "\\t" $12 "\\t" $13 "\\t" $14 }' "${BLAST_OUT}" \
                | sort -t $'\t' -k1,1 -k6,6gr -k5,5g -k4,4gr -k3,3gr \
                | awk -F'\t' '!seen[$1]++' \
                | taxonkit lineage -i 7 \
                | taxonkit reformat -f '{s}' -F -t \
                | awk -F'\t' '\
                    BEGIN { OFS="\\t" } \
                    { \
                        q=$1; c=1; \
                        if (match(q, /_count=([0-9]+)/, a)) c=a[1]; \
                        species=$NF; \
                        if (species=="" || species=="0" || species=="NA") species="Unclassified"; \
                        print q,c,$2,$3,$4,$5,$6,$7,$8,$9,$10,species \
                    }'
            fi
        } > "${BEST_HITS}"

        awk -F'\t' '\
            BEGIN { OFS="\\t" } \
            NR > 1 { sum[$12] += $2; total += $2 } \
            END { \
                print "species", "reads", "percent"; \
                for (s in sum) printf "%s\\t%d\\t%.6f\\n", s, sum[s], (total ? 100 * sum[s] / total : 0) \
            }' "${BEST_HITS}" \
            | sort -t $'\t' -k2,2nr > "${SPECIES_COUNTS}"

        awk -F'\t' -v sample="${SAMPLE}" -v genus="${GENUS}" '\
            BEGIN { OFS="\\t" } NR > 1 { print sample, genus, $1, $2, $3 }' \
            "${SPECIES_COUNTS}" >> "${SUMMARY_DIR}/blast_species_counts_all_results.tsv"
    done
done

# ------------------------------------------------------------------------------
# Etape 3 : un plot circulaire PDF par Cop*/genre + une figure globale.
# ------------------------------------------------------------------------------
RSCRIPT="${SUMMARY_DIR}/plot_blast_donuts.R"
cat > "${RSCRIPT}" <<'RSCRIPT_EOF'
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggforce)
})

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outdir <- args[2]

d <- read.delim(infile, check.names = FALSE, stringsAsFactors = FALSE)
d$reads <- as.numeric(d$reads)
d <- d[is.finite(d$reads) & d$reads > 0, , drop = FALSE]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

make_donut <- function(x, outfile, title) {
  if (nrow(x) == 0L) {
    pdf(outfile, width = 7, height = 6)
    plot.new()
    text(0.5, 0.5, "Aucun hit BLAST qualifie")
    title(main = title)
    dev.off()
    return(invisible(NULL))
  }
  x <- x[order(x$reads, decreasing = TRUE), , drop = FALSE]
  x$species_plot <- x$species
  if (nrow(x) > 12L) x$species_plot[13:nrow(x)] <- "Other species"
  x <- aggregate(reads ~ species_plot, x, sum)
  x <- x[order(x$reads, decreasing = TRUE), , drop = FALSE]
  x$end <- cumsum(x$reads)
  x$start <- x$end - x$reads
  total <- sum(x$reads)
  x$label <- ifelse(100 * x$reads / total >= 3,
                    paste0(x$species_plot, "\n", sprintf("%.1f%%", 100 * x$reads / total)), "")

  p <- ggplot(x) +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.55, r = 1,
                     start = 2 * pi * start / total,
                     end = 2 * pi * end / total,
                     fill = species_plot), colour = "white", linewidth = 0.3) +
    geom_text(aes(x = 1.27 * sin(2 * pi * (start + reads / 2) / total),
                  y = 1.27 * cos(2 * pi * (start + reads / 2) / total),
                  label = label), size = 2.8, lineheight = 0.9, check_overlap = TRUE) +
    coord_fixed(clip = "off") +
    labs(title = title,
         subtitle = paste0("Reads avec meilleur hit BLAST qualifie; n = ", total),
         fill = "Species") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9),
          legend.position = "bottom",
          legend.text = element_text(size = 7),
          plot.margin = margin(10, 35, 10, 35))
  ggsave(outfile, p, width = 9, height = 7, units = "in", limitsize = FALSE)
}

# Un PDF individuel pour chacune des 3 x N combinaisons Cop*/genre.
combos <- unique(d[, c("sample", "genus")])
for (i in seq_len(nrow(combos))) {
  s <- combos$sample[i]
  g <- combos$genus[i]
  x <- d[d$sample == s & d$genus == g, c("species", "reads"), drop = FALSE]
  make_donut(x, file.path(outdir, paste0(s, "_", g, "_blast_species_donut.pdf")),
             paste0(s, " — ", g))
}

# Figure globale : un panel par combinaison Cop*/genre.
if (nrow(d) > 0L) {
  d$panel <- paste(d$sample, d$genus, sep = " — ")
  d$species_plot <- d$species
  pdat <- do.call(rbind, lapply(split(d, d$panel), function(x) {
    x <- x[order(x$reads, decreasing = TRUE), , drop = FALSE]
    if (nrow(x) > 12L) x$species_plot[13:nrow(x)] <- "Other species"
    z <- aggregate(reads ~ sample + genus + panel + species_plot, x, sum)
    z <- z[order(z$reads, decreasing = TRUE), , drop = FALSE]
    z$end <- cumsum(z$reads)
    z$start <- z$end - z$reads
    z$total <- sum(z$reads)
    z
  }))
  p <- ggplot(pdat) +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.55, r = 1,
                     start = 2 * pi * start / total,
                     end = 2 * pi * end / total,
                     fill = species_plot), colour = "white", linewidth = 0.2) +
    coord_fixed() +
    facet_wrap(~ panel) +
    labs(title = "Validation BLASTn des assignations Kraken par echantillon et genre",
         fill = "Species") +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(size = 6))
  ggsave(file.path(outdir, "all_samples_all_genera_blast_species_donuts.pdf"),
         p, width = 15, height = 13, units = "in", limitsize = FALSE)
}
RSCRIPT_EOF

Rscript "${RSCRIPT}" \
    "${SUMMARY_DIR}/blast_species_counts_all_results.tsv" \
    "${SUMMARY_DIR}"

echo "=============================================================="
echo "PIPELINE TERMINE"
echo "Resultats individuels : ${OUT_ROOT}/<copXXX>/<Genus>/"
echo "Table Kraken globale  : ${SUMMARY_DIR}/kraken_extraction_counts.tsv"
echo "Table BLAST globale   : ${SUMMARY_DIR}/blast_species_counts_all_results.tsv"
echo "Plots individuels     : ${SUMMARY_DIR}/copXXX_<Genus>_blast_species_donut.pdf"
echo "Plot global           : ${SUMMARY_DIR}/all_samples_all_genera_blast_species_donuts.pdf"
echo "Log                   : ${LOGFILE}"
echo "=============================================================="

conda deactivate
