#!/usr/bin/env bash
#SBATCH --job-name=12_genus_remote_blast
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=/home/plstenge/coprolites_comparison/00_scripts/12_kraken_genus_remote_ncbi_blast_per_sample.err
#SBATCH --output=/home/plstenge/coprolites_comparison/00_scripts/12_kraken_genus_remote_ncbi_blast_per_sample.out

# ==============================================================================
# Extraction de reads Kraken2 par Cop* et genre, puis BLASTn REMOTE sur nt NCBI.
#
# Genres analyses (taxon + descendants) :
#   Capra    : 9924
#   Alces    : 9851
#   Rangifer : 9870
#
# IMPORTANT : cette version interroge les serveurs NCBI avec blastn -remote.
# Elle ne telecharge ni ne construit de base nt locale.
#
# NCBI demande de ne lancer QU'UNE recherche remote a la fois : le script est
# donc volontairement sequentiel, avec 1 CPU SLURM. Une soumission par bloc de
# sequences est effectuee pour limiter la taille des requetes. Ne pas lancer
# plusieurs jobs de ce script simultanement.
#
#  Les tables finales sont dans :
#   ${WORKDIR}/12_kraken_genus_remote_ncbi_blast_nt/summary/
# et les resultats individuels dans :
#   ${WORKDIR}/12_kraken_genus_remote_ncbi_blast_nt/<copXXX>/<Genus>/
#
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

# Resultats entierement separes de 11_kraken_genus_blast_nt.
OUT_ROOT="${WORKDIR}/12_kraken_genus_remote_ncbi_blast_nt"
TMP_ROOT="${OUT_ROOT}/tmp_cascade_fastq"
SUMMARY_DIR="${OUT_ROOT}/summary"
mkdir -p "${OUT_ROOT}" "${TMP_ROOT}" "${SUMMARY_DIR}"

LOGFILE="${WORKDIR}/00_scripts/12_kraken_genus_remote_ncbi_blast_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOGFILE}") 2>&1

# Modifier cette liste si de nouveaux Cop* doivent etre analyses.
SAMPLES=("cop408" "cop410" "cop412" "cop414" "cop417")

# TaxIDs de GENRE. --include-children inclut toutes les especes/sous-taxons.
declare -A GENERA=(
    ["Capra"]="9924"
    ["Alces"]="9851"
    ["Rangifer"]="9870"
)

# Parametres BLAST remote. Pour les reads courts/aDNA, blastn-short est adapte.
EVALUE="1e-5"
WORD_SIZE="7"
MIN_PIDENT="90"
MIN_ALN_LEN="30"
MAX_TARGET_SEQS="25"

# Nombre de sequences uniques par soumission BLAST remote. Diminuer si les
# soumissions NCBI echouent; augmenter avec prudence. Les blocs sont soumis
# sequentiellement, jamais en parallele.
REMOTE_BATCH_SIZE="100"
REMOTE_PAUSE_SECONDS="10"

for CMD in blastn taxonkit seqkit split; do
    command -v "${CMD}" >/dev/null || {
        echo "ERREUR: ${CMD} est absent de l'environnement metagenomics." >&2
        exit 1
    }
done

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

# Soumet un fichier FASTA deja decoupe. Une seule requete remote est en cours a
# tout moment; les sorties existent aussi bloc par bloc, ce qui permet de
# reprendre manuellement une analyse interrompue.
run_remote_blast_batches() {
    local fasta="$1" output_tsv="$2" batch_dir="$3"
    local batch out

    mkdir -p "${batch_dir}"
    rm -f "${batch_dir}"/*.fasta "${batch_dir}"/*.tsv

    # split -l doit respecter les enregistrements FASTA : ici 2 lignes/read,
    # car le FASTA a ete genere avec seqkit -w 0 (une sequence sur une ligne).
    split -d -a 5 -l "$((REMOTE_BATCH_SIZE * 2))" \
        --additional-suffix=.fasta "${fasta}" "${batch_dir}/batch_"

    : > "${output_tsv}"
    for batch in "${batch_dir}"/*.fasta; do
        [[ -s "${batch}" ]] || continue
        out="${batch%.fasta}.tsv"
        echo "NCBI remote blastn : $(basename "${batch}") ($(grep -c '^>' "${batch}") sequences)"

        # -remote execute la recherche sur nt chez NCBI. -num_threads n'est pas
        # utilise : il est incompatible/non pertinent pour le service distant.
        blastn \
            -remote \
            -task blastn-short \
            -query "${batch}" \
            -db nt \
            -evalue "${EVALUE}" \
            -word_size "${WORD_SIZE}" \
            -dust no \
            -soft_masking false \
            -max_target_seqs "${MAX_TARGET_SEQS}" \
            -outfmt '6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms' \
            -out "${out}"

        cat "${out}" >> "${output_tsv}"
        sleep "${REMOTE_PAUSE_SECONDS}"
    done
}

# ------------------------------------------------------------------------------
# ETAPE 1. Extraction cascade Kraken : un FASTQ par Cop* x genre.
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

    # k29 contient uniquement les U de k35; k25 les U de k29.
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

        echo "Extraction : ${SAMPLE} / ${GENUS} (TaxID ${TAXID}; descendants inclus)"
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
# ETAPE 2. BLAST remote NCBI + taxonomie, par Cop* x genre.
# ------------------------------------------------------------------------------
printf 'sample\tgenus\tspecies\treads\tpercent_within_sample_genus\n' \
    > "${SUMMARY_DIR}/blast_remote_species_counts_all_results.tsv"
printf 'sample\tgenus\traw_kraken_reads\tunique_sequences_sent_to_ncbi\tqualified_reads\n' \
    > "${SUMMARY_DIR}/blast_remote_run_summary.tsv"

for SAMPLE in "${SAMPLES[@]}"; do
    for GENUS in "${!GENERA[@]}"; do
        OUTDIR="${OUT_ROOT}/${SAMPLE}/${GENUS}"
        FQ_ALL="${OUTDIR}/${SAMPLE}_${GENUS}_all_kraken.fastq"
        BLAST_INPUT="${OUTDIR}/${SAMPLE}_${GENUS}.ncbi_remote_blast_input.fasta"
        BLAST_OUT="${OUTDIR}/${SAMPLE}_${GENUS}.ncbi_remote_blastn.tsv"
        BATCH_DIR="${OUTDIR}/ncbi_remote_batches"
        BEST_HITS="${OUTDIR}/${SAMPLE}_${GENUS}.ncbi_remote_best_hits.tsv"
        SPECIES_COUNTS="${OUTDIR}/${SAMPLE}_${GENUS}.ncbi_remote_species_counts.tsv"

        if [[ ! -s "${FQ_ALL}" ]]; then
            echo "${SAMPLE}/${GENUS}: aucun read Kraken; BLAST remote ignore."
            printf 'species\treads\tpercent\n' > "${SPECIES_COUNTS}"
            printf '%s\t%s\t0\t0\t0\n' "${SAMPLE}" "${GENUS}" \
                >> "${SUMMARY_DIR}/blast_remote_run_summary.tsv"
            continue
        fi

        # Dereplication exacte pour limiter les recherches distantes; la
        # multiplicite est conservee dans les headers (_count=N).
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

        echo "=============================================================="
        echo "BLAST NCBI REMOTE : ${SAMPLE} / ${GENUS}"
        echo "Sequences uniques : ${N_UNIQUE}; reads Kraken : ${N_RAW}"
        echo "=============================================================="

        run_remote_blast_batches "${BLAST_INPUT}" "${BLAST_OUT}" "${BATCH_DIR}"
        gzip -c "${BLAST_OUT}" > "${BLAST_OUT}.gz"

        # Filtre des alignements, meilleur hit qualifie par sequence, ajout du
        # lineage et du rang species via TaxonKit.
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
                        q=$1; count=1; \
                        if (match(q, /_count=([0-9]+)/, a)) count=a[1]; \
                        species=$NF; \
                        if (species=="" || species=="0" || species=="NA") species="Unclassified"; \
                        print q,count,$2,$3,$4,$5,$6,$7,$8,$9,$10,species \
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

        N_QUALIFIED=$(awk -F'\t' 'NR > 1 { n += $2 } END { print n+0 }' "${BEST_HITS}")
        printf '%s\t%s\t%s\t%s\t%s\n' \
            "${SAMPLE}" "${GENUS}" "${N_RAW}" "${N_UNIQUE}" "${N_QUALIFIED}" \
            >> "${SUMMARY_DIR}/blast_remote_run_summary.tsv"

        awk -F'\t' -v sample="${SAMPLE}" -v genus="${GENUS}" '\
            BEGIN { OFS="\\t" } NR > 1 { print sample, genus, $1, $2, $3 }' \
            "${SPECIES_COUNTS}" >> "${SUMMARY_DIR}/blast_remote_species_counts_all_results.tsv"
    done
done

echo "=============================================================="
echo "PIPELINE BLAST REMOTE TERMINE"
echo "Resultats individuels : ${OUT_ROOT}/<copXXX>/<Genus>/"
echo "Table extraction Kraken : ${SUMMARY_DIR}/kraken_extraction_counts.tsv"
echo "Table BLAST finale      : ${SUMMARY_DIR}/blast_remote_species_counts_all_results.tsv"
echo "Resume des soumissions  : ${SUMMARY_DIR}/blast_remote_run_summary.tsv"
echo "Log                     : ${LOGFILE}"
echo "=============================================================="

conda deactivate
