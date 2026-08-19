#!/usr/bin/env bash
#SBATCH --job-name=13_blast_kraken_fastq
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=64
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=/home/plstenge/coprolites_comparison/00_scripts/13_blast_kraken_fastq_%j.err
#SBATCH --output=/home/plstenge/coprolites_comparison/00_scripts/13_blast_kraken_fastq_%j.out

set -euo pipefail

module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics

# Ne touche pas au job remote actuel : ce script ne lit que les FASTQ deja produits.
WORKDIR="/home/plstenge/coprolites_comparison"
IN_ROOT="${WORKDIR}/12_kraken_genus_remote_ncbi_blast_nt"
OUT_ROOT="${WORKDIR}/13_blast_local_existing_kraken_fastq"

# Adapter seulement si le chemin/date de nt local est different.
BLAST_DB="/storage/biodatabanks/ncbi/NT/ncbi_blast_nt_2024-8-24/flat/nt"
THREADS="${SLURM_CPUS_PER_TASK:-64}"
EVALUE="1e-3"
MAX_TARGET_SEQS="30"

mkdir -p "${OUT_ROOT}"

for CMD in blastn seqkit; do
    command -v "${CMD}" >/dev/null 2>&1 || {
        echo "ERREUR: ${CMD} est absent de l'environnement metagenomics." >&2
        exit 1
    }
done

# IMPORTANT:
# - chaque FASTQ est converti directement en FASTA;
# - chaque read est conserve comme une sequence distincte;
# - aucune dereplication/concaténation n'est faite.
# Ainsi, un FASTA valide contient strictement deux lignes par read: >header puis sequence.

find "${IN_ROOT}" -type f -name '*_all_kraken.fastq' -size +0c | sort | while IFS= read -r FASTQ; do
    RELATIVE="${FASTQ#${IN_ROOT}/}"
    SAMPLE="$(basename "$(dirname "$(dirname "${FASTQ}")")")"
    GENUS="$(basename "$(dirname "${FASTQ}")")"

    OUTDIR="${OUT_ROOT}/${SAMPLE}/${GENUS}"
    mkdir -p "${OUTDIR}"

    FASTA="${OUTDIR}/${SAMPLE}_${GENUS}_all_kraken.fasta"
    BLAST_OUT="${OUTDIR}/${SAMPLE}_${GENUS}_all_kraken.blastn.tsv"
    DONE="${OUTDIR}/${SAMPLE}_${GENUS}_all_kraken.blastn.done"

    if [[ -s "${DONE}" ]]; then
        echo "DEJA TERMINE: ${RELATIVE}"
        continue
    fi

    echo "============================================================"
    echo "Traitement: ${SAMPLE} / ${GENUS}"
    echo "FASTQ: ${FASTQ}"

    # Conversion FASTQ -> FASTA, sequence sur une seule ligne.
    seqkit fq2fa -w 0 "${FASTQ}" -o "${FASTA}"

    N_FASTQ=$(awk 'END { print int(NR/4) }' "${FASTQ}")
    N_FASTA=$(grep -c '^>' "${FASTA}" || true)

    if [[ "${N_FASTQ}" -eq 0 || "${N_FASTA}" -eq 0 ]]; then
        echo "ATTENTION: FASTQ/FASTA vide pour ${SAMPLE}/${GENUS}; ignore."
        rm -f "${FASTA}"
        continue
    fi

    if [[ "${N_FASTQ}" -ne "${N_FASTA}" ]]; then
        echo "ERREUR: ${N_FASTQ} reads FASTQ mais ${N_FASTA} sequences FASTA pour ${SAMPLE}/${GENUS}." >&2
        exit 1
    fi

    echo "Reads envoyes au BLAST: ${N_FASTA}"

    blastn \
        -task blastn-short \
        -num_threads "${THREADS}" \
        -query "${FASTA}" \
        -db "${BLAST_DB}" \
        -evalue "${EVALUE}" \
        -max_target_seqs "${MAX_TARGET_SEQS}" \
        -outfmt '6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames' \
        -out "${BLAST_OUT}"

    touch "${DONE}"
    echo "TERMINE: ${SAMPLE}/${GENUS}"
done

echo "============================================================"
echo "BLAST local termine. Resultats: ${OUT_ROOT}"
conda deactivate
