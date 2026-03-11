#!/bin/bash

#SBATCH --job-name=Olea_vs_Fraxinus_mapping
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/olea_vs_fraxinus.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/olea_vs_fraxinus.out"

echo ""
echo "=========================================="
echo "Mapping des reads Kraken 'Olea' sur Olea_europaea et Fraxinus_excelsior"
echo "=========================================="
echo ""

# ========= ENVIRONNEMENT =========
module load conda/4.12.0
source ~/.bashrc
conda activate mapdamage_py39

BASE_DIR="/home/plstenge/coprolites_comparison"
KRAKENTOOLS_DIR="${BASE_DIR}/08_kraken2/KrakenTools"
KRAKEN_UNMERGED="${BASE_DIR}/08_kraken2/unmerged_reads"
FINAL_UNMERGED_DIR="${BASE_DIR}/06_fastp/unmerged_reads"
OUT_BASE="${BASE_DIR}/13_olea_vs_fraxinus"
THREADS=16

mkdir -p "${OUT_BASE}"
LOGFILE="${OUT_BASE}/olea_vs_fraxinus_$(date +%Y%m%d_%H%M%S).log"
touch "$LOGFILE"

# ========= GÉNOMES =========
# TaxID Olea = 4146
OLEA_TAXID="4146"

REF_OLEA="/home/plstenge/genomes/Olea_europaea/CACTIH01.fasta"
REF_FRAX="/home/plstenge/genomes/Fraxinus_excelsior/GCA_965226085.2_daFraExce3.hap1.2_genomic.fna"

for ref in "$REF_OLEA" "$REF_FRAX"; do
    if [[ ! -f "$ref" ]]; then
        echo "ERREUR: génome manquant: $ref" | tee -a "$LOGFILE"
        exit 1
    fi
    if [[ ! -f "${ref}.bwt" ]]; then
        echo "Index BWA manquant, indexation: $ref" | tee -a "$LOGFILE"
        bwa index "$ref" 2>>"$LOGFILE"
    fi
    if [[ ! -f "${ref}.fai" ]]; then
        echo "Index samtools manquant, faidx: $ref" | tee -a "$LOGFILE"
        samtools faidx "$ref" 2>>"$LOGFILE"
    fi
done

# ========= FONCTIONS =========

mapping_stats_file="${OUT_BASE}/mapping_rates_olea_vs_fraxinus.tsv"
echo -e "Sample\tRefGenome\tTotal_Reads\tMapped_Reads\tMapping_Rate(%)" > "$mapping_stats_file"

# stats globales
calculate_mapping_rate_simple() {
    local bam="$1"
    local sample="$2"
    local refname="$3"

    if [[ ! -f "$bam" ]]; then
        echo "  ⚠ BAM manquant: $bam" | tee -a "$LOGFILE"
        return
    fi

    local total mapped rate
    total=$(samtools view -c "$bam" 2>/dev/null || echo 0)
    mapped=$(samtools view -c -F 4 "$bam" 2>/dev/null || echo 0)

    if [[ "$total" -gt 0 ]]; then
        rate=$(echo "scale=3; $mapped * 100 / $total" | bc 2>/dev/null || echo "0")
    else
        rate="0"
    fi

    echo -e "${sample}\t${refname}\t${total}\t${mapped}\t${rate}" >> "$mapping_stats_file"
    echo "  Stats ${refname}: ${mapped}/${total} (${rate}%)" | tee -a "$LOGFILE"
}

# métriques par read (MAPQ, flag, insert size, etc.)
# => un TSV facilement exploitable en R
per_read_metrics() {
    local bam="$1"
    local outfile="$2"

    if [[ ! -f "$bam" ]]; then
        echo "  ⚠ BAM manquant pour per-read metrics: $bam" | tee -a "$LOGFILE"
        return
    fi

    # Colonnes SAM: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL + tags [web:6][web:9]
    # On extrait les colonnes utiles pour l’analyse
    echo -e "QNAME\tFLAG\tRNAME\tPOS\tMAPQ\tCIGAR\tRNEXT\tPNEXT\tTLEN" > "$outfile"
    samtools view "$bam" 2>>"$LOGFILE" | \
        awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' >> "$outfile"
}

# ========= BOUCLE SUR LES ÉCHANTILLONS =========

shopt -s nullglob

for KRAKEN_FILE in "${KRAKEN_UNMERGED}"/*.kraken; do
    KRAKEN_BASE=$(basename "$KRAKEN_FILE" .kraken)
    SAMPLE="${KRAKEN_BASE%_unmerged}"

    R1="${FINAL_UNMERGED_DIR}/${SAMPLE}_final_unmerged_R1.fastq.gz"
    R2="${FINAL_UNMERGED_DIR}/${SAMPLE}_final_unmerged_R2.fastq.gz"

    if [[ ! -s "$R1" || ! -s "$R2" ]]; then
        echo "⚠ FASTQ unmerged manquants ou vides pour $SAMPLE, on saute." | tee -a "$LOGFILE"
        continue
    fi

    echo ""
    echo "=========================================="
    echo "Échantillon: $SAMPLE"
    echo "==========================================" | tee -a "$LOGFILE"

    SAMPLE_OUT="${OUT_BASE}/${SAMPLE}"
    mkdir -p "$SAMPLE_OUT"

    # -------- 1) EXTRACTION READS KRAKEN ASSIGNÉS A OLEA -------- [web:1][web:8][web:11]
    OLEA_R1="${SAMPLE_OUT}/${SAMPLE}_Olea_R1.fastq"
    OLEA_R2="${SAMPLE_OUT}/${SAMPLE}_Olea_R2.fastq"

    if [[ ! -s "$OLEA_R1" || ! -s "$OLEA_R2" ]]; then
        echo "  Extraction des reads assignés à Olea (TaxID=${OLEA_TAXID})..." | tee -a "$LOGFILE"

        python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
            -k "$KRAKEN_FILE" \
            -s "$R1" -s2 "$R2" \
            -t "$OLEA_TAXID" \
            -o "$OLEA_R1" -o2 "$OLEA_R2" \
            --fastq-output 2>>"$LOGFILE"

        if [[ ! -s "$OLEA_R1" || ! -s "$OLEA_R2" ]]; then
            echo "  ⚠ Aucun read assigné à Olea pour $SAMPLE, on passe au suivant." | tee -a "$LOGFILE"
            rm -f "$OLEA_R1" "$OLEA_R2" 2>/dev/null
            continue
        fi
    else
        echo "  Reads Olea déjà extraits pour $SAMPLE, on réutilise." | tee -a "$LOGFILE"
    fi

    READS_OLEA_COUNT=$(grep -c "^@" "$OLEA_R1" 2>/dev/null || echo 0)
    echo "  ${READS_OLEA_COUNT} paires assignées à Olea extraites." | tee -a "$LOGFILE"

    # -------- 2) MAPPING SUR CHAQUE GÉNOME --------
    for REF in "$REF_OLEA" "$REF_FRAX"; do
        ref_base=$(basename "$REF")
        ref_name=""
        if [[ "$REF" == "$REF_OLEA" ]]; then
            ref_name="Olea_europaea"
        else
            ref_name="Fraxinus_excelsior"
        fi

        echo ""
        echo "  Mapping sur ${ref_name}..." | tee -a "$LOGFILE"

        MAP_OUT="${SAMPLE_OUT}/${ref_name}"
        mkdir -p "$MAP_OUT"

        SAI1="${MAP_OUT}/${SAMPLE}_${ref_name}_R1.sai"
        SAI2="${MAP_OUT}/${SAMPLE}_${ref_name}_R2.sai"
        BAM_SORTED="${MAP_OUT}/${SAMPLE}_${ref_name}.sorted.bam"

        if [[ ! -s "$BAM_SORTED" ]]; then
            # BWA aln + sampe [web:6]
            bwa aln -n 0.08 -l 24 -k 2 -q 20 -t "$THREADS" "$REF" "$OLEA_R1" \
                > "$SAI1" 2>>"$LOGFILE"
            bwa aln -n 0.08 -l 24 -k 2 -q 20 -t "$THREADS" "$REF" "$OLEA_R2" \
                > "$SAI2" 2>>"$LOGFILE"

            bwa sampe "$REF" "$SAI1" "$SAI2" "$OLEA_R1" "$OLEA_R2" 2>>"$LOGFILE" | \
                samtools view -bS - 2>>"$LOGFILE" | \
                samtools sort -@ "$THREADS" -o "$BAM_SORTED" - 2>>"$LOGFILE"

            if [[ ! -s "$BAM_SORTED" ]]; then
                echo "  ✗ Erreur lors du mapping sur ${ref_name}" | tee -a "$LOGFILE"
                continue
            fi
            samtools index "$BAM_SORTED" 2>>"$LOGFILE"
        else
            echo "  BAM déjà existant pour ${ref_name}, on réutilise." | tee -a "$LOGFILE"
        fi

        # -------- 3) STATISTIQUES GLOBALES --------
        calculate_mapping_rate_simple "$BAM_SORTED" "$SAMPLE" "$ref_name"

        # -------- 4) MÉTRIQUES PAR READ --------
        PER_READ_TSV="${MAP_OUT}/${SAMPLE}_${ref_name}_per_read_metrics.tsv"
        if [[ ! -s "$PER_READ_TSV" ]]; then
            echo "  Extraction des métriques par read pour ${ref_name}..." | tee -a "$LOGFILE"
            per_read_metrics "$BAM_SORTED" "$PER_READ_TSV"
        else
            echo "  Fichier de métriques par read déjà présent pour ${ref_name}." | tee -a "$LOGFILE"
        fi

        # (optionnel) stats samtools globales (utile aussi):
        # samtools stats "$BAM_SORTED" > "${MAP_OUT}/${SAMPLE}_${ref_name}_samtools_stats.txt" 2>>"$LOGFILE"
    done

    # Si tu veux économiser de la place, tu peux décommenter :
    # rm -f "$OLEA_R1" "$OLEA_R2"
done

shopt -u nullglob

echo ""
echo "=========================================="
echo "Terminé. Résumé: $mapping_stats_file"
echo "Log: $LOGFILE"
echo "Per-read metrics: ${OUT_BASE}/*/*_per_read_metrics.tsv"
echo "=========================================="
echo ""

conda deactivate
exit 0
