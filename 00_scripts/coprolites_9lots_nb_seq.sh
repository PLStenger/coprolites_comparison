#!/bin/bash

#SBATCH --job-name=coprolites_9lots_nb_seq
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/coprolites_9lots_nb_seq.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/coprolites_9lots_nb_seq.out"

################################################################################
# Script de récapitulation des statistiques MultiQC pour 9 lots
# Author: Pierre-Louis Stenger
# Date: December 2025
#
# Ce script parcourt les 9 lots et calcule :
# - Somme des Total Sequences
# - Moyenne des Total Sequences
# - Écart-type des Total Sequences
# - Moyenne de avg_sequence_length
# - Moyenne de median_sequence_length
################################################################################

set -eo pipefail

# Répertoire de base
BASE_DIR="/home/plstenge/coprolites_comparison/07_quality_check_clean"

# Fichier de sortie
OUTPUT_FILE="/home/plstenge/coprolites_comparison/11_summary_tables/recapitulatif_multiqc_9lots.tsv"

echo "=========================================="
echo "Récapitulatif MultiQC - 9 LOTS"
echo "Date: $(date)"
echo "=========================================="

# Définition des 9 lots
declare -a LOTS=(
  "Lot1_Illumina"
  "Lot2_Run1_R1_Ps4_150_no_filtered"
  "Lot3_Run1_R1_Ps6_150_filtered"
  "Lot4_Run2_R1_Ps6_150_no_filtered"
  "Lot5_Run2_R1_Ps6_150_filtered"
  "Lot6_Run3_R2_Ps6_150_no_filtered"
  "Lot7_Run3_R2_Ps8_150_filtered"
  "Lot8_Run4_R2_Ps8_150_filtered"
  "Lot9_Run4_R2_Ps8_75_no_filtered"
)

# Créer l'en-tête du fichier de sortie
echo -e "Lot\tTotal_Sequences_Sum\tTotal_Sequences_Mean\tTotal_Sequences_SD\tAvg_Sequence_Length_Mean\tMedian_Sequence_Length_Mean" > "$OUTPUT_FILE"

# Fonction pour calculer l'écart-type avec awk
calculate_sd() {
    local values=("$@")
    local n=${#values[@]}
    
    if [[ $n -lt 2 ]]; then
        echo "0"
        return
    fi
    
    # Utiliser awk pour calculer la moyenne et l'écart-type
    printf "%s\n" "${values[@]}" | awk -v n="$n" '
    {
        sum += $1
        sumsq += $1 * $1
    }
    END {
        mean = sum / n
        variance = (sumsq - n * mean * mean) / (n - 1)
        if (variance < 0) variance = 0
        sd = sqrt(variance)
        printf "%.2f", sd
    }'
}

# Parcourir chaque lot
for lot in "${LOTS[@]}"; do
    echo ""
    echo "------------------------------------------"
    echo "Traitement du ${lot}"
    echo "------------------------------------------"
    
    # Chemin vers le fichier MultiQC
    MULTIQC_FILE="${BASE_DIR}/${lot}/multiqc_clean_${lot}_data/multiqc_fastqc.txt"
    
    # Vérifier si le fichier existe
    if [[ ! -f "$MULTIQC_FILE" ]]; then
        echo "⚠ Fichier introuvable: ${MULTIQC_FILE}"
        echo -e "${lot}\tNA\tNA\tNA\tNA\tNA" >> "$OUTPUT_FILE"
        continue
    fi
    
    # Extraire les colonnes nécessaires (skip header, garder lignes avec cop)
    DATA_FILE=$(mktemp)
    grep -E "^clean_cop[0-9]+_dedup" "$MULTIQC_FILE" | cut -f1,5,11,12 > "$DATA_FILE"
    
    # Vérifier si on a des données
    if [[ ! -s "$DATA_FILE" ]]; then
        echo "⚠ Aucune donnée valide trouvée pour ${lot}"
        echo -e "${lot}\tNA\tNA\tNA\tNA\tNA" >> "$OUTPUT_FILE"
        rm -f "$DATA_FILE"
        continue
    fi
    
    # Extraire les valeurs dans des tableaux
    mapfile -t total_sequences < <(cut -f2 "$DATA_FILE")
    mapfile -t avg_lengths < <(cut -f3 "$DATA_FILE")
    mapfile -t median_lengths < <(cut -f4 "$DATA_FILE")
    
    # Calculer la somme des Total Sequences
    sum_ts=0
    for val in "${total_sequences[@]}"; do
        sum_ts=$(echo "$sum_ts + $val" | bc)
    done
    
    # Calculer la moyenne des Total Sequences
    n_samples=${#total_sequences[@]}
    mean_ts=$(echo "scale=2; $sum_ts / $n_samples" | bc)
    
    # Calculer l'écart-type des Total Sequences
    sd_ts=$(calculate_sd "${total_sequences[@]}")
    
    # Calculer la moyenne de avg_sequence_length
    sum_avg=0
    for val in "${avg_lengths[@]}"; do
        sum_avg=$(echo "$sum_avg + $val" | bc)
    done
    mean_avg=$(echo "scale=2; $sum_avg / $n_samples" | bc)
    
    # Calculer la moyenne de median_sequence_length
    sum_median=0
    for val in "${median_lengths[@]}"; do
        sum_median=$(echo "$sum_median + $val" | bc)
    done
    mean_median=$(echo "scale=2; $sum_median / $n_samples" | bc)
    
    # Afficher les résultats
    echo "  ✓ ${n_samples} échantillons trouvés"
    echo "  ✓ Somme Total Sequences: $sum_ts"
    echo "  ✓ Moyenne Total Sequences: $mean_ts"
    echo "  ✓ Écart-type Total Sequences: $sd_ts"
    echo "  ✓ Moyenne Avg Length: $mean_avg"
    echo "  ✓ Moyenne Median Length: $mean_median"
    
    # Sauvegarder dans le fichier de sortie
    echo -e "${lot}\t${sum_ts}\t${mean_ts}\t${sd_ts}\t${mean_avg}\t${mean_median}" >> "$OUTPUT_FILE"
    
    # Nettoyage
    rm -f "$DATA_FILE"
done

echo ""
echo "=========================================="
echo "Récapitulatif terminé!"
echo "Date: $(date)"
echo "=========================================="
echo ""
echo "Fichier de sortie: ${OUTPUT_FILE}"
echo ""

# Afficher le résultat
echo "Aperçu du résultat:"
head -n 10 "$OUTPUT_FILE"

# Notification
echo "Récapitulatif MultiQC terminé le $(date)" | \
    mail -s "Pipeline Coprolites - Récapitulatif MultiQC terminé" pierrelouis.stenger@gmail.com 2>/dev/null || true

exit 0
