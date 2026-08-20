#!/usr/bin/env bash
#SBATCH --job-name=14_sort_fasta
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=/home/plstenge/coprolites_comparison/00_scripts/14_sort_fasta.err
#SBATCH --output=/home/plstenge/coprolites_comparison/00_scripts/14_sort_fasta.out

input="/home/plstenge/coprolites_comparison/13_blast_local_existing_kraken_fastq/cop408/Alces/cop408_Alces_all_kraken.fasta"
output="/home/plstenge/coprolites_comparison/13_blast_local_existing_kraken_fastq/cop408/Alces/cop408_Alces_all_kraken.sorted_by_length_desc.fasta"

awk '
  /^>/ {
    if (header != "") {
      print length(seq) "\t" header "\t" seq
    }
    header = $0
    seq = ""
    next
  }
  {
    seq = seq $0
  }
  END {
    if (header != "") {
      print length(seq) "\t" header "\t" seq
    }
  }
' "$input" |
sort -rn -k1,1 |
cut -f2- |
awk -F "\t" '{print $1 "\n" $2}' \
> "$output"
