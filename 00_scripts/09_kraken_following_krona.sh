#!/usr/bin/env bash

cd /home/plstenge/coprolites_comparison/09_kraken_following/combined_kraken_for_mapdamage


module load conda/4.12.0
source ~/.bashrc
conda activate bioinformatic


awk 'BEGIN {OFS="\t"}
{
    best_taxid = 0;
    best_n = -1;

    for (i = 2; i <= NF; i++) {
        split($i, x, ":");

        if (x[1] ~ /^[0-9]+$/ && x[1] != 0 && x[2] ~ /^[0-9]+$/) {
            if (x[2] > best_n) {
                best_taxid = x[1];
                best_n = x[2];
            }
        }
    }

    if (best_taxid > 0) {
        print $1, best_taxid, best_n;
    }
}' cop408_following.kraken > cop408_following.krona.tsv

ktImportTaxonomy \
  -q 1 \
  -t 2 \
  -m 3 \
  -o cop408_following.combined.krona.html \
  cop408_following.krona.tsv
