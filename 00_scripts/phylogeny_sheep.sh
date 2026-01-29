#!/bin/bash

################################################################################
# PHYLOGÉNIE ADN ANCIEN MOUTON - IDENTIFICATION DE RACE
# Récupération reads → Consensus poolé → Phylogénie + IDENTIFICATION
# 
# VERSION: Optimisée pour identification de race (bootstrap robuste)
# OBJECTIF: Identifier la race de mouton la plus proche des échantillons anciens
# AUTEUR: Bioinformatique Paléo
# DATE: 2026-01-29
################################################################################

#SBATCH --job-name=sheep_race_identification
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=200G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/phylogeny_identification.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/phylogeny_identification.out"

set -euo pipefail

echo ""
echo "=========================================="
echo "IDENTIFICATION DE RACE - MOUTON ANCIEN"
echo "=========================================="
echo ""
echo "Timestamp: $(date)"
echo ""

# ========== INITIALISATION CONDA ==========
echo "Initialisation des environnements conda..."

module load conda/4.12.0
source ~/.bashrc

# Créer environnement s'il n'existe pas
if ! conda env list | grep -q "phylogeny_env"; then
    echo "Création environnement phylogeny_env..."
    conda create -y -n phylogeny_env \
        -c bioconda -c conda-forge \
        bcftools samtools mafft raxml \
        python=3.10 biopython numpy pandas scipy
fi

conda activate phylogeny_env
echo "✓ Environnement phylogeny_env activé"
echo ""

# ========== CONFIGURATION GLOBALE ==========
BASE_DIR="/home/plstenge/coprolites_comparison"
BAM_EXTRACTION_DIR="${BASE_DIR}/13_phylogeny/01_bam_extraction"
CONSENSUS_DIR="${BASE_DIR}/13_phylogeny/02_consensus"
ALIGN_DIR="${BASE_DIR}/13_phylogeny/03_alignment"
PHYLOGENY_DIR="${BASE_DIR}/13_phylogeny/04_phylogeny"
VISUALIZATION_DIR="${BASE_DIR}/13_phylogeny/05_visualization"
REFERENCES_DIR="${BASE_DIR}/13_phylogeny/00_references"
LOG_DIR="${BASE_DIR}/13_phylogeny/logs"
IDENTIFICATION_DIR="${BASE_DIR}/13_phylogeny/identification_results"

THREADS=8
MOUFLON_REF="/home/plstenge/genomes/Ovis_aries_mitochondrion/Ovis_aries_mitochondrion.fa"

echo "Configuration:"
echo "  BASE_DIR: $BASE_DIR"
echo "  REFERENCES: $REFERENCES_DIR"
echo "  IDENTIFICATION: $IDENTIFICATION_DIR"
echo "  THREADS: $THREADS"
echo ""

# ========== CRÉATION RÉPERTOIRES ==========
mkdir -p "$BAM_EXTRACTION_DIR" "$CONSENSUS_DIR" "$ALIGN_DIR" "$PHYLOGENY_DIR" \
         "$VISUALIZATION_DIR" "$REFERENCES_DIR" "$LOG_DIR" "$IDENTIFICATION_DIR"

LOGFILE="${LOG_DIR}/identification_$(date +%Y%m%d_%H%M%S).log"
touch "$LOGFILE"

# ========== VÉRIFICATION DES OUTILS ==========
echo "=========================================="
echo "Vérification des outils"
echo "=========================================="
echo ""

for tool in samtools bcftools mafft raxmlHPC-PTHREADS python3; do
    if command -v "$tool" &>/dev/null; then
        echo "  ✓ $tool"
    else
        echo "  ✗ $tool MANQUANT"
        if [[ "$tool" == "mafft" || "$tool" == "raxmlHPC-PTHREADS" ]]; then
            exit 1
        fi
    fi
done
echo ""

################################################################################
# ÉTAPE 1: CRÉER LE CONSENSUS POOLÉ (UNIQUE)
################################################################################

echo "=========================================="
echo "ÉTAPE 1: Création consensus poolé UNIQUE"
echo "=========================================="
echo ""

# Chercher tous les BAM mouton individuels
SHEEP_BAMS=()
for bam_file in "${BAM_EXTRACTION_DIR}"/*Ovis_aries*.bam; do
    if [[ -f "$bam_file" ]]; then
        SHEEP_BAMS+=("$bam_file")
    fi
done

if [[ ${#SHEEP_BAMS[@]} -eq 0 ]]; then
    echo "✗ ERREUR: Aucun BAM mouton trouvé dans $BAM_EXTRACTION_DIR"
    exit 1
fi

echo "BAM mouton trouvés: ${#SHEEP_BAMS[@]}"
TOTAL_READS=0
for bam in "${SHEEP_BAMS[@]}"; do
    SAMPLE=$(basename "$bam" .bam)
    READS=$(samtools view -c "$bam" 2>/dev/null || echo "0")
    TOTAL_READS=$((TOTAL_READS + READS))
    echo "  → $SAMPLE: $READS reads"
done
echo "  Reads total: $TOTAL_READS"
echo ""

# Merger tous les BAM en un seul
POOLED_BAM="${BAM_EXTRACTION_DIR}/POOLED_ancient_sheep.bam"

if [[ -f "$POOLED_BAM" ]]; then
    echo "Réutilisation du BAM poolé existant: $POOLED_BAM"
else
    echo "Fusion des BAM en un seul échantillon poolé..."
    samtools merge -@ $THREADS -f "$POOLED_BAM" "${SHEEP_BAMS[@]}"
    samtools index "$POOLED_BAM"
    echo "✓ BAM poolé créé et indexé"
fi

# Statistiques du pooled
TOTAL_READS=$(samtools view -c "$POOLED_BAM")
MAPPED_READS=$(samtools view -c -F 4 "$POOLED_BAM")
COV_MEAN=$(samtools depth "$POOLED_BAM" 2>/dev/null | \
    awk '{sum+=$3; count++} END {if (count>0) printf "%.1f", sum/count; else print "0"}')

echo ""
echo "Statistiques BAM poolé final:"
echo "  Total reads: $TOTAL_READS"
echo "  Mapped reads: $MAPPED_READS"
echo "  Couverture moyenne: ${COV_MEAN}X"
echo ""

################################################################################
# ÉTAPE 2: GÉNÉRER CONSENSUS DEPUIS LE BAM POOLÉ
################################################################################

echo "=========================================="
echo "ÉTAPE 2: Génération consensus"
echo "=========================================="
echo ""

CONSENSUS_RAW="${CONSENSUS_DIR}/pooled_consensus_raw.txt"
CONSENSUS_FA="${CONSENSUS_DIR}/POOLED_ancient_sheep.fa"

if [[ -f "$CONSENSUS_FA" ]]; then
    echo "Réutilisation du consensus existant"
else
    echo "Génération du consensus par appel de majorité..."
    
    samtools mpileup -f "$MOUFLON_REF" "$POOLED_BAM" 2>/dev/null | \
        awk '{
            depth=$4
            bases=$5
            if (depth>=1) {
                gsub(/[^ACGTacgt]/,"",bases)
                a=gsub(/[aA]/,"",bases)
                c=gsub(/[cC]/,"",bases)  
                g=gsub(/[gG]/,"",bases)
                t=gsub(/[tT]/,"",bases)
                
                if (a>=c && a>=g && a>=t) printf "A"
                else if (c>=a && c>=g && c>=t) printf "C"
                else if (g>=a && g>=c && g>=t) printf "G"
                else printf "T"
            } else {
                printf "N"
            }
        }' > "$CONSENSUS_RAW"
    
    # Créer FASTA formaté
    cat > "$CONSENSUS_FA" << EOF
>POOLED_ancient_sheep
EOF
    
    cat "$CONSENSUS_RAW" | sed 's/\(.\{80\}\)/\1\n/g' >> "$CONSENSUS_FA"
fi

# Statistiques consensus
CONS_LENGTH=$(awk '/^>/ {next} {total+=length($0)} END {print total}' "$CONSENSUS_FA")
N_COUNT=$(grep -v "^>" "$CONSENSUS_FA" | grep -o 'N' | wc -l)
PERCENT_N=$(echo "scale=2; $N_COUNT * 100 / $CONS_LENGTH" | bc 2>/dev/null || echo "0")

echo "Consensus généré:"
echo "  Longueur: $CONS_LENGTH bp"
echo "  Bases N: $N_COUNT ($PERCENT_N%)"
echo "  Déterminées: $((CONS_LENGTH - N_COUNT)) bp"
echo ""

################################################################################
# ÉTAPE 3: CRÉER L'ALIGNEMENT MULTIPLE
################################################################################

echo "=========================================="
echo "ÉTAPE 3: Alignement multiple MAFFT"
echo "=========================================="
echo ""

COMBINED_FA="${ALIGN_DIR}/combined_references_and_ancient.fasta"
ALIGNED_FA="${ALIGN_DIR}/aligned_identification.fasta"

if [[ ! -f "$COMBINED_FA" ]]; then
    echo "Création fichier combiné..."
    
    if [[ -f "${REFERENCES_DIR}/sheep_reference_breeds_mtDNA.fasta" ]]; then
        # Combiner POOLED + références
        cat "$CONSENSUS_FA" "${REFERENCES_DIR}/sheep_reference_breeds_mtDNA.fasta" > "$COMBINED_FA"
        N_SEQS=$(grep -c '^>' "$COMBINED_FA")
        echo "✓ Combiné: 1 ancien + $((N_SEQS-1)) références = $N_SEQS séquences"
    else
        echo "✗ ERREUR: Fichier références manquant"
        exit 1
    fi
    
    echo ""
    echo "Composition:"
    grep "^>" "$COMBINED_FA" | head -5
    echo "..."
    echo ""
fi

# Aligner avec MAFFT
if [[ ! -f "$ALIGNED_FA" ]]; then
    echo "Alignement MAFFT (--auto --thread $THREADS)..."
    
    mafft --auto --thread $THREADS "$COMBINED_FA" > "$ALIGNED_FA"
    
    ALIGN_SEQS=$(grep -c "^>" "$ALIGNED_FA")
    ALIGN_LEN=$(head -n 2 "$ALIGNED_FA" | tail -n 1 | wc -c)
    
    echo "✓ Alignement généré"
    echo "  Séquences: $ALIGN_SEQS"
    echo "  Longueur: $ALIGN_LEN bp"
else
    echo "Réutilisation de l'alignement existant"
fi

echo ""

################################################################################
# ÉTAPE 4: CONSTRUIRE L'ARBRE PHYLOGÉNÉTIQUE RAxML
################################################################################

echo "=========================================="
echo "ÉTAPE 4: Phylogénie RAxML (100 bootstraps)"
echo "=========================================="
echo ""

cd "$PHYLOGENY_DIR"

# Nettoyer anciens fichiers
rm -f RAxML_* 2>/dev/null || true

echo "Paramètres RAxML:"
echo "  Modèle: GTRGAMMA (modèle standard pour mtDNA)"
echo "  Bootstraps: 100 (bootstrap robuste)"
echo "  Threads: $THREADS"
echo "  Outgroup: Capra_hircus_mitochondrion (chèvre)"
echo ""
echo "Durée estimée: 30-60 minutes"
echo ""

RUN_ID="sheep_identification_$(date +%Y%m%d_%H%M%S)"

# Lancer RAxML avec outgroup explicite
raxmlHPC-PTHREADS \
    -f a \
    -x 12345 \
    -p 12345 \
    -# 100 \
    -m GTRGAMMA \
    -s "$ALIGNED_FA" \
    -n "$RUN_ID" \
    -T $THREADS \
    -o Capra_hircus_mitochondrion \
    2>&1 | tee "${LOG_DIR}/raxml_${RUN_ID}.log"

echo ""
echo "✓ RAxML terminé"
echo ""

# Vérifier résultats
BEST_TREE="RAxML_bipartitions.${RUN_ID}"
if [[ ! -f "$BEST_TREE" ]]; then
    echo "✗ ERREUR: Arbre RAxML non généré"
    exit 1
fi

echo "Fichiers générés:"
ls -lh RAxML_* | awk '{print "  " $9 " (" $5 ")"}'
echo ""
echo "✓ Arbre phylogénétique: $BEST_TREE"
echo ""

################################################################################
# ÉTAPE 5: ANALYSE DE DISTANCES GÉNÉTIQUES ET IDENTIFICATION
################################################################################

echo "=========================================="
echo "ÉTAPE 5: Analyse distances et identification"
echo "=========================================="
echo ""

# Créer script Python pour l'analyse
ANALYSIS_SCRIPT="${IDENTIFICATION_DIR}/analyze_phylogeny.py"

cat > "$ANALYSIS_SCRIPT" << 'PYTHON_SCRIPT'
#!/usr/bin/env python3
"""
Script d'identification de race de mouton ancien
Analyse distances génétiques K2P et bootstrap support
"""

import sys
import re
from collections import defaultdict
from pathlib import Path

# Paramètres
ALIGNED_FA = sys.argv[1] if len(sys.argv) > 1 else None
TREE_FILE = sys.argv[2] if len(sys.argv) > 2 else None
OUTPUT_DIR = sys.argv[3] if len(sys.argv) > 3 else "."

if not ALIGNED_FA or not TREE_FILE:
    print("Usage: python analyze_phylogeny.py <aligned.fasta> <tree.nwk> <output_dir>")
    sys.exit(1)

print("\n" + "="*60)
print("IDENTIFICATION DE RACE - MOUTON ANCIEN")
print("="*60 + "\n")

# ========== ÉTAPE 1: CHARGER L'ALIGNEMENT ==========

from Bio import SeqIO, AlignIO

print("1️⃣  ÉTAPE 1: Chargement de l'alignement")
print("-" * 60)

sequences = {}
try:
    alignment = AlignIO.read(ALIGNED_FA, "fasta")
    for record in alignment:
        sequences[record.id] = str(record.seq).upper()
    print(f"✓ Alignement chargé: {len(sequences)} séquences")
    print(f"  Longueur: {len(alignment[0])} bp")
except Exception as e:
    print(f"✗ Erreur lecture alignement: {e}")
    sys.exit(1)

# Identifier la séquence ancienne
ancient_id = None
for seq_id in sequences.keys():
    if "ancient" in seq_id.lower() or "pooled" in seq_id.lower():
        ancient_id = seq_id
        break

if not ancient_id:
    print("✗ Séquence ancienne non trouvée (cherche 'ancient' ou 'pooled')")
    sys.exit(1)

ancient_seq = sequences[ancient_id]
print(f"✓ Séquence ancienne trouvée: {ancient_id}")
print(f"  Longueur: {len(ancient_seq)} bp")
print()

# ========== ÉTAPE 2: CALCULER DISTANCES K2P ==========

print("2️⃣  ÉTAPE 2: Calcul distances génétiques K2P")
print("-" * 60)

def calculate_k2p_distance(seq1, seq2):
    """Calcul distance Kimura 2-parameter"""
    if len(seq1) != len(seq2):
        return None
    
    n = len(seq1)
    if n == 0:
        return None
    
    # Compter transitions et transversions
    transitions = 0  # A↔G, C↔T
    transversions = 0  # A↔C, A↔T, G↔C, G↔T
    differences = 0
    
    for i in range(n):
        if seq1[i] == seq2[i]:
            continue
        differences += 1
        
        # Purines: A, G; Pyrimidines: C, T
        if (seq1[i] in 'AG' and seq2[i] in 'AG') or (seq1[i] in 'CT' and seq2[i] in 'CT'):
            transitions += 1
        elif seq1[i] in 'ACGT' and seq2[i] in 'ACGT':
            transversions += 1
    
    p = transitions / n  # Fréquence transitions
    q = transversions / n  # Fréquence transversions
    
    # Formule K2P
    if p >= 0.5 or q >= 0.25:
        return None  # Divergence trop élevée
    
    k2p = -0.5 * math.log(1 - 2*p - q) - 0.25 * math.log(1 - 2*q)
    return k2p

import math

# Calculer distances
distances = {}
for seq_id, seq in sequences.items():
    if seq_id == ancient_id:
        continue
    
    # Ignorer outgroup (chèvre)
    if "capra" in seq_id.lower():
        continue
    
    # Ignorer séquences trop courtes
    if len(seq) < 1000:
        continue
    
    dist = calculate_k2p_distance(ancient_seq, seq)
    if dist is not None:
        distances[seq_id] = dist

# Trier par distance
sorted_distances = sorted(distances.items(), key=lambda x: x[1])

print(f"✓ Distances K2P calculées: {len(sorted_distances)} races")
print()

# ========== ÉTAPE 3: ANALYSER BOOTSTRAP ==========

print("3️⃣  ÉTAPE 3: Analyse support bootstrap")
print("-" * 60)

# Lire l'arbre
try:
    with open(TREE_FILE, 'r') as f:
        tree_content = f.read().strip()
    print(f"✓ Arbre phylogénétique chargé: {TREE_FILE}")
except Exception as e:
    print(f"⚠ Erreur lecture arbre: {e}")
    tree_content = None

print()

# ========== ÉTAPE 4: RAPPORT D'IDENTIFICATION ==========

print("4️⃣  ÉTAPE 4: RAPPORT D'IDENTIFICATION")
print("="*60)
print()

# Créer rapport
report = []
report.append("RAPPORT D'IDENTIFICATION DE RACE")
report.append("="*60)
report.append("")
report.append(f"Échantillon: {ancient_id}")
report.append(f"Longueur génome: {len(ancient_seq)} bp")
report.append("")
report.append("RÉSULTATS PHYLOGÉNÉTIQUES")
report.append("-"*60)
report.append("")
report.append("TOP 10 RACES LES PLUS PROCHES (basé sur distance génétique K2P)")
report.append("")

top_races = sorted_distances[:10]

for rank, (race_id, dist) in enumerate(top_races, 1):
    # Nettoyer le nom
    race_name = race_id.replace('Ovis_aries_breed_', '').replace('_mitochondrion', '').replace('_', ' ')
    
    # Interprétation distance
    if dist < 0.001:
        similarity = "IDENTIQUE/Très proche"
    elif dist < 0.005:
        similarity = "Très proche"
    elif dist < 0.01:
        similarity = "Proche"
    elif dist < 0.02:
        similarity = "Distance modérée"
    else:
        similarity = "Distant"
    
    print(f"Rang {rank}: {race_name}")
    print(f"  Distance K2P: {dist:.6f}")
    print(f"  Interprétation: {similarity}")
    print()
    
    report.append(f"Rang {rank}: {race_name}")
    report.append(f"  Distance K2P: {dist:.6f} ({similarity})")
    report.append("")

# Identification finale
if sorted_distances:
    closest_race = sorted_distances[0]
    race_name = closest_race[0].replace('Ovis_aries_breed_', '').replace('_mitochondrion', '').replace('_', ' ')
    dist = closest_race[1]
    
    report.append("")
    report.append("IDENTIFICATION FINALE")
    report.append("="*60)
    report.append(f"Race la plus proche: {race_name}")
    report.append(f"Distance génétique: {dist:.6f}")
    report.append("")
    
    if dist < 0.001:
        confidence = "TRÈS ÉLEVÉE"
        interpretation = "Identique ou quasi-identique à la race moderne"
    elif dist < 0.005:
        confidence = "ÉLEVÉE"
        interpretation = "Très proche, vraisemblablement même race ou race ancestrale"
    elif dist < 0.01:
        confidence = "MODÉRÉE"
        interpretation = "Proche, lignée maternelle connexe"
    else:
        confidence = "FAIBLE"
        interpretation = "Distance génétique modérée, race différente"
    
    report.append(f"Confiance: {confidence}")
    report.append(f"Interprétation: {interpretation}")
    report.append("")

# Contexte biologique
report.append("")
report.append("CONTEXTE BIOLOGIQUE")
report.append("-"*60)
report.append("L'ADN mitochondrial représente seulement l'hérédité MATERNELLE")
report.append("(~16.5 kb de gène cytochrome b et ARNt)")
report.append("")
report.append("Les distances K2P reflètent:")
report.append("• Temps de divergence maternelle")
report.append("• Isolement reproductif historique")
report.append("• Flux génétique entre races")
report.append("")

# Sauvegarder rapport
report_text = "\n".join(report)
report_file = Path(OUTPUT_DIR) / "identification_report.txt"
with open(report_file, 'w') as f:
    f.write(report_text)

print(report_text)
print("")
print(f"✓ Rapport sauvegardé: {report_file}")

# Sauvegarder résultats complets
tsv_file = Path(OUTPUT_DIR) / "closest_relatives.tsv"
with open(tsv_file, 'w') as f:
    f.write("Rank\tRace_Name\tDistance_K2P\tInterpretation\n")
    for rank, (race_id, dist) in enumerate(sorted_distances[:20], 1):
        race_name = race_id.replace('Ovis_aries_breed_', '').replace('_mitochondrion', '').replace('_', ' ')
        
        if dist < 0.001:
            sim = "Identique/Très proche"
        elif dist < 0.005:
            sim = "Très proche"
        elif dist < 0.01:
            sim = "Proche"
        elif dist < 0.02:
            sim = "Distance modérée"
        else:
            sim = "Distant"
        
        f.write(f"{rank}\t{race_name}\t{dist:.6f}\t{sim}\n")

print(f"✓ Résultats TSV: {tsv_file}")
print("")

PYTHON_SCRIPT

# Exécuter l'analyse
echo "Exécution de l'analyse..."
python3 "$ANALYSIS_SCRIPT" "$ALIGNED_FA" "$PHYLOGENY_DIR/$BEST_TREE" "$IDENTIFICATION_DIR"

echo ""

################################################################################
# ÉTAPE 6: PRÉPARATION VISUALISATION
################################################################################

echo "=========================================="
echo "ÉTAPE 6: Préparation visualisation FigTree"
echo "=========================================="
echo ""

mkdir -p "${VISUALIZATION_DIR}/phylogeny"

# Copier l'arbre
cp "$PHYLOGENY_DIR/$BEST_TREE" "${VISUALIZATION_DIR}/phylogeny/sheep_identification.nwk"

echo "Pour visualiser l'arbre:"
echo ""
echo "1. Ouvrir FigTree:"
echo "   figtree ${VISUALIZATION_DIR}/phylogeny/sheep_identification.nwk"
echo ""
echo "2. Dans FigTree:"
echo "   - Node Labels → ☑ Label (voir bootstrap)"
echo "   - Appearance → ☑ Colour by partition"
echo "   - Scale Bar → ☑ Show (ajouter l'échelle)"
echo "   - Layout → Radial (vue circulaire)"
echo ""
echo "3. Identifier POOLED_ancient_sheep dans l'arbre"
echo "   - Chercher ses voisins les plus proches"
echo "   - Vérifier leur bootstrap support (≥70 = confiant)"
echo ""

################################################################################
# RÉSUMÉ FINAL
################################################################################

echo ""
echo "=========================================="
echo "✓ IDENTIFICATION TERMINÉE"
echo "=========================================="
echo ""

echo "📊 RÉSULTATS GÉNÉRÉS:"
echo ""
echo "1. Rapport texte complet:"
echo "   ${IDENTIFICATION_DIR}/identification_report.txt"
echo ""
echo "2. Tableau races classées:"
echo "   ${IDENTIFICATION_DIR}/closest_relatives.tsv"
echo ""
echo "3. Arbre phylogénétique (FigTree):"
echo "   ${VISUALIZATION_DIR}/phylogeny/sheep_identification.nwk"
echo ""
echo "4. Fichiers d'alignement:"
echo "   ${ALIGN_DIR}/aligned_identification.fasta"
echo ""
echo "5. Log complet:"
echo "   $LOGFILE"
echo ""

echo "📈 PROCHAINES ÉTAPES:"
echo ""
echo "1. Lire le rapport d'identification:"
echo "   cat ${IDENTIFICATION_DIR}/identification_report.txt"
echo ""
echo "2. Visualiser l'arbre:"
echo "   figtree ${VISUALIZATION_DIR}/phylogeny/sheep_identification.nwk"
echo ""
echo "3. Comparer avec autres données (isotopes, archéologie)"
echo ""
echo "4. Publication des résultats"
echo ""

echo "Timestamp: $(date)"
echo "=========================================="
echo ""

conda deactivate

exit 0
