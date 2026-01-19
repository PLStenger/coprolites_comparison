#!/bin/bash

################################################################################
# DIAGNOSTIC SCRIPT - MapDamage BWA Issue
# Teste chaque étape d'indexation et mapping
################################################################################

set -x  # Affiche chaque commande exécutée

BASE_DIR="/home/plstenge/coprolites_comparison"
GENOME="/home/plstenge/genomes/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa"

# Vérifier que le génome existe
echo "=== 1. Vérification du génome ==="
if [[ ! -f "$GENOME" ]]; then
    echo "✗ ERREUR: Génome non trouvé: $GENOME"
    exit 1
fi

ls -lh "$GENOME"
file "$GENOME" | head -1

# Vérifier sa taille
SIZE=$(wc -c < "$GENOME")
echo "Taille du fichier: $SIZE bytes"

# Vérifier le nombre de séquences
echo ""
echo "=== 2. Nombre de séquences ==="
grep -c "^>" "$GENOME"

# Vérifier les 10 premiers headers
echo ""
echo "=== 3. Premiers headers ==="
grep "^>" "$GENOME" | head -10

# Vérifier les espaces dans les headers
echo ""
echo "=== 4. Headers avec espaces problématiques ==="
grep "^>" "$GENOME" | grep " " | head -5

# Tester l'indexation BWA
echo ""
echo "=== 5. Test indexation BWA ==="
FIXED_GENOME="/home/plstenge/genomes/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fixed.fasta"

# Créer un symlink ou copie
if [[ ! -f "$FIXED_GENOME" ]]; then
    echo "Création de lien pour génome fixé..."
    ln -sf "$GENOME" "$FIXED_GENOME" || cp "$GENOME" "$FIXED_GENOME"
fi

# Activer l'environnement mapdamage
echo ""
echo "=== 6. Activation environnement mapdamage ==="
conda activate mapdamage_py39

which bwa
bwa 2>&1 | head -3

which samtools
samtools --version | head -1

# Essayer l'indexation
echo ""
echo "=== 7. Test indexation (dry run) ==="
echo "Commande: bwa index \"$FIXED_GENOME\""

# Tester avec une petite partie du génome d'abord
echo ""
echo "=== 8. Test avec sous-ensemble ==="
TEST_DIR="/tmp/bwa_test_$$"
mkdir -p "$TEST_DIR"

# Créer un petit fichier de test (100kb)
head -c 100000 "$FIXED_GENOME" > "$TEST_DIR/test.fa"

echo "Fichier test créé: $TEST_DIR/test.fa"
ls -lh "$TEST_DIR/test.fa"

# Essayer l'indexation sur le petit fichier
echo "Test d'indexation sur petit fichier..."
cd "$TEST_DIR"
bwa index test.fa 2>&1 | head -20

if [[ $? -eq 0 ]]; then
    echo "✓ Indexation BWA fonctionne!"
else
    echo "✗ Indexation BWA échoue!"
fi

# Vérifier les index créés
echo ""
echo "=== 9. Index créés ==="
ls -lh "$TEST_DIR"/test.fa* 2>/dev/null || echo "Aucun index créé!"

# Test sur le génome complet
echo ""
echo "=== 10. Test indexation génome complet ==="
cd "$BASE_DIR"
echo "Tentative d'indexation du génome complet (cela peut prendre du temps)..."
timeout 300 bwa index "$FIXED_GENOME" 2>&1 | tail -20

if [[ $? -eq 124 ]]; then
    echo "⚠ TIMEOUT: l'indexation prend trop longtemps"
elif [[ $? -eq 0 ]]; then
    echo "✓ Indexation réussie!"
    ls -lh "${FIXED_GENOME}"* | head -10
else
    echo "✗ Indexation échouée!"
fi

# Vérifier l'espace disque
echo ""
echo "=== 11. Espace disque disponible ==="
df -h /home/plstenge/genomes/ | tail -2

# Vérifier les permissions
echo ""
echo "=== 12. Permissions ==="
ls -ld /home/plstenge/genomes/
ls -l "$FIXED_GENOME"

# Nettoyage
echo ""
echo "=== 13. Nettoyage ==="
rm -rf "$TEST_DIR"

echo ""
echo "=== FIN DIAGNOSTIC ==="
