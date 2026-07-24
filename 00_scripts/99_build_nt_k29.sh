#!/bin/bash

#SBATCH --job-name=build_nt_k29_PL
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=36
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/build_nt_k29_PL.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/build_nt_k29_PL.out"

#===============================================================================
# build_nt_k29_PL.sh
# Construction d'une base Kraken2 (+ fichiers requis pour Bracken) k=29
# à partir du fichier core_nt local de l'utilisateur.
# Basé sur la méthodologie d'Olivier Rué (Migale, INRAE)
#===============================================================================
#set -euo pipefail
#set -o errtrace

#-------------------------------------------------------------------------------
# CONFIGURATION (à adapter si besoin)
#-------------------------------------------------------------------------------
DB_NAME="nt_k29_PL"
FINAL_DEST="/storage/groups/gdec/shared/Kraken_database/${DB_NAME}"

SOURCE_FASTA="/storage/biodatabanks/ncbi/NT/current/fasta/All/all.fasta"

KMER_LEN=29
MINIMIZER_LEN=23
MINIMIZER_SPACES=2

THREADS=${SLURM_CPUS_PER_TASK:-36}
# Espace de travail rapide (RAM disque). Vérifié dynamiquement plus bas.
WORKDIR="/dev/shm/${DB_NAME}_build"
# Répertoire de secours si /dev/shm est trop petit (disque rapide local)
FALLBACK_WORKDIR="/tmp/${DB_NAME}_build"

module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics

CONDA_ENV="metagenomics"
BRACKEN_ENV="metagenomics"          # env conda contenant bracken-build, sinon on l'installe dans kraken2 env
BRACKEN_READ_LENS=(50 75 100 150)   # longueurs de reads à préparer pour Bracken

LOGFILE="build_${DB_NAME}_$(date +%Y%m%d_%H%M%S).log"

#-------------------------------------------------------------------------------
# UTILS
#-------------------------------------------------------------------------------
log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"; }
die()  { log "ERREUR: $*"; exit 1; }
trap 'log "Échec à la ligne $LINENO. Voir $LOGFILE."' ERR

human_bytes() { numfmt --to=iec-i --suffix=B "$1" 2>/dev/null || echo "$1"; }

check_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Commande requise introuvable: $1"
}

#-------------------------------------------------------------------------------
# 0. VÉRIFICATIONS PRÉALABLES
#-------------------------------------------------------------------------------
log "=== Étape 0: Vérifications préalables ==="

[[ -f "$SOURCE_FASTA" ]] || die "Fichier FASTA source introuvable: $SOURCE_FASTA"

check_cmd awk
check_cmd conda
check_cmd df
check_cmd free
check_cmd numfmt

# Activation conda (nécessite conda déjà initialisé dans le shell)
source "$(conda info --base)/etc/profile.d/conda.sh"

conda activate "$CONDA_ENV" || die "Impossible d'activer l'environnement conda '$CONDA_ENV'"
check_cmd k2

# --- Vérification RAM disponible ---
# Le tuto rapporte ~341 GB pour k=29 sur core_nt complet avec ces paramètres.
# On prévoit une marge de sécurité (facteur 1.3) car la RAM est aussi utilisée
# pour /dev/shm ET pour le processus k2 build lui-même.
MEM_ESTIMATE_GB=341
MEM_SAFETY_GB=$(( MEM_ESTIMATE_GB * 13 / 10 ))   # ~443 GB requis recommandé

AVAIL_MEM_KB=$(awk '/MemAvailable/ {print $2}' /proc/meminfo)
AVAIL_MEM_GB=$(( AVAIL_MEM_KB / 1024 / 1024 ))

log "Mémoire disponible détectée: ${AVAIL_MEM_GB} GB"
log "Estimation mémoire nécessaire pour k=${KMER_LEN} (d'après le tuto sur core_nt complet): ~${MEM_ESTIMATE_GB} GB (hash table), marge recommandée: ~${MEM_SAFETY_GB} GB"

if (( AVAIL_MEM_GB < MEM_ESTIMATE_GB )); then
    die "Mémoire disponible (${AVAIL_MEM_GB} GB) insuffisante par rapport à l'estimation (${MEM_ESTIMATE_GB} GB). Le fichier all.fasta de l'utilisateur étant probablement proche de la taille du core_nt, réduisez --minimizer-spaces, utilisez --max-db-size, ou lancez sur un nœud avec plus de RAM."
elif (( AVAIL_MEM_GB < MEM_SAFETY_GB )); then
    log "ATTENTION: mémoire disponible sous la marge de sécurité recommandée. Le build peut échouer (OOM) selon la taille réelle de all.fasta. Poursuite mais surveillez la RAM (voir Étape 5)."
fi

# --- Choix du répertoire de travail (RAM disque vs disque classique) ---
SHM_AVAIL_KB=$(df --output=avail /dev/shm 2>/dev/null | tail -1 || echo 0)
SHM_AVAIL_GB=$(( SHM_AVAIL_KB / 1024 / 1024 ))
SOURCE_SIZE_BYTES=$(stat -c%s "$SOURCE_FASTA")
SOURCE_SIZE_GB=$(( SOURCE_SIZE_BYTES / 1024 / 1024 / 1024 ))

log "Taille du FASTA source: ${SOURCE_SIZE_GB} GB. Espace disponible /dev/shm: ${SHM_AVAIL_GB} GB."

# On a besoin d'espace pour: taxonomy (~15GB dézippé) + hash table temporaire + FASTA (lien symbolique, négligeable)
NEEDED_SHM_GB=$(( MEM_ESTIMATE_GB + 30 ))

if (( SHM_AVAIL_GB >= NEEDED_SHM_GB )); then
    log "Utilisation de /dev/shm comme répertoire de travail (build ~2x plus rapide selon le tuto)."
else
    log "Espace /dev/shm insuffisant (${SHM_AVAIL_GB} GB < ${NEEDED_SHM_GB} GB requis). Bascule sur ${FALLBACK_WORKDIR}."
    WORKDIR="$FALLBACK_WORKDIR"
    DISK_AVAIL_KB=$(df --output=avail "$(dirname "$WORKDIR")" | tail -1)
    DISK_AVAIL_GB=$(( DISK_AVAIL_KB / 1024 / 1024 ))
    (( DISK_AVAIL_GB >= NEEDED_SHM_GB )) || die "Espace disque insuffisant sur $(dirname "$WORKDIR") (${DISK_AVAIL_GB} GB dispo, ${NEEDED_SHM_GB} GB requis)."
fi

mkdir -p "$WORKDIR"
DB="${WORKDIR}/${DB_NAME}"
mkdir -p "${DB}/library/added" "${DB}/taxonomy"

log "Répertoire de travail: $WORKDIR"
log "Base en construction: $DB"

#-------------------------------------------------------------------------------
# 1. TÉLÉCHARGEMENT DE LA TAXONOMIE
#-------------------------------------------------------------------------------
log "=== Étape 1: Téléchargement de la taxonomie NCBI ==="

if [[ ! -s "${DB}/taxonomy/nodes.dmp" || ! -s "${DB}/taxonomy/names.dmp" ]]; then
    die "Taxonomie absente dans ${DB}/taxonomy, mais le nœud n'a pas Internet. Copie la taxonomie depuis un téléchargement effectué en amont."
fi

[[ -f "${DB}/taxonomy/nucl_gb.accession2taxid" ]] || die "nucl_gb.accession2taxid manquant après download-taxonomy"
[[ -f "${DB}/taxonomy/nucl_wgs.accession2taxid" ]] || die "nucl_wgs.accession2taxid manquant après download-taxonomy"

#-------------------------------------------------------------------------------
# 2. GÉNÉRATION DU FICHIER seqid2taxid.map
#-------------------------------------------------------------------------------
log "=== Étape 2: Génération de seqid2taxid.map ==="

MAP_FILE="${DB}/seqid2taxid.map"

generate_map() {
    local GB="$1" WGS="$2" FASTA="$3" OUT="$4"
    awk '
    BEGIN { FS="\t"; OFS="\t"; }
    ARGIND <= 2 {
        if ($2 != "accession.version") { map[$2] = $3; }
        next;
    }
    /^>/ {
        full_id = $0;
        sub(/^>/, "", full_id);
        split(full_id, parts, /[ |]/);
        id = parts[1];
        if (id in map) { print id, map[id]; }
        else { print id, "1"; }   # défaut = root si accession absente des mappings
    }' "$GB" "$WGS" "$FASTA" > "$OUT"
}

if [[ -s "$MAP_FILE" ]]; then
    log "seqid2taxid.map déjà présent, skip génération."
else
    generate_map "${DB}/taxonomy/nucl_gb.accession2taxid" \
                 "${DB}/taxonomy/nucl_wgs.accession2taxid" \
                 "$SOURCE_FASTA" \
                 "$MAP_FILE"
fi

# --- Vérification cohérence nb séquences vs nb lignes du mapping ---
N_SEQ=$(grep -c '^>' "$SOURCE_FASTA")
N_MAP=$(wc -l < "$MAP_FILE")

log "Nombre de séquences dans le FASTA source: $N_SEQ"
log "Nombre de lignes dans seqid2taxid.map: $N_MAP"

if (( N_SEQ != N_MAP )); then
    die "Incohérence entre le nombre de séquences ($N_SEQ) et le nombre de lignes du mapping ($N_MAP). Vérifiez le FASTA et le script awk avant de continuer."
fi
log "Cohérence FASTA / mapping validée ($N_SEQ séquences)."

#-------------------------------------------------------------------------------
# 3. PRÉPARATION DU FASTA (lien .fna requis par k2 build)
#-------------------------------------------------------------------------------
log "=== Étape 3: Préparation du fichier .fna ==="

LINK_PATH="${DB}/library/added/file.fna"
if [[ -L "$LINK_PATH" || -f "$LINK_PATH" ]]; then
    log "Lien .fna déjà présent."
else
    ln -s "$SOURCE_FASTA" "$LINK_PATH"
fi
echo "file.fna" > "${DB}/library/added/manifest.txt"

#-------------------------------------------------------------------------------
# 4. BUILD KRAKEN2 (k=29)
#-------------------------------------------------------------------------------
log "=== Étape 4: Build Kraken2 k=${KMER_LEN} (threads=${THREADS}) ==="
log "Ceci peut prendre 1 à 4 jours selon la taille réelle de all.fasta (le tuto rapporte ~3 jours pour ~900GB avec 64 threads)."

START_TIME=$(date +%s)

# Lancement en arrière-plan avec nohup pour résister à une déconnexion SSH,
# tout en gardant le suivi dans ce script via wait.
k2 build \
    --threads "$THREADS" \
    --db "$DB" \
    --kmer-len "$KMER_LEN" \
    --minimizer-len "$MINIMIZER_LEN" \
    --minimizer-spaces "$MINIMIZER_SPACES" \
    2>&1 | tee -a "$LOGFILE" || {
        # Le tuto montre qu'une exception "OSError: Errno 5" sur le thread
        # de lecture stderr apparaît systématiquement en fin de build mais
        # N'EST PAS bloquante : on vérifie donc la présence effective des
        # fichiers .k2d avant de considérer l'étape en échec.
        log "k2 build a retourné un code non-nul (souvent le faux-positif 'OSError Errno 5' documenté dans le tuto). Vérification des fichiers de sortie..."
    }

END_TIME=$(date +%s)
log "Durée du build: $(( (END_TIME - START_TIME) / 3600 ))h $(( ((END_TIME - START_TIME) % 3600) / 60 ))min"

for f in hash.k2d opts.k2d taxo.k2d; do
    [[ -s "${DB}/${f}" ]] || die "Fichier attendu manquant après build: ${DB}/${f}. Le build a probablement réellement échoué (voir $LOGFILE)."
done
log "Build Kraken2 terminé avec succès (hash.k2d, opts.k2d, taxo.k2d présents)."

#-------------------------------------------------------------------------------
# 5. FICHIERS POUR BRACKEN
#-------------------------------------------------------------------------------
log "=== Étape 5: Génération des fichiers Bracken ==="

check_cmd bracken-build || log "ATTENTION: bracken-build introuvable dans l'env courant, tentative d'activation de l'env '$BRACKEN_ENV'."
if ! command -v bracken-build >/dev/null 2>&1; then
    conda deactivate
    conda activate "$BRACKEN_ENV" 2>/dev/null || die "bracken-build introuvable. Installez Bracken (conda install -c bioconda bracken) et relancez cette étape manuellement."
fi

for RLEN in "${BRACKEN_READ_LENS[@]}"; do
    log "-- bracken-build pour read length = ${RLEN} --"
    bracken-build \
        -d "$DB" \
        -t "$THREADS" \
        -k "$KMER_LEN" \
        -l "$RLEN" \
        2>&1 | tee -a "$LOGFILE"
done

# Retour à l'env kraken2 pour la suite si besoin
conda activate "$CONDA_ENV" 2>/dev/null || true

#-------------------------------------------------------------------------------
# 6. COPIE VERS LA DESTINATION FINALE
#-------------------------------------------------------------------------------
log "=== Étape 6: Copie vers ${FINAL_DEST} ==="

mkdir -p "$(dirname "$FINAL_DEST")" || die "Impossible de créer le répertoire parent de destination"

DEST_PARENT_AVAIL_KB=$(df --output=avail "$(dirname "$FINAL_DEST")" | tail -1)
DEST_PARENT_AVAIL_GB=$(( DEST_PARENT_AVAIL_KB / 1024 / 1024 ))
DB_SIZE_GB=$(du -sBG "$DB" | awk '{gsub("G","",$1); print $1}')

log "Taille de la base construite: ${DB_SIZE_GB} GB. Espace disponible sur destination: ${DEST_PARENT_AVAIL_GB} GB."
(( DEST_PARENT_AVAIL_GB > DB_SIZE_GB )) || die "Espace insuffisant sur ${FINAL_DEST} pour copier la base (${DB_SIZE_GB} GB nécessaires)."

mkdir -p "$FINAL_DEST"
rsync -avh --progress "${DB}/" "${FINAL_DEST}/" 2>&1 | tail -n 50 | tee -a "$LOGFILE"

log "Copie terminée. Contenu final:"
ls -lh "$FINAL_DEST" | tee -a "$LOGFILE"

#-------------------------------------------------------------------------------
# 7. NETTOYAGE DE L'ESPACE DE TRAVAIL RAM/DISQUE
#-------------------------------------------------------------------------------
log "=== Étape 7: Nettoyage ==="
read -p "Supprimer le répertoire de travail temporaire ${WORKDIR} ? [y/N] " -n 1 -r REPLY
echo
if [[ "$REPLY" =~ ^[Yy]$ ]]; then
    rm -rf "$WORKDIR"
    log "Répertoire de travail supprimé."
else
    log "Répertoire de travail conservé: $WORKDIR (à supprimer manuellement pour libérer la RAM/disque)."
fi

log "=== TERMINÉ: base ${DB_NAME} disponible dans ${FINAL_DEST} ==="
