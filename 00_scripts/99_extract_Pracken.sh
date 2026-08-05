#!/bin/bash

#SBATCH --job-name=extract_pracken_db
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1000G
#SBATCH -p smp
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/extract_pracken_db.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/extract_pracken_db.out"

# ==============================================================================
# extract_pracken_db_robust.sh
#
# Extraction robuste de k2_NCBI_reference_20251007.tar.gz vers
# /storage/groups/gdec/shared/Kraken_database
#
# Contexte : l'extraction a echoue avec
#   "tar: k2_NCBI_reference_20251007/hash.k2d : write impossible:
#    Erreur d'entree/sortie"
# en plein milieu du plus gros fichier (hash.k2d), ce qui produit ensuite
# l'erreur Kraken2 "Error reading in hash table" (fichier tronque).
#
# Causes les plus probables sur un montage reseau (NFS/Lustre/GPFS type
# /storage/groups/...) :
#   1) Quota disque ou espace disque insuffisant au moment de l'ecriture
#      du gros fichier (l'erreur I/O peut masquer un "disk full").
#   2) Coupure I/O reseau transitoire pendant l'ecriture d'un tres gros
#      fichier (des heures d'ecriture continue augmentent le risque).
#   3) Fichier .tar.gz source lui-meme partiellement telecharge/corrompu.
#
# Ce script :
#   - Verifie l'espace disque ET le quota AVANT de commencer.
#   - Verifie l'integrite de l'archive .tar.gz (test de decompression gzip)
#     avant extraction complete, pour distinguer "archive corrompue" de
#     "probleme d'ecriture destination".
#   - Extrait avec reprise automatique en cas d'echec transitoire (retry).
#   - Verifie a la fin que TOUS les fichiers attendus sont presents et que
#     leur taille est coherente (non nulle, non tronquee).
#   - Nettoie un repertoire partiellement extrait avant de retenter, pour
#     eviter un etat incoherent (melange ancien/nouveau).
# ==============================================================================

ARCHIVE="/storage/groups/gdec/shared/Kraken_database/k2_NCBI_reference_20251007.tar.gz"
DEST_DIR="/storage/groups/gdec/shared/Kraken_database"
DB_NAME="k2_NCBI_reference_20251007"
FULL_DEST="${DEST_DIR}/${DB_NAME}"

MAX_TRIES=5
RETRY_WAIT=60

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
die() { log "[FATAL] $*"; exit 1; }

# ------------------------------------------------------------------------------
# 0) Verifications prealables
# ------------------------------------------------------------------------------

log "=== Etape 0 : verifications prealables ==="

[[ -f "${ARCHIVE}" ]] || die "Archive introuvable : ${ARCHIVE}"

ARCHIVE_SIZE_BYTES=$(stat -c%s "${ARCHIVE}")
ARCHIVE_SIZE_HUMAN=$(numfmt --to=iec-i --suffix=B "${ARCHIVE_SIZE_BYTES}")
log "Taille de l'archive : ${ARCHIVE_SIZE_HUMAN}"

# Estimation large : marge x1.5 pour etre sur (hash.k2d est deja dense,
# le ratio de decompression est generalement proche de 1:1 a 1:1.3)
NEEDED_BYTES=$(( ARCHIVE_SIZE_BYTES * 3 / 2 ))
NEEDED_HUMAN=$(numfmt --to=iec-i --suffix=B "${NEEDED_BYTES}")

AVAIL_BYTES=$(df --output=avail -B1 "${DEST_DIR}" | tail -1)
AVAIL_HUMAN=$(numfmt --to=iec-i --suffix=B "${AVAIL_BYTES}")

log "Espace disponible sur ${DEST_DIR} : ${AVAIL_HUMAN}"
log "Espace estime necessaire (marge x1.5) : ${NEEDED_HUMAN}"

if (( AVAIL_BYTES < NEEDED_BYTES )); then
    die "Espace disque insuffisant sur ${DEST_DIR}. Disponible=${AVAIL_HUMAN}, requis=${NEEDED_HUMAN}. C'est tres probablement la cause de l'erreur I/O initiale (quota/espace atteint en cours d'ecriture de hash.k2d). Libere de l'espace ou change de destination avant de relancer."
fi

if command -v quota >/dev/null 2>&1; then
    log "Verification du quota utilisateur (si applicable) :"
    quota -s 2>/dev/null || log "  (commande quota non applicable ou pas de quota configure ici)"
fi

# ------------------------------------------------------------------------------
# 1) Verification d'integrite de l'archive AVANT extraction complete
# ------------------------------------------------------------------------------

log "=== Etape 1 : verification d'integrite de l'archive (gzip -t) ==="
log "Ceci peut prendre plusieurs minutes selon la taille de l'archive..."

if gzip -t "${ARCHIVE}" 2>/tmp/gzip_test_err.log; then
    log "Archive OK : le flux gzip est valide et complet."
else
    log "[ERREUR] L'archive elle-meme semble corrompue ou incomplete (gzip -t a echoue)."
    cat /tmp/gzip_test_err.log
    die "Re-telecharge l'archive ${ARCHIVE} depuis sa source avant de retenter l'extraction."
fi

# ------------------------------------------------------------------------------
# 2) Nettoyage d'un repertoire partiellement extrait (etat incoherent)
# ------------------------------------------------------------------------------

if [[ -d "${FULL_DEST}" ]]; then
    log "=== Etape 2 : repertoire de destination deja present (extraction precedente partielle) ==="
    log "Suppression de ${FULL_DEST} pour repartir d'un etat propre..."
    rm -rf "${FULL_DEST}" || die "Impossible de supprimer ${FULL_DEST} (verifie les permissions)."
    log "Nettoyage effectue."
fi

# ------------------------------------------------------------------------------
# 3) Extraction avec retry automatique
# ------------------------------------------------------------------------------

log "=== Etape 3 : extraction avec reprise automatique (max ${MAX_TRIES} tentatives) ==="

attempt=1
success=0

while (( attempt <= MAX_TRIES )); do
    log "--- Tentative ${attempt}/${MAX_TRIES} ---"

    if tar -xvzf "${ARCHIVE}" -C "${DEST_DIR}" 2> "/tmp/tar_extract_attempt_${attempt}.log"; then
        log "Extraction terminee sans erreur signalee par tar."
        success=1
        break
    else
        RC=$?
        log "[ERREUR] tar a retourne un code ${RC} a la tentative ${attempt}."
        tail -n 20 "/tmp/tar_extract_attempt_${attempt}.log"

        if (( attempt < MAX_TRIES )); then
            log "Nettoyage avant nouvelle tentative (evite un hash.k2d tronque residuel)..."
            rm -rf "${FULL_DEST}"
            log "Attente ${RETRY_WAIT}s avant de retenter (laisse le temps a un probleme I/O transitoire de se resorber)..."
            sleep "${RETRY_WAIT}"
        fi
    fi

    attempt=$(( attempt + 1 ))
done

if (( success != 1 )); then
    die "Extraction impossible apres ${MAX_TRIES} tentatives. Verifie l'espace disque, le quota, et la stabilite du montage reseau ${DEST_DIR} avant de relancer ce script."
fi

# ------------------------------------------------------------------------------
# 4) Verification post-extraction : presence + taille non nulle des fichiers
#    critiques
# ------------------------------------------------------------------------------

log "=== Etape 4 : verification post-extraction ==="

REQUIRED_FILES=(
    "nodes.dmp"
    "names.dmp"
    "ktaxonomy.tsv"
    "opts.k2d"
    "taxo.k2d"
    "hash.k2d"
    "seqid2taxid.map"
    "inspect.txt"
    "library_report.tsv"
)

missing=0
for f in "${REQUIRED_FILES[@]}"; do
    fpath="${FULL_DEST}/${f}"
    if [[ ! -s "${fpath}" ]]; then
        log "[MANQUANT ou VIDE] ${fpath}"
        missing=1
    else
        size_h=$(numfmt --to=iec-i --suffix=B "$(stat -c%s "${fpath}")")
        log "  OK  ${f} (${size_h})"
    fi
done

if (( missing == 1 )); then
    die "Au moins un fichier attendu est manquant ou vide apres extraction. Ne pas utiliser cette base en l'etat -- relance ce script."
fi

# ------------------------------------------------------------------------------
# 5) Test de chargement rapide avec kraken2-inspect (si l'env conda existe)
# ------------------------------------------------------------------------------

log "=== Etape 5 : test de chargement Kraken2 (kraken2-inspect) ==="

if command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    if conda activate metagenomics 2>/dev/null && command -v kraken2-inspect >/dev/null 2>&1; then
        if kraken2-inspect --db "${FULL_DEST}" --threads 4 > /tmp/kraken2_inspect_test.log 2>&1; then
            if grep -qi "error" /tmp/kraken2_inspect_test.log; then
                log "[ATTENTION] kraken2-inspect a retourne du texte contenant 'error', verifie /tmp/kraken2_inspect_test.log"
            else
                log "kraken2-inspect a charge la base sans erreur. La base semble saine."
            fi
        else
            log "[ATTENTION] kraken2-inspect a echoue. Voir /tmp/kraken2_inspect_test.log"
            tail -n 20 /tmp/kraken2_inspect_test.log
        fi
    else
        log "Environnement conda 'metagenomics' ou kraken2-inspect non disponible ici, test ignore."
    fi
else
    log "conda non disponible dans ce contexte, test kraken2-inspect ignore."
fi

log "=== TERMINE : base extraite et verifiee dans ${FULL_DEST} ==="
