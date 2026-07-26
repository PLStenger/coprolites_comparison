#!/bin/bash
#SBATCH --job-name=build_nt_k29_PL
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=36
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=/home/plstenge/coprolites_comparison/00_scripts/build_nt_k29_PL.err
#SBATCH --output=/home/plstenge/coprolites_comparison/00_scripts/build_nt_k29_PL.out

#===============================================================================
# build_nt_k29_PL.sh
# Construction d'une base Kraken2 k=29 + fichiers Bracken
# Source FASTA: /storage/biodatabanks/ncbi/NT/current/fasta/All/all.fasta
# Destination finale: /storage/groups/gdec/shared/Kraken_database/nt_k29_PL
#===============================================================================

#-------------------------------------------------------------------------------
# CONFIGURATION
#-------------------------------------------------------------------------------
DB_NAME="nt_k29_PL"
FINAL_DEST="/storage/groups/gdec/shared/Kraken_database/${DB_NAME}"
SOURCE_FASTA="/storage/biodatabanks/ncbi/NT/current/fasta/All/all.fasta"

KMER_LEN=29
MINIMIZER_LEN=23
MINIMIZER_SPACES=2

THREADS="${SLURM_CPUS_PER_TASK:-36}"
MAX_DB_SIZE_GB=800

WORKDIR_BASE="/dev/shm/${DB_NAME}_build"
FALLBACK_WORKDIR_BASE="/tmp/${DB_NAME}_build"

CONDA_ENV="metagenomics"
BRACKEN_ENV="metagenomics"
BRACKEN_READ_LENS=(50 75 100 150)

SCRIPT_LOG_DIR="/home/plstenge/coprolites_comparison/00_scripts"
mkdir -p "${SCRIPT_LOG_DIR}"
LOGFILE="${SCRIPT_LOG_DIR}/build_${DB_NAME}_$(date +%Y%m%d_%H%M%S).log"

export LC_ALL=C
export LANG=C
export OMP_NUM_THREADS="${THREADS}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

#-------------------------------------------------------------------------------
# OUTILS
#-------------------------------------------------------------------------------
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"
}

die() {
    log "ERREUR: $*"
    exit 1
}

check_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Commande requise introuvable: $1"
}

human_gb_from_kb() {
    local kb="$1"
    echo $(( kb / 1024 / 1024 ))
}

cleanup_on_error() {
    log "Échec à la ligne ${1}. Voir ${LOGFILE}"
}
trap 'cleanup_on_error $LINENO' ERR

#-------------------------------------------------------------------------------
# ENV CONDA
#-------------------------------------------------------------------------------
module load conda/4.12.0
source ~/.bashrc
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}" || die "Impossible d'activer l'environnement ${CONDA_ENV}"

#-------------------------------------------------------------------------------
# VÉRIFICATIONS
#-------------------------------------------------------------------------------
log "=== Étape 0: Vérifications préalables ==="

check_cmd awk
check_cmd grep
check_cmd sed
check_cmd cut
check_cmd sort
check_cmd uniq
check_cmd wc
check_cmd df
check_cmd stat
check_cmd tar
check_cmd gzip
check_cmd wget
check_cmd rsync
check_cmd numfmt
check_cmd conda
check_cmd k2
check_cmd python

[[ -f "${SOURCE_FASTA}" ]] || die "FASTA source introuvable: ${SOURCE_FASTA}"

log "Conda env actif: ${CONDA_DEFAULT_ENV:-NA}"
log "Executable k2: $(command -v k2)"
log "Python: $(python --version 2>&1)"
log "SLURM job id: ${SLURM_JOB_ID:-NA}"
log "SLURM cpus-per-task: ${THREADS}"
log "SLURM mem per node: ${SLURM_MEM_PER_NODE:-NA}"

MEM_ESTIMATE_GB=341
MEM_SAFETY_GB=443

AVAIL_MEM_KB=$(awk '/MemAvailable/ {print $2}' /proc/meminfo)
AVAIL_MEM_GB=$(human_gb_from_kb "${AVAIL_MEM_KB}")

log "Mémoire disponible détectée: ${AVAIL_MEM_GB} GB"
log "Estimation mémoire nécessaire pour k=${KMER_LEN}: ~${MEM_ESTIMATE_GB} GB ; marge recommandée: ~${MEM_SAFETY_GB} GB"

if (( AVAIL_MEM_GB < MEM_ESTIMATE_GB )); then
    die "Mémoire disponible insuffisante (${AVAIL_MEM_GB} GB < ${MEM_ESTIMATE_GB} GB)"
fi

SHM_AVAIL_KB=$(df --output=avail /dev/shm 2>/dev/null | tail -1 || echo 0)
SHM_AVAIL_GB=$(human_gb_from_kb "${SHM_AVAIL_KB}")
SOURCE_SIZE_BYTES=$(stat -c%s "${SOURCE_FASTA}")
SOURCE_SIZE_GB=$(( SOURCE_SIZE_BYTES / 1024 / 1024 / 1024 ))

log "Taille du FASTA source: ${SOURCE_SIZE_GB} GB"
log "Espace disponible /dev/shm: ${SHM_AVAIL_GB} GB"

NEEDED_SHM_GB=$(( MEM_ESTIMATE_GB + 30 ))

WORKDIR="${WORKDIR_BASE}"
if (( SHM_AVAIL_GB >= NEEDED_SHM_GB )); then
    log "Utilisation de /dev/shm comme répertoire de travail"
else
    WORKDIR="${FALLBACK_WORKDIR_BASE}"
    log "Bascule sur ${WORKDIR} car /dev/shm insuffisant"
fi

mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

DB="${WORKDIR}/${DB_NAME}"
mkdir -p "${DB}/library/added"
mkdir -p "${DB}/taxonomy"

log "Répertoire de travail: ${WORKDIR}"
log "Répertoire DB: ${DB}"
log "TMPDIR: ${TMPDIR}"

#-------------------------------------------------------------------------------
# TÉLÉCHARGEMENT TAXONOMIE NCBI
#-------------------------------------------------------------------------------
log "=== Étape 1: Téléchargement robuste de la taxonomie NCBI ==="

NCBI_TAXONOMY_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
NCBI_GB_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"
NCBI_WGS_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz"

download_file_retry() {
    local url="$1"
    local outfile="$2"
    local max_tries=5
    local attempt=1
    local sleep_time=60

    while (( attempt <= max_tries )); do
        log "Téléchargement ${outfile} - tentative ${attempt}/${max_tries}"
        if wget -4 --tries=1 --timeout=120 --waitretry=5 -O "${outfile}" "${url}" >> "${LOGFILE}" 2>&1; then
            [[ -s "${outfile}" ]] && return 0
        fi
        rm -f "${outfile}"
        log "Échec téléchargement ${outfile}"
        if (( attempt < max_tries )); then
            log "Nouvelle tentative dans ${sleep_time}s"
            sleep "${sleep_time}"
        fi
        attempt=$(( attempt + 1 ))
    done
    return 1
}

if [[ -s "${DB}/taxonomy/nodes.dmp" && -s "${DB}/taxonomy/names.dmp" && -s "${DB}/taxonomy/nucl_gb.accession2taxid" && -s "${DB}/taxonomy/nucl_wgs.accession2taxid" ]]; then
    log "Taxonomie déjà présente, skip"
else
    rm -rf "${DB}/taxonomy"
    mkdir -p "${DB}/taxonomy"
    cd "${DB}/taxonomy"

    download_file_retry "${NCBI_TAXONOMY_URL}" "taxdump.tar.gz" || die "Impossible de télécharger taxdump.tar.gz"
    download_file_retry "${NCBI_GB_URL}" "nucl_gb.accession2taxid.gz" || die "Impossible de télécharger nucl_gb.accession2taxid.gz"
    download_file_retry "${NCBI_WGS_URL}" "nucl_wgs.accession2taxid.gz" || die "Impossible de télécharger nucl_wgs.accession2taxid.gz"

    log "Extraction de taxdump.tar.gz"
    tar -xzf taxdump.tar.gz >> "${LOGFILE}" 2>&1

    log "Décompression des fichiers accession2taxid"
    gzip -dc nucl_gb.accession2taxid.gz > nucl_gb.accession2taxid
    gzip -dc nucl_wgs.accession2taxid.gz > nucl_wgs.accession2taxid

    [[ -s "nodes.dmp" ]] || die "nodes.dmp manquant après extraction"
    [[ -s "names.dmp" ]] || die "names.dmp manquant après extraction"
    [[ -s "nucl_gb.accession2taxid" ]] || die "nucl_gb.accession2taxid manquant"
    [[ -s "nucl_wgs.accession2taxid" ]] || die "nucl_wgs.accession2taxid manquant"

    cd - >/dev/null
    log "Taxonomie NCBI téléchargée avec succès"
fi

#-------------------------------------------------------------------------------
# GÉNÉRATION DU FICHIER seqid2taxid.map
#-------------------------------------------------------------------------------
log "=== Étape 2: Génération de seqid2taxid.map ==="

MAP_FILE="${DB}/seqid2taxid.map"

generate_map() {
    local gb="$1"
    local wgs="$2"
    local fasta="$3"
    local out="$4"

    awk '
    BEGIN { FS="\t"; OFS="\t"; }
    ARGIND <= 2 {
        if ($2 != "accession.version") {
            map[$2] = $3;
        }
        next;
    }
    /^>/ {
        full_id = $0;
        sub(/^>/, "", full_id);
        split(full_id, parts, /[ |]/);
        id = parts[1];
        if (id in map) {
            print id, map[id];
        } else {
            print id, "1";
        }
    }' "${gb}" "${wgs}" "${fasta}" > "${out}"
}

if [[ -s "${MAP_FILE}" ]]; then
    log "seqid2taxid.map déjà présent, skip"
else
    generate_map \
        "${DB}/taxonomy/nucl_gb.accession2taxid" \
        "${DB}/taxonomy/nucl_wgs.accession2taxid" \
        "${SOURCE_FASTA}" \
        "${MAP_FILE}"
fi

N_SEQ=$(grep -c '^>' "${SOURCE_FASTA}")
N_MAP=$(wc -l < "${MAP_FILE}")

log "Nombre de séquences dans le FASTA source: ${N_SEQ}"
log "Nombre de lignes dans seqid2taxid.map: ${N_MAP}"

if (( N_SEQ != N_MAP )); then
    die "Incohérence FASTA/mapping: ${N_SEQ} séquences vs ${N_MAP} lignes"
fi

N_ROOT=$(awk '$2 == 1 {c++} END {print c+0}' "${MAP_FILE}")
ROOT_PCT=$(awk -v nroot="${N_ROOT}" -v nseq="${N_SEQ}" 'BEGIN{ if(nseq==0){print "0.00"} else {printf "%.2f", (100*nroot/nseq)} }')

log "Nombre d’entrées taxid=1: ${N_ROOT} (${ROOT_PCT}%)"

#-------------------------------------------------------------------------------
# PRÉPARATION DU FASTA
#-------------------------------------------------------------------------------
log "=== Étape 3: Préparation du FASTA .fna ==="

LINK_PATH="${DB}/library/added/file.fna"
rm -f "${LINK_PATH}"
ln -s "${SOURCE_FASTA}" "${LINK_PATH}"
echo "file.fna" > "${DB}/library/added/manifest.txt"

[[ -L "${LINK_PATH}" ]] || die "Le lien symbolique .fna n'a pas été créé"

#-------------------------------------------------------------------------------
# BUILD KRAKEN2
#-------------------------------------------------------------------------------
log "=== Étape 4: Build Kraken2 k=${KMER_LEN} ==="
log "Threads: ${THREADS}"
log "max-db-size: ${MAX_DB_SIZE_GB} GB"

START_TIME=$(date +%s)

set +e
k2 build \
    --threads "${THREADS}" \
    --db "${DB}" \
    --kmer-len "${KMER_LEN}" \
    --minimizer-len "${MINIMIZER_LEN}" \
    --minimizer-spaces "${MINIMIZER_SPACES}" \
    --max-db-size "${MAX_DB_SIZE_GB}" \
    >> "${LOGFILE}" 2>&1
BUILD_EXIT_CODE=$?
set -e

END_TIME=$(date +%s)
BUILD_HOURS=$(( (END_TIME - START_TIME) / 3600 ))
BUILD_MINS=$(( ((END_TIME - START_TIME) % 3600) / 60 ))

log "Durée du build: ${BUILD_HOURS}h ${BUILD_MINS}min"
log "Code retour k2 build: ${BUILD_EXIT_CODE}"

for f in hash.k2d opts.k2d taxo.k2d; do
    [[ -s "${DB}/${f}" ]] || die "Fichier manquant après build: ${DB}/${f}"
done

log "Build Kraken2 terminé avec succès"

#-------------------------------------------------------------------------------
# BRACKEN
#-------------------------------------------------------------------------------
log "=== Étape 5: Génération des fichiers Bracken ==="

if ! command -v bracken-build >/dev/null 2>&1; then
    conda deactivate || true
    conda activate "${BRACKEN_ENV}" || die "Impossible d'activer ${BRACKEN_ENV} pour bracken-build"
fi

check_cmd bracken-build

for RLEN in "${BRACKEN_READ_LENS[@]}"; do
    log "bracken-build pour read length=${RLEN}"
    bracken-build \
        -d "${DB}" \
        -t "${THREADS}" \
        -k "${KMER_LEN}" \
        -l "${RLEN}" \
        >> "${LOGFILE}" 2>&1
done

conda deactivate || true
conda activate "${CONDA_ENV}" || true

#-------------------------------------------------------------------------------
# COPIE FINALE
#-------------------------------------------------------------------------------
log "=== Étape 6: Copie vers la destination finale ==="

mkdir -p "$(dirname "${FINAL_DEST}")"
mkdir -p "${FINAL_DEST}"

DEST_PARENT_AVAIL_KB=$(df --output=avail "$(dirname "${FINAL_DEST}")" | tail -1)
DEST_PARENT_AVAIL_GB=$(human_gb_from_kb "${DEST_PARENT_AVAIL_KB}")
DB_SIZE_GB=$(du -sBG "${DB}" | awk '{gsub("G","",$1); print $1}')

log "Taille estimée de la base: ${DB_SIZE_GB} GB"
log "Espace disponible sur destination: ${DEST_PARENT_AVAIL_GB} GB"

(( DEST_PARENT_AVAIL_GB > DB_SIZE_GB )) || die "Espace insuffisant sur destination"

rsync -aH --info=progress2 "${DB}/" "${FINAL_DEST}/" >> "${LOGFILE}" 2>&1

log "Copie terminée"
ls -lh "${FINAL_DEST}" >> "${LOGFILE}" 2>&1

#-------------------------------------------------------------------------------
# NETTOYAGE
#-------------------------------------------------------------------------------
log "=== Étape 7: Nettoyage ==="

if [[ "${WORKDIR}" == /dev/shm/* || "${WORKDIR}" == /tmp/* ]]; then
    rm -rf "${WORKDIR}"
    log "Répertoire de travail supprimé: ${WORKDIR}"
else
    log "WORKDIR non supprimé par sécurité: ${WORKDIR}"
fi

log "=== TERMINÉ: base ${DB_NAME} disponible dans ${FINAL_DEST} ==="
