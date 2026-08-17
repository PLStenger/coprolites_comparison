#!/bin/bash
#SBATCH --job-name=09_kraken_goat_sheep
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites_comparison/00_scripts/09_kraken_goat_sheep.err"
#SBATCH --output="/home/plstenge/coprolites_comparison/00_scripts/09_kraken_goat_sheep.out"

# ==============================================================================
# Script : 09_kraken_following_goat_sheep_k35_29_25_mpa.sh
# But :
# 1. Regrouper les FASTQ merged apres fastp en deux pools :
#      goat  = cop410 + cop414 + cop417
#      sheep = cop408 + cop412
# 2. Lancer la cascade Kraken2 k35 -> k29 -> k25 sur chacun des deux pools.
# 3. Produire les tables MPA par niveau et les tables MPA finales :
#      - une table finale par pool ;
#      - une table combinant goat et sheep.
#
# Le script utilise exclusivement les reads merged. Les reads non-merges R1/R2
# ne sont volontairement pas inclus.
# ==============================================================================

module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics

# ==============================================================================
# PATHS ET PARAMETRES
# ==============================================================================

WORKDIR="/home/plstenge/coprolites_comparison"
ENTREE_FASTP="${WORKDIR}/06_fastp"
ENTREE_GROUPED="${WORKDIR}/06_fastp_grouped_species"

KRAKEN2_DB_K35="/storage/groups/gdec/shared/Kraken_database/k2_core_nt_20250609"
KRAKEN2_DB_K29="/storage/groups/gdec/shared/Kraken_database/core_nt_k29"
KRAKEN2_DB_K25="/storage/groups/gdec/shared/Kraken_database/core_nt_k25"

SORTIE="${WORKDIR}/09_kraken_following_goat_sheep"
SORTIE_K35="${SORTIE}/k35"
SORTIE_K29="${SORTIE}/k29"
SORTIE_K25="${SORTIE}/k25"
SORTIE_TMP="${SORTIE}/tmp_unassigned_fastq"
SORTIE_MPA="${SORTIE}/mpa"
SORTIE_MPA_K35="${SORTIE_MPA}/k35"
SORTIE_MPA_K29="${SORTIE_MPA}/k29"
SORTIE_MPA_K25="${SORTIE_MPA}/k25"
SORTIE_MPA_PER_SAMPLE="${SORTIE_MPA}/per_sample_combined"
SORTIE_COMBINED_KRAKEN="${SORTIE}/combined_kraken_for_mapdamage"

KRAKENTOOLS_DIR="${WORKDIR}/08_krakentools/KrakenTools"
THREADS="${SLURM_CPUS_PER_TASK:-36}"

SAMPLES=("goat" "sheep")

mkdir -p \
  "${ENTREE_GROUPED}" \
  "${SORTIE_K35}" \
  "${SORTIE_K29}" \
  "${SORTIE_K25}" \
  "${SORTIE_TMP}" \
  "${SORTIE_MPA_K35}" \
  "${SORTIE_MPA_K29}" \
  "${SORTIE_MPA_K25}" \
  "${SORTIE_MPA_PER_SAMPLE}" \
  "${SORTIE_COMBINED_KRAKEN}"

# ==============================================================================
# FONCTIONS
# ==============================================================================

require_nonempty_file() {
  local file="$1"
  if [[ ! -s "${file}" ]]; then
    echo "ERREUR : fichier absent ou vide : ${file}" >&2
    exit 1
  fi
}

extract_unassigned() {
  local kraken_file="$1"
  local seq_file="$2"
  local out_fastq="$3"

  python3 "${KRAKENTOOLS_DIR}/extract_kraken_reads.py" \
    -k "${kraken_file}" \
    -s "${seq_file}" \
    -t 0 \
    -o "${out_fastq}" \
    --fastq-output
}

convert_reports_to_mpa() {
  local in_dir="$1"
  local out_dir="$2"
  local label="$3"
  local report
  local base
  local mpa_file

  echo ""
  echo "Conversion des rapports en MPA : niveau ${label}"

  for report in "${in_dir}"/*.report; do
    [[ -s "${report}" ]] || continue
    base=$(basename "${report}" .report)
    mpa_file="${out_dir}/${base}.mpa"
    echo "  ${base}"
    python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" -r "${report}" -o "${mpa_file}"
  done
}

# ==============================================================================
# STEP 0 : VERIFICATION DE KRAKENTOOLS
# ==============================================================================

echo "======================================================================="
echo "STEP 0 : Verification de KrakenTools"
echo "======================================================================="

if [[ ! -d "${KRAKENTOOLS_DIR}" ]]; then
  echo "KrakenTools absent : installation dans ${KRAKENTOOLS_DIR}"
  mkdir -p "${WORKDIR}/08_krakentools"
  git clone https://github.com/jenniferlu717/KrakenTools.git "${KRAKENTOOLS_DIR}"
else
  echo "KrakenTools trouve : ${KRAKENTOOLS_DIR}"
fi

require_nonempty_file "${KRAKENTOOLS_DIR}/extract_kraken_reads.py"
require_nonempty_file "${KRAKENTOOLS_DIR}/kreport2mpa.py"
require_nonempty_file "${KRAKENTOOLS_DIR}/combine_mpa.py"

# ==============================================================================
# STEP 1 : CREATION DES POOLS GOAT ET SHEEP A PARTIR DES READS MERGED
# ==============================================================================

echo "======================================================================="
echo "STEP 1 : Regroupement des FASTQ merged"
echo "======================================================================="

COP410="${ENTREE_FASTP}/clean_cop410_grouped_dedup_fastp_merged.fastq.gz"
COP414="${ENTREE_FASTP}/clean_cop414_grouped_dedup_fastp_merged.fastq.gz"
COP417="${ENTREE_FASTP}/clean_cop417_grouped_dedup_fastp_merged.fastq.gz"
COP408="${ENTREE_FASTP}/clean_cop408_grouped_dedup_fastp_merged.fastq.gz"
COP412="${ENTREE_FASTP}/clean_cop412_grouped_dedup_fastp_merged.fastq.gz"

for f in "${COP410}" "${COP414}" "${COP417}" "${COP408}" "${COP412}"; do
  require_nonempty_file "${f}"
done

GOAT_FASTQ="${ENTREE_GROUPED}/clean_goat_grouped_dedup_fastp_merged.fastq.gz"
SHEEP_FASTQ="${ENTREE_GROUPED}/clean_sheep_grouped_dedup_fastp_merged.fastq.gz"

# Les fichiers gzip concatenes sont des flux gzip valides.
# Les fichiers de pools sont recrees a chaque lancement.
echo "Creation du pool goat : cop410 + cop414 + cop417"
cat "${COP410}" "${COP414}" "${COP417}" > "${GOAT_FASTQ}"

echo "Creation du pool sheep : cop408 + cop412"
cat "${COP408}" "${COP412}" > "${SHEEP_FASTQ}"

for f in "${GOAT_FASTQ}" "${SHEEP_FASTQ}"; do
  gzip -t "${f}"
  n_reads=$(gzip -cd "${f}" | awk 'END {print NR / 4}')
  echo "$(basename "${f}") : ${n_reads} reads merged"
done

# ==============================================================================
# STEP 2 : CASCADE KRAKEN2 K35 -> K29 -> K25
# ==============================================================================

echo "======================================================================="
echo "STEP 2 : Cascade Kraken2 k35 -> k29 -> k25"
echo "Pools : ${SAMPLES[*]}"
echo "======================================================================="

for sample in "${SAMPLES[@]}"; do
  echo ""
  echo "-----------------------------------------------------------------------"
  echo "Pool : ${sample}"
  echo "-----------------------------------------------------------------------"

  merged_fastq="${ENTREE_GROUPED}/clean_${sample}_grouped_dedup_fastp_merged.fastq.gz"
  require_nonempty_file "${merged_fastq}"

  k35_kraken="${SORTIE_K35}/${sample}_merged.kraken"
  k35_report="${SORTIE_K35}/${sample}_merged.report"

  k29_kraken="${SORTIE_K29}/${sample}_merged_unassigned_k35.kraken"
  k29_report="${SORTIE_K29}/${sample}_merged_unassigned_k35.report"
  unassigned_k35_fastq="${SORTIE_TMP}/${sample}_unassigned_after_k35.fastq"

  k25_kraken="${SORTIE_K25}/${sample}_merged_unassigned_k29.kraken"
  k25_report="${SORTIE_K25}/${sample}_merged_unassigned_k29.report"
  unassigned_k29_fastq="${SORTIE_TMP}/${sample}_unassigned_after_k29.fastq"

  combined_kraken="${SORTIE_COMBINED_KRAKEN}/${sample}_following.kraken"

  rm -f \
    "${k35_kraken}" "${k35_report}" \
    "${k29_kraken}" "${k29_report}" \
    "${k25_kraken}" "${k25_report}" \
    "${unassigned_k35_fastq}" "${unassigned_k29_fastq}" \
    "${combined_kraken}"

  echo "[k35] Classification du pool complet"
  kraken2 \
    --conf 0.2 \
    --db "${KRAKEN2_DB_K35}" \
    --threads "${THREADS}" \
    --output "${k35_kraken}" \
    --report "${k35_report}" \
    "${merged_fastq}"

  n_unclass_k35=$(awk -F'\t' '$1 == "U" {n++} END {print n + 0}' "${k35_kraken}")
  echo "[k35] Reads non assignes : ${n_unclass_k35}"

  if [[ "${n_unclass_k35}" -gt 0 ]]; then
    echo "[k29] Extraction des reads non assignes par k35"
    extract_unassigned "${k35_kraken}" "${merged_fastq}" "${unassigned_k35_fastq}"

    if [[ -s "${unassigned_k35_fastq}" ]]; then
      echo "[k29] Classification des reads non assignes par k35"
      kraken2 \
        --conf 0.2 \
        --db "${KRAKEN2_DB_K29}" \
        --threads "${THREADS}" \
        --output "${k29_kraken}" \
        --report "${k29_report}" \
        "${unassigned_k35_fastq}"
    else
      : > "${k29_kraken}"
      : > "${k29_report}"
    fi
  else
    : > "${k29_kraken}"
    : > "${k29_report}"
  fi

  n_unclass_k29=0
  if [[ -s "${k29_kraken}" ]]; then
    n_unclass_k29=$(awk -F'\t' '$1 == "U" {n++} END {print n + 0}' "${k29_kraken}")
  fi
  echo "[k29] Reads non assignes : ${n_unclass_k29}"

  if [[ "${n_unclass_k29}" -gt 0 ]]; then
    echo "[k25] Extraction des reads non assignes par k29"
    extract_unassigned "${k29_kraken}" "${unassigned_k35_fastq}" "${unassigned_k29_fastq}"

    if [[ -s "${unassigned_k29_fastq}" ]]; then
      echo "[k25] Classification des reads non assignes par k35 et k29"
      kraken2 \
        --conf 0.2 \
        --db "${KRAKEN2_DB_K25}" \
        --threads "${THREADS}" \
        --output "${k25_kraken}" \
        --report "${k25_report}" \
        "${unassigned_k29_fastq}"
    else
      : > "${k25_kraken}"
      : > "${k25_report}"
    fi
  else
    : > "${k25_kraken}"
    : > "${k25_report}"
  fi

  n_unclass_k25=0
  if [[ -s "${k25_kraken}" ]]; then
    n_unclass_k25=$(awk -F'\t' '$1 == "U" {n++} END {print n + 0}' "${k25_kraken}")
  fi
  echo "[k25] Reads jamais assignes : ${n_unclass_k25}"

  cat "${k35_kraken}" "${k29_kraken}" "${k25_kraken}" > "${combined_kraken}"
  echo "Fichier Kraken combine : ${combined_kraken}"

  rm -f "${unassigned_k35_fastq}" "${unassigned_k29_fastq}"
done

# ==============================================================================
# STEP 3 : CONVERSION DES RAPPORTS KRAKEN VERS MPA
# ==============================================================================

echo "======================================================================="
echo "STEP 3 : Conversion des rapports Kraken en tables MPA"
echo "======================================================================="

convert_reports_to_mpa "${SORTIE_K35}" "${SORTIE_MPA_K35}" "k35"
convert_reports_to_mpa "${SORTIE_K29}" "${SORTIE_MPA_K29}" "k29"
convert_reports_to_mpa "${SORTIE_K25}" "${SORTIE_MPA_K25}" "k25"

# ==============================================================================
# STEP 4 : FUSION K35 + K29 + K25 PAR POOL
# ==============================================================================

echo "======================================================================="
echo "STEP 4 : Fusion des comptages MPA k35 + k29 + k25 par pool"
echo "======================================================================="

SORTIE_MPA_K35_ENV="${SORTIE_MPA_K35}" \
SORTIE_MPA_K29_ENV="${SORTIE_MPA_K29}" \
SORTIE_MPA_K25_ENV="${SORTIE_MPA_K25}" \
SORTIE_MPA_PER_SAMPLE_ENV="${SORTIE_MPA_PER_SAMPLE}" \
python3 << 'PYEOF'
import os

samples = ["goat", "sheep"]
mpa_dirs = {
    "k35": os.environ["SORTIE_MPA_K35_ENV"],
    "k29": os.environ["SORTIE_MPA_K29_ENV"],
    "k25": os.environ["SORTIE_MPA_K25_ENV"],
}
out_dir = os.environ["SORTIE_MPA_PER_SAMPLE_ENV"]


def load_mpa(path):
    counts = {}
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return counts
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            taxon = fields[0]
            try:
                value = float(fields[1])
            except ValueError:
                continue
            counts[taxon] = counts.get(taxon, 0.0) + value
    return counts


for sample in samples:
    filenames = {
        "k35": f"{sample}_merged.mpa",
        "k29": f"{sample}_merged_unassigned_k35.mpa",
        "k25": f"{sample}_merged_unassigned_k29.mpa",
    }
    combined = {}
    found = False

    for level, directory in mpa_dirs.items():
        path = os.path.join(directory, filenames[level])
        if os.path.isfile(path) and os.path.getsize(path) > 0:
            found = True
        for taxon, value in load_mpa(path).items():
            combined[taxon] = combined.get(taxon, 0.0) + value

    if not found:
        raise RuntimeError(f"Aucune table MPA non vide trouvee pour {sample}")

    out_path = os.path.join(out_dir, f"{sample}_following.mpa")
    with open(out_path, "w", encoding="utf-8") as out:
        for taxon in sorted(combined):
            value = combined[taxon]
            formatted = str(int(value)) if value.is_integer() else str(value)
            out.write(f"{taxon}\t{formatted}\n")
    print(f"Table MPA finale ecrite : {out_path}")
PYEOF

# ==============================================================================
# STEP 5 : TABLES MPA COMBINEES
# ==============================================================================

echo "======================================================================="
echo "STEP 5 : Tables MPA combinees goat + sheep"
echo "======================================================================="

final_mpa_files=()
for sample in "${SAMPLES[@]}"; do
  file="${SORTIE_MPA_PER_SAMPLE}/${sample}_following.mpa"
  [[ -s "${file}" ]] && final_mpa_files+=("${file}")
done

if [[ ${#final_mpa_files[@]} -ne ${#SAMPLES[@]} ]]; then
  echo "ERREUR : au moins une table MPA finale est manquante." >&2
  exit 1
fi

python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" \
  -i "${final_mpa_files[@]}" \
  -o "${SORTIE_MPA}/combined_goat_sheep_following.tsv"

echo "Table finale goat + sheep : ${SORTIE_MPA}/combined_goat_sheep_following.tsv"

for level_dir_label in \
  "${SORTIE_MPA_K35}:k35" \
  "${SORTIE_MPA_K29}:k29" \
  "${SORTIE_MPA_K25}:k25"; do

  IFS=':' read -r level_dir label <<< "${level_dir_label}"
  mpa_files=("${level_dir}"/*.mpa)

  if [[ ${#mpa_files[@]} -gt 0 ]]; then
    python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" \
      -i "${mpa_files[@]}" \
      -o "${level_dir}/combined_goat_sheep.tsv"
    echo "Table combinee ${label} : ${level_dir}/combined_goat_sheep.tsv"
  fi
done

# ==============================================================================
# RESUME
# ==============================================================================

echo ""
echo "======================================================================="
echo "TERMINE : Kraken following goat / sheep"
echo "======================================================================="
echo "FASTQ pools : ${ENTREE_GROUPED}"
echo "Kraken k35 : ${SORTIE_K35}"
echo "Kraken k29 : ${SORTIE_K29}"
echo "Kraken k25 : ${SORTIE_K25}"
echo "MPA par niveau : ${SORTIE_MPA}/{k35,k29,k25}"
echo "MPA finale goat : ${SORTIE_MPA_PER_SAMPLE}/goat_following.mpa"
echo "MPA finale sheep : ${SORTIE_MPA_PER_SAMPLE}/sheep_following.mpa"
echo "MPA combinee : ${SORTIE_MPA}/combined_goat_sheep_following.tsv"
echo "======================================================================="
