#!/usr/bin/env bash
# RNA-seq pipeline (mouse mm39), Hao's environment
# - Binaries:   /home/h392x566/RNA_Seq/bin/{STAR,samtools,featureCounts,trimmomatic-0.40.jar}
# - References: /home/h392x566/JZ/ref/{mm39.fa, mm39.ncbiRefSeq.gtf, STAR_index_mm39}
# Usage:
#   1) Put sample IDs in samples.txt (e.g., 16_S5, 17_S6, ...)
#   2) FASTQs named like: <SAMPLE>_L001_R1_001.fastq(.gz) and ..._R2_...
#   3) ./run_rnaseq_mm39.sh
set -euo pipefail

############################################
# 0) CONFIG
############################################
# Threads
THREADS=32

# Binaries / tools
BIN_DIR="/home/h392x566/RNA_Seq/bin"
STAR="${BIN_DIR}/STAR"
SAMTOOLS="${BIN_DIR}/samtools"
FEATURECOUNTS="${BIN_DIR}/featureCounts"
TRIM_JAR="${BIN_DIR}/trimmomatic-0.40.jar"

# References (mm39)
GENOME_FASTA="/home/h392x566/JZ/ref/mm39.fa"
GTF="/home/h392x566/JZ/ref/mm39.ncbiRefSeq.gtf"
STAR_INDEX="/home/h392x566/JZ/ref/STAR_index_mm39"

# Read length minus 1 for STAR splice junctions (set to 149 for 150bp reads; adjust if needed)
SJDB_OVERHANG=149

# Strandedness for featureCounts: 0=unstranded, 1=forward, 2=reverse (TruSeq Stranded is usually 2)
STRAND=2

# Annotation fields for featureCounts
FC_GTF_FEATURE="exon"
FC_GTF_GROUP="gene_id"   # change to gene_name if you prefer symbols

# Adapter file for Trimmomatic (update if needed)
ADAPTERS="./adapters/TruSeq3-PE.fa"

# Lanes to merge (add more like L003, L004 if present)
LANES=("L001" "L002")

############################################
# 1) PREP
############################################
mkdir -p 00_merged 01_trimmed 02_star 03_counts qc

echo "[INFO] Using tools:"
echo "  STAR:          $(${STAR} --version 2>&1 | head -n1 || echo 'OK')"
echo "  samtools:      $(${SAMTOOLS} --version 2>&1 | head -n1)"
echo "  featureCounts: $(${FEATURECOUNTS} -v 2>&1 | head -n1)"
echo "  Trimmomatic:   $(java -jar "${TRIM_JAR}" -version 2>&1 | head -n1 || true)"

if [[ ! -s "${GENOME_FASTA}" || ! -s "${GTF}" ]]; then
  echo "[ERROR] Genome FASTA or GTF not found. Check paths." >&2
  exit 1
fi

############################################
# 2) BUILD STAR INDEX (once)
############################################
if [[ ! -d "${STAR_INDEX}" || -z "$(ls -A "${STAR_INDEX}" 2>/dev/null)" ]]; then
  echo "[INFO] STAR index not found at ${STAR_INDEX}. Building..."
  mkdir -p "${STAR_INDEX}"
  "${STAR}" --runMode genomeGenerate \
    --runThreadN "${THREADS}" \
    --genomeDir "${STAR_INDEX}" \
    --genomeFastaFiles "${GENOME_FASTA}" \
    --sjdbGTFfile "${GTF}" \
    --sjdbOverhang "${SJDB_OVERHANG}"
  echo "[INFO] STAR index built."
else
  echo "[INFO] Using existing STAR index at ${STAR_INDEX}"
fi

############################################
# 3) FUNCTIONS
############################################
merge_lanes() {
  local sample="$1"
  local r1_out="00_merged/${sample}_R1.fastq"
  local r2_out="00_merged/${sample}_R2.fastq"

  # Collect R1/R2 lane files (gz and non-gz)
  local r1_files=()
  local r2_files=()
  for lane in "${LANES[@]}"; do
    # Try gz first, then plain fastq
    for ext in ".fastq.gz" ".fastq"; do
      local r1="${sample}_${lane}_R1_001${ext}"
      local r2="${sample}_${lane}_R2_001${ext}"
      [[ -f "$r1" ]] && r1_files+=("$r1")
      [[ -f "$r2" ]] && r2_files+=("$r2")
    done
  done

  if [[ ${#r1_files[@]} -eq 0 || ${#r2_files[@]} -eq 0 ]]; then
    echo "[ERROR] No FASTQ lanes found for ${sample}. Expected files like ${sample}_L001_R1_001.fastq(.gz)." >&2
    exit 1
  fi

  echo "[INFO] Merging lanes for ${sample} (R1: ${#r1_files[@]} files, R2: ${#r2_files[@]} files)"
  # Detect gz by extension
  if [[ "${r1_files[0]}" == *.gz ]]; then
    zcat "${r1_files[@]}" > "${r1_out}"
    zcat "${r2_files[@]}" > "${r2_out}"
  else
    cat  "${r1_files[@]}" > "${r1_out}"
    cat  "${r2_files[@]}" > "${r2_out}"
  fi
}

trim_trimmomatic() {
  local sample="$1"
  local r1_in="00_merged/${sample}_R1.fastq"
  local r2_in="00_merged/${sample}_R2.fastq"
  local r1_p="01_trimmed/${sample}_R1_paired.fastq"
  local r1_u="01_trimmed/${sample}_R1_unpaired.fastq"
  local r2_p="01_trimmed/${sample}_R2_paired.fastq"
  local r2_u="01_trimmed/${sample}_R2_unpaired.fastq"

  echo "[INFO] Trimming ${sample} with Trimmomatic..."
  java -jar "${TRIM_JAR}" PE -threads "${THREADS}" -phred33 \
    "${r1_in}" "${r2_in}" \
    "${r1_p}" "${r1_u}" \
    "${r2_p}" "${r2_u}" \
    ILLUMINACLIP:${ADAPTERS}:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
}

map_star() {
  local sample="$1"
  local r1_p="01_trimmed/${sample}_R1_paired.fastq"
  local r2_p="01_trimmed/${sample}_R2_paired.fastq"
  local prefix="02_star/${sample}_"

  echo "[INFO] Mapping ${sample} with STAR..."
  "${STAR}" \
    --runThreadN "${THREADS}" \
    --genomeDir "${STAR_INDEX}" \
    --readFilesIn "${r1_p}" "${r2_p}" \
    --outFileNamePrefix "${prefix}" \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --twopassMode Basic \
    --outFilterMultimapNmax 20 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outSAMattrRGline ID:${sample} SM:${sample} PL:ILLUMINA LB:${sample}_lib PU:lane_merge

  "${SAMTOOLS}" index "02_star/${sample}_Aligned.sortedByCoord.out.bam"
}

############################################
# 4) MAIN LOOP
############################################
if [[ ! -f samples.txt ]]; then
  echo "[ERROR] samples.txt not found. Create it with one sample ID per line (e.g., 16_S5)." >&2
  exit 1
fi

while read -r SAMPLE; do
  [[ -z "${SAMPLE}" ]] && continue
  echo "================ ${SAMPLE} ================"
  merge_lanes "${SAMPLE}"
  trim_trimmomatic "${SAMPLE}"
  map_star "${SAMPLE}"
done < samples.txt

############################################
# 5) FEATURECOUNTS
############################################
echo "[INFO] Running featureCounts on all BAMs..."
ls 02_star/*_Aligned.sortedByCoord.out.bam > 03_counts/bam.list

${FEATURECOUNTS} -T "${THREADS}" \
  -a "${GTF}" -t "${FC_GTF_FEATURE}" -g "${FC_GTF_GROUP}" \
  -p -B -C -s "${STRAND}" \
  -o 03_counts/read_counts.mm39.txt \
  $(cat 03_counts/bam.list)

echo "[INFO] All done. Outputs:"
echo "  - BAMs:      02_star/<SAMPLE>_Aligned.sortedByCoord.out.bam(.bai)"
echo "  - Counts:    03_counts/read_counts.mm39.txt"
echo "  - STAR logs: 02_star/<SAMPLE>_Log.final.out"
