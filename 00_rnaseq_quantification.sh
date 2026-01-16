#!/bin/bash
# 00_rnaseq_quantification.sh
# RNA-Seq Analysis Protocol: Raw Reads -> FastQC -> fastp -> SortMeRNA -> BBMap(Repair) -> STAR -> featureCounts

set -e
set -o pipefail

# --- Configuration ---
THREADS=12

# Input/Output Directories
RAW_DIR="../../raw_data"
QC_DIR="../../qc"
FASTP_DIR="../../fastp"
CLEAN_DIR="../../fastq_clean"
MAPPED_DIR="../../mapped"

# Reference Paths (Adjust these!)
SORTMERNA_DB="$HOME/sortmerna/database/smr_v4.3_default_db.fasta"
SORTMERNA_WORK_DIR="$HOME/sortmerna/run"
GENOME_DIR="../../ref"            # STAR index directory
GTF_FILE="../../ref/annotation.gtf" # Genomic annotation

# Output
FEATURECOUNTS_OUT="count.out"

# Create directories
mkdir -p "$QC_DIR/pre_fastqc" "$QC_DIR/post_fastp_fastqc"
mkdir -p "$FASTP_DIR/trimmed_fastq" "$FASTP_DIR/reports"
mkdir -p "$CLEAN_DIR" "$MAPPED_DIR" "$SORTMERNA_WORK_DIR"

# --- 1. Quality Control (Pre-Trim) ---
echo "Step 1: Pre-Trim QC..."
if compgen -G "$RAW_DIR/*.fastq.gz" > /dev/null; then
    fastqc "$RAW_DIR"/*.fastq.gz -o "$QC_DIR/pre_fastqc" -t $THREADS
    multiqc "$QC_DIR/pre_fastqc" -o "$QC_DIR/pre_fastqc"
else
    echo "Warning: No files found in $RAW_DIR"
fi

# --- 2. Trimming (fastp) ---
echo "Step 2: Trimming with fastp..."
for r1 in "$RAW_DIR"/*_1.fastq.gz; do
    [ -e "$r1" ] || continue
    base=$(basename "$r1" _1.fastq.gz)
    r2="$RAW_DIR/${base}_2.fastq.gz"
    
    echo "  Processing $base..."
    fastp -i "$r1" -I "$r2" \
          -o "$FASTP_DIR/trimmed_fastq/${base}_1.trimmed.fastq.gz" \
          -O "$FASTP_DIR/trimmed_fastq/${base}_2.trimmed.fastq.gz" \
          -h "$FASTP_DIR/reports/${base}_fastp.html" \
          -j "$FASTP_DIR/reports/${base}_fastp.json" \
          --thread 4
done

# Post-Trim QC
fastqc "$FASTP_DIR/trimmed_fastq"/*.gz -o "$QC_DIR/post_fastp_fastqc" -t $THREADS
multiqc "$QC_DIR/post_fastp_fastqc" -o "$QC_DIR/post_fastp_fastqc"

# --- 3. rRNA Removal (SortMeRNA) ---
echo "Step 3: rRNA Removal with SortMeRNA..."
if [ ! -f "$SORTMERNA_DB" ]; then echo "Error: SortMeRNA DB not found!"; exit 1; fi

for r1 in "$FASTP_DIR/trimmed_fastq"/*_1.trimmed.fastq.gz; do
    [ -e "$r1" ] || continue
    base=$(basename "$r1" _1.trimmed.fastq.gz)
    r2="$FASTP_DIR/trimmed_fastq/${base}_2.trimmed.fastq.gz"
    
    echo "  Filtering $base..."
    
    rm -rf "$SORTMERNA_WORK_DIR"/*
    
    sortmerna \
        --ref "$SORTMERNA_DB" \
        --reads "$r1" \
        --reads "$r2" \
        --paired_in \
        --fastx \
        --other "$SORTMERNA_WORK_DIR/${base}_clean" \
        --aligned "$SORTMERNA_WORK_DIR/${base}_rrna" \
        --threads $THREADS \
        --workdir "$SORTMERNA_WORK_DIR"
    
    mv "$SORTMERNA_WORK_DIR"/${base}_clean*.fq "$CLEAN_DIR/" 2>/dev/null || true
done

# --- 4. Repair / De-interleave (BBMap) ---
echo "Step 4: De-interleaving with repair.sh..."
for f in "$CLEAN_DIR"/*clean*.fq; do
    [ -e "$f" ] || continue
    if [[ "$f" == *"_R1.fq"* ]] || [[ "$f" == *"_R2.fq"* ]]; then continue; fi

    base=$(basename "$f" .fq) 
    
    if [[ "$base" == *"_1"* ]] || [[ "$base" == *"_2"* ]]; then
        continue 
    fi

    echo "  Repairing $base..."
    repair.sh in="$f" \
              out1="$CLEAN_DIR/${base}_R1.fq.gz" \
              out2="$CLEAN_DIR/${base}_R2.fq.gz" \
              overwrite=t
done

# --- 5. Mapping (STAR) ---
echo "Step 5: Mapping with STAR..."
for r1 in "$CLEAN_DIR"/*_R1.fq.gz; do
    [ -e "$r1" ] || continue
    base=$(basename "$r1" _R1.fq.gz)
    r2="$CLEAN_DIR/${base}_R2.fq.gz"
    
    echo "  Mapping $base..."
    STAR --runMode alignReads \
         --genomeDir "$GENOME_DIR" \
         --readFilesIn "$r1" "$r2" \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN $THREADS \
         --outFileNamePrefix "$MAPPED_DIR/${base}_"
done

# --- 6. Quantification (featureCounts) ---
echo "Step 6: Counting features..."
if [ -n "$(ls -A $MAPPED_DIR/*.bam 2>/dev/null)" ]; then
    featureCounts -a "$GTF_FILE" \
                  -o "$FEATURECOUNTS_OUT" \
                  -T $THREADS \
                  -p \
                  "$MAPPED_DIR"/*.bam
else
    echo "  No BAM files found to count."
fi

echo "Pipeline complete. Output: $FEATURECOUNTS_OUT"
