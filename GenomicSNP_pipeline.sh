#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Define directories
BASE_DIR="$HOME/Genomic_Feature_Analysis"
DATA_DIR="$BASE_DIR/data"
READS_DIR="$DATA_DIR/reads"
GENOME_DIR="$BASE_DIR/genome"
ANNOTATION_DIR="$GENOME_DIR/annotation"
OUTPUT_DIR="$BASE_DIR/output"
ALIGNMENT_DIR="$OUTPUT_DIR/alignment"
VARIANT_DIR="$OUTPUT_DIR/variants"
QC_DIR="$OUTPUT_DIR/qc"

# Create necessary directories
mkdir -p "$DATA_DIR" "$READS_DIR" "$GENOME_DIR" "$ANNOTATION_DIR" "$OUTPUT_DIR" "$ALIGNMENT_DIR" "$VARIANT_DIR" "$QC_DIR"

# Define sample IDs
SAMPLES=("SRR8906861" "SRR8906666")

# Step 1: Download sequencing data from NCBI
echo "Downloading raw sequencing data..."
for SAMPLE in "${SAMPLES[@]}"; do
    fasterq-dump --split-files --outdir "$READS_DIR" "$SAMPLE"
done

# Step 2: Quality Control using FastQC
echo "Performing quality control with FastQC..."
fastqc -o "$QC_DIR" "$READS_DIR"/*.fastq

# Step 3: Trimming reads using Trimmomatic
echo "Trimming reads with Trimmomatic..."
for FILE in "$READS_DIR"/*_1.fastq; do
    SAMPLE=$(basename "$FILE" _1.fastq)
    trimmomatic PE         -threads 4         "$READS_DIR/${SAMPLE}_1.fastq" "$READS_DIR/${SAMPLE}_2.fastq"         "$READS_DIR/${SAMPLE}_1.trimmed.fastq" "$READS_DIR/${SAMPLE}_1.unpaired.fastq"         "$READS_DIR/${SAMPLE}_2.trimmed.fastq" "$READS_DIR/${SAMPLE}_2.unpaired.fastq"         ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

# Step 4: Index maize reference genome using BWA
echo "Indexing the maize reference genome..."
bwa index "$GENOME_DIR/maize_reference.fa"

# Step 5: Align reads to the reference genome using BWA
echo "Aligning reads to the reference genome..."
for FILE in "$READS_DIR"/*_1.trimmed.fastq; do
    SAMPLE=$(basename "$FILE" _1.trimmed.fastq)
    bwa mem -t 4 "$GENOME_DIR/maize_reference.fa"         "$READS_DIR/${SAMPLE}_1.trimmed.fastq" "$READS_DIR/${SAMPLE}_2.trimmed.fastq"         > "$ALIGNMENT_DIR/${SAMPLE}.sam"
done

# Step 6: Convert SAM to BAM, sort, and index using SAMtools
echo "Converting SAM to BAM and indexing..."
for FILE in "$ALIGNMENT_DIR"/*.sam; do
    SAMPLE=$(basename "$FILE" .sam)
    samtools view -bS "$FILE" | samtools sort -o "$ALIGNMENT_DIR/${SAMPLE}.sorted.bam"
    samtools index "$ALIGNMENT_DIR/${SAMPLE}.sorted.bam"
    rm "$FILE"  # Remove SAM files to save space
done

# Step 7: Variant Calling using BCFtools
echo "Calling SNPs and small indels using BCFtools..."
bcftools mpileup -Ou -f "$GENOME_DIR/maize_reference.fa" "$ALIGNMENT_DIR"/*.sorted.bam | bcftools call -mv -Ob -o "$VARIANT_DIR/called_variants.bcf"

# Convert BCF to VCF format
bcftools view "$VARIANT_DIR/called_variants.bcf" > "$VARIANT_DIR/called_variants.vcf"

# Step 8: Filter SNPs based on quality
echo "Filtering SNPs based on quality..."
bcftools filter -i 'QUAL>20' "$VARIANT_DIR/called_variants.vcf" > "$VARIANT_DIR/filtered_variants.vcf"

# Step 9: Annotate SNPs with gene regions using BCFTools
echo "Annotating SNPs with gene regions..."
bcftools annotate --annotations "$ANNOTATION_DIR/maize_annotations.gff" --columns CHROM,POS,GENE "$VARIANT_DIR/filtered_variants.vcf" > "$VARIANT_DIR/annotated_variants.vcf"

# Step 10: Identify overlaps with regulatory elements
echo "Identifying overlaps with regulatory elements..."
bedtools intersect -a "$VARIANT_DIR/annotated_variants.vcf" -b "$ANNOTATION_DIR/regulatory_regions.bed" > "$VARIANT_DIR/overlapping_variants.vcf"

# Step 11: Generate summary statistics for SNP distribution
echo "Generating summary statistics for SNP distribution..."
bcftools stats "$VARIANT_DIR/filtered_variants.vcf" > "$VARIANT_DIR/variant_stats.txt"

# Step 12: Compare SNP distributions between resistant and susceptible genotypes
echo "Comparing SNP distributions between genotypes..."
vcftools --vcf "$VARIANT_DIR/filtered_variants.vcf" --weir-fst-pop "$DATA_DIR/resistant_samples.txt" --weir-fst-pop "$DATA_DIR/susceptible_samples.txt" --out "$VARIANT_DIR/fst_comparison"

echo "Genomic Feature Analysis workflow completed successfully."
