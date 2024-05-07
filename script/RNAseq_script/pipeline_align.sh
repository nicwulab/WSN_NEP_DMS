#!/bin/bash

# Set variables
PROJECT_PATH='/home/yiquan2/NEP'
CODE_PATH="$PROJECT_PATH/script"
REF="$PROJECT_PATH/data/ref"
ANNOT="$PROJECT_PATH/data/ref/flu_WSN_annotation.gff3"
RESULT_PATH="$PROJECT_PATH/result"

# Create an array of sample names
SAMPLES=( $(ls "$PROJECT_PATH/data/raw" | grep "_L00M_R1_001.fastq.gz" | sed 's/_L00M_R1_001.fastq.gz//') )

# Iterate over sample names
for SAMPLENAME in "${SAMPLES[@]}"; do
    OUTPUT_BAM="$RESULT_PATH/${SAMPLENAME}/${SAMPLENAME}Aligned.sortedByCoord.out.bam"

    # Check if the output file already exists, and skip if it does
    if [ -e "$OUTPUT_BAM" ]; then
        echo "Output file for $SAMPLENAME already exists. Skipping..."
        continue
    fi

    # Create result directory
    mkdir -p "$RESULT_PATH/${SAMPLENAME}"

    # Unzip fastq files
    zcat "$PROJECT_PATH/data/raw/${SAMPLENAME}_L00M_R1_001.fastq.gz" > "$PROJECT_PATH/data/raw/${SAMPLENAME}_L00M_R1_001.fastq"
    zcat "$PROJECT_PATH/data/raw/${SAMPLENAME}_L00M_R2_001.fastq.gz" > "$PROJECT_PATH/data/raw/${SAMPLENAME}_L00M_R2_001.fastq"

    # Align with STAR
    STAR --genomeDir "$REF" \
         --readFilesIn "$PROJECT_PATH/data/raw/${SAMPLENAME}_L00M_R1_001.fastq" "$PROJECT_PATH/data/raw/${SAMPLENAME}_L00M_R2_001.fastq" \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix "$RESULT_PATH/${SAMPLENAME}/${SAMPLENAME}" \
         --sjdbGTFfile "$ANNOT" \
         --runThreadN 50
done
