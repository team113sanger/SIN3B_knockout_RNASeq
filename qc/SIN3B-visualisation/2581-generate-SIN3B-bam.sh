#!/usr/bin/env bash

# Define the project BASE directory:
BASE_DIR="/lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression"

# Define the output directories:
BAM_DIR="${BASE_DIR}/bam/bam"
RES_DIR="${BASE_DIR}/qc/SIN3B-visualisation/bam"

# Create the output directories:
mkdir -p ${RES_DIR}

# Define the SIN3B region, rounded and padded by 1k:
GENE_REGION="chr19:16830000-16882000"
GENE_LABEL="SIN3B"

# Process each sample in turn:
while read SAMPLE
do
    echo "processing sample ${SAMPLE}..."
    OUT_BAM="${RES_DIR}/${SAMPLE}-${GENE_LABEL}.bam"
    samtools view -o ${OUT_BAM} ${BAM_DIR}/${SAMPLE}.bam ${GENE_REGION}
    samtools index ${OUT_BAM}
done < ${BASE_DIR}/metadata/2581-samples.txt
echo "done."
