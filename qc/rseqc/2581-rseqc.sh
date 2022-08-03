#!/usr/bin/env bash

# Define the project BASE directory:
BASE_DIR="/lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression"

# Define the output directories:
BAM_DIR="${BASE_DIR}/bam/bam"
RES_DIR="${BASE_DIR}/qc/rseqc/results"

# Define the reference:
REF="${BASE_DIR}/metadata/reference/GRCh38_ERCC92_RNA2021-refseq.bed"

# Create the output directories:
mkdir -p ${RES_DIR}

# Process each sample in turn:
while read SAMPLE
do
    echo "processing sample ${SAMPLE}..."
    BAM_FILE="${BAM_DIR}/${SAMPLE}.bam"
    RES_FILE="${RES_DIR}/${SAMPLE}-rseqc.txt"
    infer_experiment.py -r ${REF} -i ${BAM_FILE} > ${RES_FILE}
done < ${BASE_DIR}/metadata/2581-samples.txt
echo "done."
