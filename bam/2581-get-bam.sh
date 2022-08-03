#!/usr/bin/env bash

# Define the project BASE directory:
BASE_DIR="/lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression"

# Define the output directories:
BAM_DIR="${BASE_DIR}/bam/bam"
# METACOUNT_DIR="${BASE_DIR}/transcript-counts/metacounts"

# Create the output directories:
mkdir -p ${BAM_DIR}

# Process each sample in turn:
while read SAMPLE
do
    echo "processing sample ${SAMPLE}..."
    BAM_FILE=$(jq -r ".samples.${SAMPLE}.dupmarked.bam" "${BASE_DIR}/metadata/2581-files.json")
    BAI_FILE=$(jq -r ".samples.${SAMPLE}.dupmarked.bai" "${BASE_DIR}/metadata/2581-files.json")
    cp ${BAM_FILE} ${BAM_DIR}/${SAMPLE}.bam
    cp ${BAI_FILE} ${BAM_DIR}/${SAMPLE}.bam.bai
done < ${BASE_DIR}/metadata/2581-samples.txt
echop "done."
