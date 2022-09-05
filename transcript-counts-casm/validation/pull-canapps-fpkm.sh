#!/usr/bin/env bash

# Define the project BASE directory:
BASE_DIR="/lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression"

# Define the output directories:
OUTPUT_DIR="${BASE_DIR}/transcript-counts/validation/canapps-fpkm"

# Process each sample in turn:
while read SAMPLE
do
    echo "processing sample ${SAMPLE}..."
    HTSEQ_FILE=$(jq -r ".samples.${SAMPLE}.counts.counts" "${BASE_DIR}/metadata/2581-files.json")
    if test -f ${HTSEQ_FILE}
    then
        cp ${HTSEQ_FILE} ${OUTPUT_DIR}/${SAMPLE}-converted-counts.txt.gz
        gunzip ${OUTPUT_DIR}/${SAMPLE}-converted-counts.txt.gz
    else
        echo "WARNING: no file for sample ${SAMPLE}"
    fi
done < ${BASE_DIR}/metadata/2581-samples.txt
