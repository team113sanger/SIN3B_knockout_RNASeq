#!/usr/bin/env bash

# Define the project BASE directory:
BASE_DIR="/lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression"

# Define the input gene count file:
LENGTH_FILE="/nfs/users/nfs_a/ad33/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression/metadata/genes/gene-length-GRCh38_103.txt"

# Define the output directories:
COUNT_DIR="${BASE_DIR}/transcript-counts-casm/counts"
METACOUNT_DIR="${BASE_DIR}/transcript-counts-casm/metacounts"

# Create the output directories:
mkdir -p ${COUNT_DIR} ${METACOUNT_DIR}

# Process each sample in turn:
while read SAMPLE
do
    echo "processing sample ${SAMPLE}..."
    HTSEQ_FILE=$(jq -r ".samples.${SAMPLE}.counts.htseq" "${BASE_DIR}/metadata/2581-files.json")
    if test -f ${HTSEQ_FILE}
    then
        sample-tpm --gene-length-column 6 --fpkm -m ${METACOUNT_DIR}/${SAMPLE}-metacounts.txt ${LENGTH_FILE} ${HTSEQ_FILE} > ${COUNT_DIR}/${SAMPLE}-counts.txt
    fi
done < ${BASE_DIR}/metadata/2581-samples.txt
