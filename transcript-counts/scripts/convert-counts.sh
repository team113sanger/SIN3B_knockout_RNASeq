#/usr/bin/env bash

# Define the project BASE directory:
BASE_DIR="/lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression"

# Define the input gene count file:
LENGTH_FILE="/nfs/users/nfs_a/ad33/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression/genes/gene-length-GRCh38_103.txt"

# Define the output directories:
HTSEQ_DIR="${BASE_DIR}/transcript-counts/counts"
TPM_DIR="${BASE_DIR}/transcript-counts/tpm"
METACOUNT_DIR="${BASE_DIR}/transcript-counts/metacounts"

# Create the output directories:
mkdir -p ${TPM_DIR} ${METACOUNT_DIR}

# Process each counts file individually:
while read SAMPLE
do
    HTSEQ_FILE="${HTSEQ_DIR}/${SAMPLE}-counts.txt.gz"
    echo "processing sample ${SAMPLE} from ${HTSEQ_FILE}..."
    if test -f ${HTSEQ_FILE}
    then
        sample-tpm --gene-length-column 6 --fpkm -m ${METACOUNT_DIR}/${SAMPLE}-metacounts.txt ${LENGTH_FILE} ${HTSEQ_FILE} > ${TPM_DIR}/${SAMPLE}-counts.txt
    fi
done < ${BASE_DIR}/metadata/2581-samples.txt
