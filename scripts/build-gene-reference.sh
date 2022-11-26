#!/usr/bin/env bash
# This script builds a reference gene dataset from the GRCh38_ERCC92_RNA2021 reference GTF file
# using the [`build-ref`](https://gitlab.internal.sanger.ac.uk/ad33/build-ref) script.

# Define the directories:
BASE_DIR="/lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression"
METADATA_DIR="${BASE_DIR}/metadata/genes"
DERMATLAS_ANNOT_DIR="/lustre/scratch119/realdata/mdt1/team113/projects/DERMATLAS/metadata/gene-metadata/annotations"

# Make sure the output folder exists:
mkdir -p ${METADATA_DIR}

# Define the input datasets:
REF_FILE="/nfs/cancer_ref02/human/GRCh38_ERCC92_RNA2021/cgpRna/103/ensembl.gtf"
CHR_FILE="${METADATA_DIR}/chr-map.txt"
ACHILLES="${DERMATLAS_ANNOT_DIR}/common-essentials/achilles-common-essentials-103.txt"
MSK_IMPACT="${DERMATLAS_ANNOT_DIR}/onco-kb/msk-impact-103.txt"

# Define the output file:
OUT_FILE="${METADATA_DIR}/gene-length-GRCh38_103.txt"

# Make sure the input files exist:
if [ ! -f ${REF_FILE} ]
then
    echo "ERROR: input reference GTF ${REF_FILE} does not exist"
    exit 1
fi

if [ ! -f ${CHR_FILE} ]
then
    echo "ERROR: input chromosome mapping ${CHR_FILE} does not exist"
    exit 1
fi

if [ ! -f ${ACHILLES} ]
then
    echo "ERROR: annotation file ${ACHILLES} does not exist"
    exit 1
fi

if [ ! -f ${MSK_IMPACT} ]
then
    echo "ERROR: annotation file ${MSK_IMPACT} does not exist"
    exit 1
fi

# Build the reference file:
build-ref -u -vv -aMSK:${MSK_IMPACT} -aACE:${ACHILLES} -c ${CHR_FILE} ${REF_FILE} > ${OUT_FILE}
