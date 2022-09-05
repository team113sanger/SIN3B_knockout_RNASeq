#/usr/bin/env bash

# Define the project BASE directory:
BASE_DIR="/lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression"

# Define the sample name regular expression to use:
# SAMPLE_REGEX="[AS][^/]+N[[:digit:]]"
SAMPLE_REGEX="(A2058|A375|SKMEL28)_(SIN3B_KO|CONTROL)_C[[:digit:]]+_N[[:digit:]]"

# Define the output directories:
TPM_DIR="${BASE_DIR}/transcript-counts/tpm"
METACOUNT_DIR="${BASE_DIR}/transcript-counts/metacounts"
OUT_DIR="${BASE_DIR}/transcript-counts/summarised"

# Create the output directories:
mkdir -p ${OUT_DIR}

collate-counts -r${SAMPLE_REGEX} -d0 -ccount ${TPM_DIR}/*-counts.txt > ${OUT_DIR}/2581-feature-count-v103.txt
collate-counts -r${SAMPLE_REGEX} -d4 -ctpm ${TPM_DIR}/*-counts.txt > ${OUT_DIR}/2581-feature-tpm-v103.txt
collate-counts -r${SAMPLE_REGEX} -d4 -cfpkm ${TPM_DIR}/*-counts.txt > ${OUT_DIR}/2581-feature-fpkm-v103.txt
collate-counts -r${SAMPLE_REGEX} -d0 -ccount ${METACOUNT_DIR}/*-metacounts.txt > ${OUT_DIR}/2581-metafeature-count-v103.txt
