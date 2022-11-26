#!/usr/bin/env bash
# This script fetches the htseq-counts data from canapps for each sample in the SIN3B project,
# splits the metadata and collates the complete counts data.
# This is done using the `sample-tpm` script (https://gitlab.internal.sanger.ac.uk/DERMATLAS/sample-tpm).
# NB: This script uses `jq` version 1.5

# Define the study ID:
STUDY_ID="2581"

# Define the directories:
BASE_DIR="/lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression"
METADATA_DIR="${BASE_DIR}/metadata"

OUTPUT_DIR="${BASE_DIR}/transcript-counts"
COUNT_DIR="${OUTPUT_DIR}/counts"
METACOUNT_DIR="${OUTPUT_DIR}/metacounts"
SUMMARISED_DIR="${OUTPUT_DIR}/summarised"

# Create the transcript counts directory:
mkdir -p ${COUNT_DIR} ${METACOUNT_DIR} ${SUMMARISED_DIR}

# Deifine the inputs we need:
INPUT_SAMPLES="${METADATA_DIR}/${STUDY_ID}-samples.txt"
INPUT_FILES="${METADATA_DIR}/${STUDY_ID}-files.json"
LENGTH_FILE="${METADATA_DIR}/genes/gene-length-GRCh38_103.txt"

# Make sure we can use jq
if [[ ! -x $(which jq) ]]
then
    echo "jq executable not found"
    exit 1
fi

# Make sure we can use the sample-tpm
if [[ ! -x $(which sample-tpm) ]]
then
    echo "sample-tpm executable not found"
    exit 1
fi

# Make sure input sample file exists:
if [ ! -f ${INPUT_SAMPLES} ]
then
    echo "ERROR: input sample file ${INPUT_SAMPLES} does not exist"
    exit 1
fi

# Make sure the canapps file data exists:
if [ ! -f ${INPUT_FILES} ]
then
    echo "ERROR: input Canapps location file ${INPUT_FILES} does not exist"
    exit 1
fi

# Make sure the length annotation file exists:
if [ ! -f ${LENGTH_FILE} ]
then
    echo "ERROR: input gene annotation file ${LENGTH_FILE} does not exist"
    exit 1
fi

# Fetch each sample in turn:
while read SAMPLE
do
    echo "processing sample ${SAMPLE}..."
    HTSEQ_FILE=$(jq -r ".samples.${SAMPLE}.counts.htseq" ${INPUT_FILES})
    if test -f ${HTSEQ_FILE}
    then
        sample-tpm --gene-length-column 6 --fpkm -m ${METACOUNT_DIR}/${SAMPLE}-metacounts.txt ${LENGTH_FILE} ${HTSEQ_FILE} > ${COUNT_DIR}/${SAMPLE}-counts.txt
    else
        echo "$sample ${SAMPLE} htseq file ${HTSEQ_FILE} not found"
    fi
done < ${INPUT_SAMPLES}

# Collate the counts:
SAMPLE_REGEX="(A2058|A375|SKMEL28).+N[[:digit:]]"
collate-counts -r${SAMPLE_REGEX} -d0 -ccount ${COUNT_DIR}/*-counts.txt > ${SUMMARISED_DIR}/${STUDY_ID}-feature-count-v103.txt
collate-counts -r${SAMPLE_REGEX} -d4 -ctpm ${COUNT_DIR}/*-counts.txt > ${SUMMARISED_DIR}/${STUDY_ID}-feature-tpm-v103.txt
collate-counts -r${SAMPLE_REGEX} -d4 -cfpkm ${COUNT_DIR}/*-counts.txt > ${SUMMARISED_DIR}/${STUDY_ID}-feature-fpkm-v103.txt
collate-counts -r${SAMPLE_REGEX} -d0 -ccount ${METACOUNT_DIR}/*-metacounts.txt > ${SUMMARISED_DIR}/${STUDY_ID}-metafeature-count-v103.txt
