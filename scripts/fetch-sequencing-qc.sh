#!/usr/bin/env bash
# This script fetches the Canapps BAM QC data from canapps for each sample in the SIN3B project
# This is done using the `canapps-bamstats` script (https://gitlab.internal.sanger.ac.uk/DERMATLAS/canapps-data).
# NB: This script uses `jq` version 1.5
# NB: This script uses `jsontable` (https://gitlab.internal.sanger.ac.uk/DERMATLAS/jsontable) version 0.1.18

# Define the study ID:
STUDY_ID="2581"

# Define the directories:
BASE_DIR="/lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression"
METADATA_DIR="${BASE_DIR}/metadata"

OUTPUT_DIR="${BASE_DIR}/qc/sequencing/bamstats"
SUMMARISED_DIR="${BASE_DIR}/qc/sequencing/summarised"
RESULTS_DIR="${BASE_DIR}/qc/sequencing/results"

# Create the transcript counts directory:
mkdir -p ${OUTPUT_DIR} ${SUMARISED_DIR} ${RESULTS_DIR}

# Define the inputs we need:
INPUT_SAMPLES="${METADATA_DIR}/${STUDY_ID}-samples.txt"
INPUT_FILES="${METADATA_DIR}/${STUDY_ID}-files.json"
SUMMARY_TOML="${METADATA_DIR}/sample-canapps-stats-RNA.toml"
SUMMARY_FILE="${SUMMARISED_DIR}/${STUDY_ID}-bamstats.txt"
THRESHOLDS_FILE="${BASE_DIR}/qc/sequencing/2581-qc-thresholds.json"
PLOT_FILE="${BASE_DIR}/scripts/plot-sequencing-qc.R"

# Make sure we can use jq
if [[ ! -x $(which jq) ]]
then
    echo "jq executable not found"
    exit 1
fi

# Make sure we can use the canapps-bamstats script:
if [[ ! -x $(which canapps-bamstats) ]]
then
    echo "canapps-bamstats executable not found"
    exit 1
fi

# Make sure we can use the jsontable script:
if [[ ! -x $(which jsontable) ]]
then
    echo "jsontable executable not found"
    exit 1
fi

# Make sure we can see the jsontable script:
if [[ ! -f ${SUMMARY_TOML} ]]
then
    echo "jsontable sample summary TOML file ${SUMMARY_TOML} not found"
    exit 1
fi

# Pull in the BAM statistics:
while read SAMPLE
do
    echo "processing sample ${SAMPLE}..."
    BAS_FILE=$(jq -r ".samples.${SAMPLE}.dupmarked.bamstats" ${INPUT_FILES})
    if test -f ${BAS_FILE}
    then
        canapps-bamstats ${BAS_FILE} > ${OUTPUT_DIR}/${SAMPLE}-bamstats.json
    else
        echo "$sample ${SAMPLE} .BAS file ${BAS_FILE} not found"
    fi
done < ${INPUT_SAMPLES}

# Summarise the data using jsontable:
jsontable -vv -mi ${SUMMARY_TOML} ${OUTPUT_DIR}/*-bamstats.json > ${SUMMARY_FILE}

# Plot the data using R:
module load r-base/4.0.2.1
${PLOT_FILE} -v -p"2581-v103-" ${SUMMARY_FILE} ${THRESHOLDS_FILE} ${RESULTS_DIR}
