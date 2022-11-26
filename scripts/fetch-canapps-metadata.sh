#/usr/bin/env bash
# This script fetches the canapps samples and associated file locations from nst_links 
# using the [`canapps-data`](https://gitlab.internal.sanger.ac.uk/ad33/canapps-data) scripts.

# Define the directories:
BASE_DIR="/lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression"
METADATA_DIR="${BASE_DIR}/metadata"

# Make sure the metadata directory exists:
mkdir -p ${METADATA_DIR}

# Set the Canapps ID to grab:
CANAPPS_ID="2581"

# Pull the metadata:
canapps-data -s ${CANAPPS_ID} > ${METADATA_DIR}/${CANAPPS_ID}-samples.txt
canapps-data ${CANAPPS_ID} > ${METADATA_DIR}/${CANAPPS_ID}-files.json
