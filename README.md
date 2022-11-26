# RNASeq analysis of project 6406 (SIN3B role in melanoma resistance and progression)

This folder contains the RNASeq analysis for Larissa's ***SIN3B*** project

## Project Identifiers

* Sequencescape ID [`6406`](http://sequencescape.psd.sanger.ac.uk/studies/6406)
* canapps ID [`2581`](https://canapps.sanger.ac.uk/action/Cancer_Pipeline_ProjectViewer?project_id=2581)

## Analysis Steps

### 1: Obtain File Locations from Canapps

* This step uses the `nst_links` system to determine which files are present for the project. This creates the `./metadata/2581-samples.txt` and `./metadata/2581-files.json` files.
* Sample IDs are extracted directly from the CANAPPS data using [`canapps-data`](https://gitlab.internal.sanger.ac.uk/ad33/canapps-data) version 2.2.6.
* See [`./scripts/fetch-canapps-metadata.sh`] for the script.

### 2: Build the Gene Metadata reference from the Alignment Genome

* In order to generate TPM counts and provide sensible gene annoattion, we need to build a reference gene file from the input reference GTF file.
* This is done using [`build-ref`](https://gitlab.internal.sanger.ac.uk/ad33/build-ref.git) version 0.1.3.
* The complete gene reference data are located in `./metadata/genes/`.
* See `./scripts/build-gene-reference.sh` for the script.
* See `./metadata/genes/README.md` for further details.

### 3: Obtain the HTSeq Counts files from Canapps

* This step uses the `2581-files.json` file to locate the Canapps HTSeq counts files for each sample and downloads them all to `./transcript-counts/counts`.
* This is done using [`sample-tpm`](https://gitlab.internal.sanger.ac.uk/DERMATLAS/sample-tpm) version 1.2.1.
* The transcript counts data are located in `./transcript-counts`.
* See `./scripts/fetch-transcript-counts.sh` for the script.
* See `./transcript-counts/README.md` for further details.

<!-- # Raw Transcript Counts via CASM

## HTSeq Data Collation

The transcript counts data are pulled directly from canapps using [sample-tpm](https://gitlab.internal.sanger.ac.uk/ad33/sample-tpm). This splits the counts and metacounts, writing the counts and metacounts to the `./counts` and `./metacounts` directories respectively.

This process is performed by the `get-counts.sh` script:

~~~bash
export BASE="/lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression"
cd ${BASE}/transcript-counts
./get-counts.sh
~~~

## Counts Matrix Collation

Once individual counts data have been generated, the counts matrices are generated as:

~~~bash
export BASE="/lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression"
cd ${BASE}/transcript-counts-casm
mkdir ${BASE}/transcript-counts-casm/summarised
SAMPLE_REGEX="[AS].+N[[:digit:]]"
collate-counts -r${SAMPLE_REGEX} -ccount ./counts/*-counts.txt > ./summarised/2581-feature-count-v103.txt
collate-counts -r${SAMPLE_REGEX} -ctpm ./counts/*-counts.txt > ./summarised/2581-feature-tpm-v103.txt
collate-counts -r${SAMPLE_REGEX} -cfpkm ./counts/*-counts.txt > ./summarised/2581-feature-fpkm-v103.txt
collate-counts -r${SAMPLE_REGEX} -ccount ./metacounts/*-metacounts.txt > ./summarised/2581-metafeature-count-v103.txt
~~~ -->
