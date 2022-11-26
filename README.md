# RNASeq analysis of project 6406 (SIN3B role in melanoma resistance and progression)

This folder contains the RNASeq analysis for Larissa's ***SIN3B*** project

## Project Identifiers

* Sequencescape ID [`6406`](http://sequencescape.psd.sanger.ac.uk/studies/6406)
* canapps ID [`2581`](https://canapps.sanger.ac.uk/action/Cancer_Pipeline_ProjectViewer?project_id=2581)

## Analysis Steps

### 1: Obtain File Locations from Canapps

* This step uses the `nst_links` system to determine which files are present for the project. This creates the `./metadata/2581-samples.txt` and `./metadata/2581-files.json` files.
* Sample IDs are extracted directly from the CANAPPS data using [`canapps-data`](https://gitlab.internal.sanger.ac.uk/ad33/canapps-data) version 2.2.6.
* See `./scripts/fetch-canapps-metadata.sh` for the script.

### 2: Run Sequencing Statistics QC

* This step uses the [`canapps-bamstats`](https://gitlab.internal.sanger.ac.uk/ad33/canapps-data) version 2.2.6 and the [`jsontable`](https://gitlab.internal.sanger.ac.uk/DERMATLAS/jsontable) version 0.1.18 scripts to collate the Canapps BAM statistics from `nst_links`.
* Once generated, the sumamrised data are plotted using R.
* The output plots are located in `./qc/sequencing/results`
* See `./scripts/fetch-sequencing-qc.sh` for the script.
* See `./qc/sequencing/README.md` for further details.

### 3: Build the Gene Metadata reference from the Alignment Genome
* In order to generate TPM counts and provide sensible gene annoattion, we need to build a reference gene file from the input reference GTF file.
* This is done using [`build-ref`](https://gitlab.internal.sanger.ac.uk/ad33/build-ref.git) version 0.1.3.
* The complete gene reference data are located in `./metadata/genes/`.
* See `./scripts/build-gene-reference.sh` for the script.
* See `./metadata/genes/README.md` for further details.

### 4: Obtain the HTSeq Counts files from Canapps

* This step uses the `2581-files.json` file to locate the Canapps HTSeq counts files for each sample and downloads them all to `./transcript-counts/counts`.
* This is done using [`sample-tpm`](https://gitlab.internal.sanger.ac.uk/DERMATLAS/sample-tpm) version 1.2.1.
* The transcript counts data are located in `./transcript-counts`.
* See `./scripts/fetch-transcript-counts.sh` for the script.
* See `./transcript-counts/README.md` for further details.

### 5: Run DESeq

* This step uses DESeq version 1.36.0 to run the downstream RNASeq analysis.
* Two sets of analyses were performed: all lines together (`all-lines`), and each line separately (`by-line`).
* The DESeq data are located in `deseq`.
* See `deseq/README.md` for more details.
