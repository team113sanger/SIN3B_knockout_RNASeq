# ***SIN3B*** Project Input Sequencing QC Metrics

This folder performs QC on the sample BAM file metrics.

## Data Generation

* QC metrics are collected for all samples from canapps using the [`canapps-bamstats`](https://gitlab.internal.sanger.ac.uk/DERMATLAS/canapps-data) script.
* This scriptuses the Canapps metadata in `./metadata/2581-files.json` to pull out all of the Canapps BAM statistics.
* The BAM statistics are located in `qc/sequencing/data`
* The data generation script is `scripts/fetch-sequencing-qc.sh`

## Data Summarisation

* Once the input bamstats files have been downloaded from Canapps and processed using `canapps-bamstats`, they are summarised into a single file using the [`jsontable`](https://gitlab.internal.sanger.ac.uk/DERMATLAS/jsontable) script.
* The `jsontable` script requires a TOML file that defines the input and output files required. This file is located at `metadata/sample-canapps-stats-RNA`.
* The summarised data are located at `qc/sequencing/summarised/2581-bamstats.txt`.
* See the [`jsontable documentation`](https://gitlab.internal.sanger.ac.uk/DERMATLAS/jsontable/-/blob/develop/README.md) for details of the TOML specification file format.

## QC Analysis

* Once the BAM statistics have been summarised, the QC plots are generated using the R script `scripts/plot-sequenceing-qc.R`.
* This script uses the thresholds in `qc/sequencing/2581-qc-thresholds.json` to plot the QC data
* The QC plots are located in `qc/sequencing/results`
