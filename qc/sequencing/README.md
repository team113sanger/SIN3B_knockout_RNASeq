# Input Sequencing and QC Metrics

This folder performs QC on the sample BAM file metrics.

## Data Generation

QC metrics are collected for all samples from canapps:

~~~bash
cd $BASE/qc/sequencing
canapps-bamstats $BASE/metadata/2581-files.json > data/2581-bamstats.txt
canapps-dupstats $BASE/metadata/2581-files.json > data/2581-dupstats.txt
~~~

## QC Analysis

The individual plots & threshold analyses are generated form the input datasets using the R script located in the `scripts/` directory. Results are saved to the `results/` directory.

~~~bash
cd ${BASE}/qc/sequencing
mkdir -p ./results
module load r-base/4.0.2.1
./scripts/plot-bamstats -v -p"2581-v103-" ./data/2581-bamstats.txt ./data/2581-qc-thresholds.json ./results > ./results/2581-v103-sequencing-qc-results.txt
~~~
