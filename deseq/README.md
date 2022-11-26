# ***SIN3B*** DESeq2 Analysis

This directory holds the scripts and results for the DESeq2 analysis of the ***SIN3B*** data.

Two analysis sets are run:

1. A complete analysis where all lines are considered, and the model `~line+condition` is fitted (`all-lines`); and
2. Three analyses, one for each line where each line is treated separately, and the model `~condition` is fitted (`by-line`)

* The two analysis scripts are located in `deseq/scripts`.
* The output datasets ate located in `deseq/data`.

Each analysis contains:

* `gsea`: The GSEA files necessary for loading the data into GSEA;
* `qc`: QC plots for the analysis;
* `rdata`: The analysis data in `.rdata` format (NB: Read into R using `readRDS()`);
* `results`: The `apeglm`-processed DESeq2 analysis in CSV format.
