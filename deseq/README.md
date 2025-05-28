# ***SIN3B*** DESeq2 Analysis

This directory holds the scripts and results for the DESeq2 analysis of the ***SIN3B*** data.

**NB**: This analysis require DESeq version 1.36.0

A complete analysis where all lines are considered was run, and the model `~line+condition` is fitted (`all-lines`).

* The analysis script is located in `deseq/scripts`.
* The output dataset is located in `deseq/data`.

A Principal Component Analysis (PCA) was run.

* The analysis script is located in `deseq/scripts`.
* The PCA plot in PDF format is located in `deseq/scripts/results`.

The full analysis contains:

* `gsea`: The GSEA files necessary for loading the data into GSEA;
* `qc`: QC plots for the analysis;
* `rdata`: The analysis data in `.rdata` format (NB: Read into R using `readRDS()`);
* `results`: The `apeglm`-processed DESeq2 analysis in CSV format.

## Methods

* Per-gene read counts were generated from the aligned BAM files using htseq-count;
* Raw per-gene counts were filtered to remove any genes with fewer than 10 counts in total across all samples;
* A Principal Component Analysis (PCA) was used to assess the global behaviour of all samples.  Two clones (C18 & C21) showed aberrant behaviour and were removed from downstream analysis;
* DESeq 1.36.0 was used to perform an analysis comparing all cell lines together;
* The all lines analysis used the model expr ~ line + condition;
* After DESeq model fitting, apeglm[1] was used to estimate posterior effect size estimators

[1]: Zhu A, Ibrahim JG, Love MI (2018). “Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences.” Bioinformatics. doi: 10.1093/bioinformatics/bty895.
