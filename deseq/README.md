# ***SIN3B*** DESeq2 Analysis

This directory holds the scripts and results for the DESeq2 analysis of the ***SIN3B*** data.

**NB**: These analyses require DESeq version 1.36.0

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

## Methods

* Per-gene read counts were generated from the aligned BAM files using htseq-count;
* Raw per-gene counts were filtered to remove any genes with fewer than 10 counts in total across all samples;
* MDS (using both a Euclidean distrance and a Spearman dissimilarity metric) was used to assess the global behaviour of all samples.  Two clones (C18 & C21) showed aberrant behaviour and were removed from downstream analysis;
* DESeq 1.36.0 was used to perform two analyses: one comparing all cell lines together, and one comparing each line individually;
* The all lines analysis used the model expr ~ line + condition, whilst the indiviual line analysis used expr ~ condition;
* After DESeq model fitting, apeglm[1] was used to estimate posterior effect size estimators

[1]: Zhu A, Ibrahim JG, Love MI (2018). “Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences.” Bioinformatics. doi: 10.1093/bioinformatics/bty895.
