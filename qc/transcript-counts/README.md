# HTSeq Transcript Count QC

The individual plots & threshold analyses are generated form the input datasets using the R script located in the `scripts/` directory. Results are saved to the `results/` directory.

The process is simply:

1. to calculate the per-sample totals from the real counts data;
2. to normalise the metacounts to the totals; and then
3. to plot the metafeatures

~~~bash
cd ${BASE}/qc/transcript-counts
mkdir -p ./results
module load r-base/4.0.2.1
./scripts/plot-transcriptstats -v -p"2581-v103-" $BASE/transcript-counts/summarised/2581-feature-count-v103.txt $BASE/transcript-counts/summarised/2581-metafeature-count-v103.txt ./results
~~~
