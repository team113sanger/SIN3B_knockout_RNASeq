# Transcript COunting Outside Canapps

THis directory contains the code to re-count the original Canapps BAM files using [`htseq-count`](https://htseq.readthedocs.io/en/master/htseqcount.html) to avoid the strand bug.

## Installation of HTSeq

HTSeq was installed using `conda`:

~~~bash
conda create -y -c bioconda -n htseq python numpy pysam htseq tabix
~~~

## Running HTSeq

~~~bash
qsubsec -ps 2581-htseq-count.qsubsec 2581-htseq-count.tff SAMPLE='FILE(../../metadata/2581-samples.txt)'
~~~

## Build New Gene Reference Dataset

Before counting genes, we can build the gene reference length dataset using [`build-ref`](https://gitlab.internal.sanger.ac.uk/ad33/build-ref).

See the `README.md` file in the `./reference` directory for details.

## Generate the Transcript Counts

Individual sample TPM counts are generated using the `./scripts/convert-counts.sh` script.

## Collate Counts

Individual sample TPM counts are collated & summarised using the `./scripts/collate-counts.sh` script.
