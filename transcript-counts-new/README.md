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
