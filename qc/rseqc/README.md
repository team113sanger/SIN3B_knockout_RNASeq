# RSeQC BAM File Strand Analysis

This directory contains the scripts necessary to check the input BAM files for strandedness.  This is necessary to determine if we need to mitigate against the Canapps RNASeq strandedness bug.

## Obtain the Local Reference

As we're running multiple tests, we copy the reference BED file locally as:

~~~bash
BASE_DIR="/lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression"
REF_DIR="${BASE_DIR}/metadata/reference"
mkdir -p ${REF_DIR}
cp /nfs/cancer_ref02/Homo_sapiens/GRCh38_ERCC92_RNA2021/rseqc/RefSeq.bed ${REF_DIR}/GRCh38_ERCC92_RNA2021-refseq.bed
~~~

## Running RSeQC

Make sure the RSeQC module is loaded, then run the `2581-rseqc.sh` script:

~~~bash
cd /lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression/qc/rseqc
conda activate rseqc
./2581-rseqc.sh
~~~

## Check Results

**NB**: See [here](http://rseqc.sourceforge.net/#infer-experiment-py) for an example if the interpretation of the results.

~~~bash
cd /lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression/qc/rseqc
grep 'Fraction of reads explained by "1++,1--,2+-,2-+"' ./results/*.txt
grep 'Fraction of reads explained by "1+-,1-+,2++,2--"' ./results/*.txt
~~~

All the results show a significant disparity between the two read orientations, so these samples are all **stranded**.
