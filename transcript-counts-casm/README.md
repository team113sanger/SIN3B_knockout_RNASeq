# Raw Transcript Counts

## HTSeq Data Collation

The transcript counts data are pulled directly from canapps using [sample-tpm](https://gitlab.internal.sanger.ac.uk/ad33/sample-tpm). This splits the counts and metacounts, writing the counts and metacounts to the `./counts` and `./metacounts` directories respectively.

This process is performed by the `get-counts.sh` script:

~~~bash
export BASE="/lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression"
cd ${BASE}/transcript-counts
./get-counts.sh
~~~

## Counts Matrix Collation

Once individual counts data have been generated, the counts matrices are generated as:

~~~bash
export BASE="/lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression"
cd ${BASE}/transcript-counts-casm
mkdir ${BASE}/transcript-counts-casm/summarised
SAMPLE_REGEX="[AS].+N[[:digit:]]"
collate-counts -r${SAMPLE_REGEX} -ccount ./counts/*-counts.txt > ./summarised/2581-feature-count-v103.txt
collate-counts -r${SAMPLE_REGEX} -ctpm ./counts/*-counts.txt > ./summarised/2581-feature-tpm-v103.txt
collate-counts -r${SAMPLE_REGEX} -cfpkm ./counts/*-counts.txt > ./summarised/2581-feature-fpkm-v103.txt
collate-counts -r${SAMPLE_REGEX} -ccount ./metacounts/*-metacounts.txt > ./summarised/2581-metafeature-count-v103.txt
~~~
