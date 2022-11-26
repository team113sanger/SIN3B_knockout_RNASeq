# Gene Reference Data for the ***SIN3B*** Transcript Analysis Project

The ***SIN3B*** project uses htseq-count files form Canapps, but we call the TPM values ourselves to overcome the shortfallings in the Canapps CHTSEQcount calling methods. This process also generates a complete annotation dataset for the counts data that we can use later in downstream analyses.

## Input Data

The reference file is built using the [`build-ref`](https://gitlab.internal.sanger.ac.uk/ad33/build-ref.git) system.  It requires:

1. An input GTF file;
2. A mapping file to update any incorrect chromosome names; and
3. A set of annotation files to include.

### Reference GTF File

The input reference GTF file is the one associated with the Canapps reference used for the project. We use the `GRCh38_ERCC92_RNA2021` reference for this project. The input GTF is therefore located at `/nfs/cancer_ref02/human/GRCh38_ERCC92_RNA2021/cgpRna/103/ensembl.gtf`.

### Chromosome Mapping File

There are several chromosome naming diffeerences between the various analysis files used at Sanger.  To overcome this, we use a hand-crafter mapping file to update the reference names to match the HTSeqCounts data. The mapping file used is located at `./metadata/genes/chr-map.txt`. **NB**: This file was originally created for the DERMATLAS project.

### Gene Annotation Files

Two annotation datasets are included in the reference dataset.  These are both originally created for the DERMATLAS project.  See the DERMATLAS project documentation (currently at `/lustre/scratch119/realdata/mdt1/team113/projects/DERMATLAS/metadata/gene-metadata/annotations/README.md`) for the complete documentation.

## Dataset Creation

* The reference dataset is generated from the above inputs using the `./scripts/build-gene-reference.sh` script.
* The output gene reference file is `./metadata/genes/gene-length-GRCh38_103.txt`.

### Generation Log Data

~~~log
INFO - reading MSK annotation set from /lustre/scratch119/realdata/mdt1/team113/projects/DERMATLAS/metadata/gene-metadata/annotations/onco-kb/msk-impact-103.txt
INFO - loaded 504 gene IDs from /lustre/scratch119/realdata/mdt1/team113/projects/DERMATLAS/metadata/gene-metadata/annotations/onco-kb/msk-impact-103.txt
INFO - reading ACE annotation set from /lustre/scratch119/realdata/mdt1/team113/projects/DERMATLAS/metadata/gene-metadata/annotations/common-essentials/achilles-common-essentials-103.txt
INFO - loaded 2388 gene IDs from /lustre/scratch119/realdata/mdt1/team113/projects/DERMATLAS/metadata/gene-metadata/annotations/common-essentials/achilles-common-essentials-103.txt
INFO - 2 annotation datasets loaded
INFO - reading chromosome mapping from /lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression/metadata/genes/chr-map.txt
INFO - reading GTFdata from /nfs/cancer_ref02/human/GRCh38_ERCC92_RNA2021/cgpRna/103/ensembl.gtf
~~~
