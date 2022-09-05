# Gene Length Refernce Generation

~~~bash
REF_FILE="/nfs/cancer_ref02/human/GRCh38_ERCC92_RNA2021/cgpRna/103/ensembl.gtf"
CHR_FILE="./chr-map.txt"
OUT_FILE="./gene-length-GRCh38_103.txt"
ACHILLES="/lustre/scratch119/realdata/mdt1/team113/projects/DERMATLAS/metadata/gene-metadata/annotations/common-essentials/achilles-common-essentials-103.txt"
MSK_IMPACT="/lustre/scratch119/realdata/mdt1/team113/projects/DERMATLAS/metadata/gene-metadata/annotations/onco-kb/msk-impact-103.txt"
build-ref -u -vv -aMSK:${MSK_IMPACT} -aACE:${ACHILLES} -c ${CHR_FILE} ${REF_FILE} > ${OUT_FILE}
~~~
