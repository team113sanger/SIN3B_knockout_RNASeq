## Transcript Metadata

The transcript metadata file `./gene-length-GRCh38_103` is created from the `ensembl.gene_length.tsv` reference file using the [`convert-genelength`](https://gitlab.internal.sanger.ac.uk/ad33/sample-tpm) script:

~~~bash
convert-genelength -vv --update-ercc --chr-map ./chr-map.txt /nfs/cancer_ref02/human/GRCh38_ERCC92_RNA2021/cgpRna/e103/ensembl.gene_length.tsv > ./gene-length-GRCh38_103.txt
~~~

The `./chr-map.txt` file provides chromomse name mapping information.
