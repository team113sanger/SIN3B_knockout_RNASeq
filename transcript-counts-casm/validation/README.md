# FPKM & TPM Data Validation

This folder pulls the canapps converted counts data and compares it to the data generated via [sample-tpm](https://gitlab.internal.sanger.ac.uk/ad33/sample-tpm). Hopefully, the only difference is rounding.

## Input Data

The input canapps FPKM data are copied from `nst_links` using `pull-canapps-fpkm.sh`:

~~~bash
export BASE="/lustre/scratch119/realdata/mdt1/team113/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression"
cd ${BASE}/transcript-counts/validation
mkdir canapps-fpkm
./pull-canapps-fpkm.sh
~~~

## Validation

Once the canapps "reference" FPKM data are downloaded, the actual validation is performed with R, by running the `validate-fpkm.R` script interactively.

The resulting abssolute maximum difference vales across all samples are recorded in the `validation-results.txt` file.

## Results

~~~plain
count   overall maximum delta: 0.000000
tpm     overall maximum delta: 0.005000
fpkm    overall maximum delta: 0.005000
~~~

As expected, all counts are the same, and all TPM & FPKM values are within a maximum of 0.005. This is because canapps rounds the data and DERMATLAS does not.
