# Project Metadata

## Sample IDs & File Locations

Sample IDs are extracted directly from the CANAPPS data using [`canapps-data`](https://gitlab.internal.sanger.ac.uk/ad33/canapps-data):

~~~bash
CANAPPS_ID="2581"
canapps-data -s ${CANAPPS_ID} > ${CANAPPS_ID}-samples.txt
canapps-data ${CANAPPS_ID} > ${CANAPPS_ID}-files.json
~~~
