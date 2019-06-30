#!/usr/bin/env bash

set -evx

SCRIPTDIR=$(dirname $(which $0))
if [ -z "${EVALROOT}" ]; then
    EVALROOT=$(cd "${SCRIPTDIR}/../" && pwd)
fi

sudo docker run -it \
    -v ${EVALROOT}/smcounter2/google-download/data/:/srv/qgen/data/ \
    -v ${EVALROOT}/smcounter2/google-download/example/:/srv/qgen/example/ \
    -v ${EVALROOT}/smcounter2/google-download/code/:/srv/qgen/code/ \
    -v ${EVALROOT}/smcounter2/data/SRP153933:/srv/qgen/SRP153933/ \
    qiaseq/qiaseq-dna

cd /srv/qgen/SRP153933/baseline/N0261-SRR7526728/
python /srv/qgen/code/qiaseq-dna/run_qiaseq_dna.py run_sm_counter_v2.params.txt v2 single N0261 > run.log 2>&1 &  

