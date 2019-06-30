SCRIPTDIR=$(dirname $(which $0))
if [ -z "${EVALROOT}" ]; then
    EVALROOT=$(cd "${SCRIPTDIR}/../" && pwd)
fi

REF="${EVALROOT}/datafiles/hs37d5.fa"

for srr in $(cat listSRR.txt); do echo $srr ; "${EVALROOT}/tools/uvc/debarcode" -i ${srr}*.fastq.gz -o "${srr}.barcoded.fq.gz" -b 0 -e 12 ; "${EVALROOT}/tools/bwa-0.7.17/bwa" mem -t 12 "$REF" ${srr}*.barcoded.fq.gz | samtools view -h1 | samtools sort -@2 - > "${srr}.bam" ; samtools index "${srr}.bam" ; "${EVALROOT}/tools/uvc/uvc1" -f "${REF}" "${srr}.bam" -o "${srr}.vcf.gz" ; done

