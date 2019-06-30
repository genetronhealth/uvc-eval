
SCRIPTDIR=$(dirname $(which $0))
if [ -z "${EVALROOT}" ]; then
    EVALROOT=$(cd "${SCRIPTDIR}/../../" && pwd)
fi

v=$(${EVALROOT}/tools/uvc/uvc1 -v | sed 's/\./-/g')
cat ${EVALROOT}//ERP015684/misc/by518_supplementary_data_s8__*.tsv \
| python ${EVALROOT}/ERP015684/tsv2quals-ERP015684.py \
${EVALROOT}/ERP015684/vcf/T*/T*_N*.${v}.vcf.gz \
> ${EVALROOT}/ERP015684/misc/tmp/${v}.by518_supplementary_data_s8.tsv

cat ${EVALROOT}/ERP015684/misc/tmp/${v}.by518_supplementary_data_s8.tsv | python ${EVALROOT}/ERP015684/rank-scatterplot.py ${EVALROOT}/ERP015684/misc/tmp/ERP015684_by518_supplementary_data_s8.pdf "${v}"

