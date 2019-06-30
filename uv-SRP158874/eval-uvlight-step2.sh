SCRIPTDIR=$(dirname $(which $0))
if [ -z "${EVALROOT}" ]; then
    EVALROOT=$(cd "${SCRIPTDIR}/../" && pwd)
fi

REF="${EVALROOT}/datafiles/hs37d5.fa"

# 62 is from https://onlinelibrary.wiley.com/doi/pdf/10.1111/ahg.12114, which says that triallelic site in exon is very unlikely.
flags=" -q 0 " # --tn-syserr-norm-devqual 0 --germ-phred-het3al-snp 62 " # --is-tumor-format-retrieved 1 --syserr-norm-devqual 0 "
outdir="${EVALROOT}/uv-SRP158874/"
# tumor normal reversed
echo "${EVALROOT}/tools/uvc/uvcTN.sh" $REF SRR7757437.bam SRR7757438.bam ${outdir}/SRR7757437_SRR7757438_TN SRR7757437_SRR7757438_TN ${flags}
echo "${EVALROOT}/tools/uvc/uvcTN.sh" $REF SRR7757439.bam SRR7757440.bam ${outdir}/SRR7757439_SRR7757440_TN SRR7757439_SRR7757440_TN ${flags}
echo "${EVALROOT}/tools/uvc/uvcTN.sh" $REF SRR7757441.bam SRR7757442.bam ${outdir}/SRR7757441_SRR7757442_TN SRR7757441_SRR7757442_TN ${flags}
echo "${EVALROOT}/tools/uvc/uvcTN.sh" $REF SRR7757443.bam SRR7757444.bam ${outdir}/SRR7757443_SRR7757444_TN SRR7757443_SRR7757444_TN ${flags}
# tumor normal non-reversed
echo "${EVALROOT}/tools/uvc/uvcTN.sh" $REF SRR7757438.bam SRR7757437.bam ${outdir}/SRR7757437_SRR7757438_NT SRR7757437_SRR7757438_NT ${flags}
echo "${EVALROOT}/tools/uvc/uvcTN.sh" $REF SRR7757440.bam SRR7757439.bam ${outdir}/SRR7757439_SRR7757440_NT SRR7757439_SRR7757440_NT ${flags}
echo "${EVALROOT}/tools/uvc/uvcTN.sh" $REF SRR7757442.bam SRR7757441.bam ${outdir}/SRR7757441_SRR7757442_NT SRR7757441_SRR7757442_NT ${flags}
echo "${EVALROOT}/tools/uvc/uvcTN.sh" $REF SRR7757444.bam SRR7757443.bam ${outdir}/SRR7757443_SRR7757444_NT SRR7757443_SRR7757444_NT ${flags}
