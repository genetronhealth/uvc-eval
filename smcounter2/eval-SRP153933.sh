#!/usr/bin/env bash
# need GNU parallel, java8, fastq-dump (with aspera), bwa, samtools, bcftools, and vcfeval (from RTG)
# SOFT, REF (b37, 1000-genome version of GRCh37), and HG19 are environment variables



SCRIPTDIR=$(dirname $(which $0))
source "${SCRIPTDIR}/../source-eval.sh"

UMI_FILTER_EXPR_1=" cVQ2M[:0] - cVQ2[:1] == 0 && QUAL >= 40 "
UMI_FILTER_EXPR_2="((APXM[:0] - (30 * APDP[:0])) <= 0 || TYPE != \"snp\") && (TYPE = \"snp\" || cMmQ[:1] >= 20)"
#UMI_FILTER_EXPR_2="((APXM[:0] - (30 * APDP[:0])) <= 0 || TYPE != \"snp\")"

SOFTNAME=smcounter2
SRP=SRP153933

if [ $(echo "${PAT}" | grep -P "SRR7526729|N13532" | wc -l) -gt 0 ]; then
    
    SRA=SRR7526729
    SAMPLE=N13532
    ROI=${SCRIPTDIR}/google-download/beds/CDHS-13532Z-10181.roi.bed # HG19
    HOMVCF_NOHDR="${EVALROOT}/smcounter2/smcounter-v2-paper/N13532/misc/na24385.plus.na12878.hom.all.vcf"
    HOMVCF="${EVALROOT}/smcounter2/smcounter-v2-paper/N13532/misc/na24385.plus.na12878.hom.all.hasheader.vcf"
    HETVCF="${EVALROOT}/smcounter2/smcounter-v2-paper/N13532/misc/na12878.uniq.all.het.vcf"
    TRUTHVCF="${HETVCF}" #"${EVALROOT}/smcounter2/smcounter-v2-paper/N13532/misc/truth.vcf"
    TRUTHVCFGZ="${TRUTHVCF}.gz"
    TRUTHBED="${EVALROOT}/smcounter2/smcounter-v2-paper/N13532/misc/N13532.hc.all.na12878.na24385.v3.3.2.bed" 
    TRUTHBED_HG19="${EVALROOT}/smcounter2/smcounter-v2-paper/N13532/misc/N13532.hc.all.na12878.na24385.v3.3.2.hg19.bed" 
    GERM_VARS="${HOMVCF}.gz"
    SOMA_VARS="${TRUTHVCFGZ}"
    CONF_REGIONS="${TRUTHBED_HG19}"
    
    if [ $(echo "${PAT}" | grep "run-dataprep" | wc -l) -gt 0 ]; then
        bcftools view "${HETVCF}.gz" | grep "^#" > "${HOMVCF}"
        cat "${HOMVCF_NOHDR}" | sed 's/^/chr/g' >> "${HOMVCF}"
        bcftools view -Oz -o "${HOMVCF}.gz" "${HOMVCF}"
        bcftools index -f -t "${HOMVCF}.gz"
        bcftools index -f -t "${TRUTHVCFGZ}"
        cat "${TRUTHBED}" | sed 's/^/chr/g' > "${TRUTHBED_HG19}"
    fi
elif [ $(echo "${PAT}" | grep -P "SRR7526728|N0261" | wc -l) -gt 0 ]; then
    SRA=SRR7526728
    SAMPLE=N0261
    ROI=${SCRIPTDIR}/google-download/beds/CDHS-13907Z-562.roi.bed # HG19
fi

DATADIR="${SCRIPTDIR}/data/${SRP}/"

CONF_REGIONS="${DATADIR}/smcounter2dir/${SAMPLE}-${SRA}.truth-regions-hg19.bed"
SOMA_VARS="${DATADIR}/smcounter2dir/${SAMPLE}-${SRA}.truth-somatic-hg19.vcf.gz"
GERM_VARS="${DATADIR}/smcounter2dir/${SAMPLE}-${SRA}.truth-germline-hg19.vcf.gz"

if [ $(echo "${PAT}" | grep run-dataprep | wc -l) -gt 0 ]; then
    bedtools intersect -a "${HG001BED}" -b "${HG002BED}" > "${DATADIR}/smcounter2dir/HG001_HG002_highconf.bed"
    CONF_ALL_REGIONS="${DATADIR}/smcounter2dir/HG001_HG002_highconf-hg19.bed"
    cat "${DATADIR}/smcounter2dir/HG001_HG002_highconf.bed" | awk 'OFS="\t"; {print "chr"$0}' > "${CONF_ALL_REGIONS}"
    bedtools intersect -a "${CONF_ALL_REGIONS}" -b "${ROI}" > "${CONF_REGIONS}"
    
    bcftools view "${HG002VCF}" | awk 'OFS="\t" { if ($0 ~ "^#") {print;} else {print "chr"$0;} }' | sed 's/##contig=<ID=/##contig=<ID=chr/g' | bcftools view -T "${CONF_REGIONS}" -Oz -o "${GERM_VARS}"
    bcftools index -f -t "${GERM_VARS}"
    
    bcftools view "${HG001VCF}" | awk 'OFS="\t" { if ($0 ~ "^#") {print;} else {print "chr"$0;} }' | sed 's/##contig=<ID=/##contig=<ID=chr/g' | bcftools view -T "${CONF_REGIONS}" -Oz -o "${SOMA_VARS}.overlap-germline.vcf.gz"
    bcftools index -f -t "${SOMA_VARS}.overlap-germline.vcf.gz"
    
    bcftools isec ${GROUNDTRUTH_ISEC_FLAGS} -Oz -p "${GERM_VARS}.isecdir/"  "${GERM_VARS}" "${SOMA_VARS}.overlap-germline.vcf.gz"
    bcftools norm ${GROUNDTRUTH_NORM_FLAGS} -Oz -o "${GERM_VARS}.isecdir/0001-norm.vcf.gz" "${GERM_VARS}.isecdir/0001.vcf.gz"
    cp "${GERM_VARS}.isecdir/0001-norm.vcf.gz" "${SOMA_VARS}"
    bcftools index -f -t "${SOMA_VARS}"
fi

SRADIR="${EVALROOT}/smcounter2/data/${SRP}/sra/" 
mkdir -p "${SRADIR}"
FQDIR="${EVALROOT}/smcounter2/data/${SRP}/fq/"
mkdir -p "${FQDIR}"

for sra in $(ls "${SRADIR}" | grep ".sra$" | grep -P "${PAT}") ; do
    if [ $(echo "${PAT}" | grep run-fastq-dump | wc -l) -gt 0 ]; then
        date && time -p "${fastqdump}" --gzip --split-files -O "${FQDIR}" "${SRADIR}/${sra}"
    fi
done

if [ $(echo "${PAT}" | grep run-extract-barcode | wc -l) -gt 0 ]; then
    for fq1 in $(ls "${FQDIR}" | grep "_1.fastq.gz$" | grep -P "${PAT}") ; do
        fq2=${fq1/%_1.fastq.gz/_2.fastq.gz}
        out1=${fq1/%_1.fastq.gz/_1_barcoded.fastq.gz}
        out2=${fq2/%_2.fastq.gz/_2_barcoded.fastq.gz}
        if [ $(echo "${PAT}" | grep run-extract-barcode | wc -l) -gt 0 ]; then
            echo time -p "${UVCBINDIR}/extract-barcodes.py" "${FQDIR}/${fq2}" "${FQDIR}/${fq1}" "${FQDIR}/${out2}" "${FQDIR}/${out1}" 0 12 # 0 11 has one missing base
        fi
    done | parallel -j 8
fi

BAMDIR="${EVALROOT}/${SOFTNAME}/data/${SRP}/bam/"
mkdir -p -m uga+rwx ${BAMDIR}
for fq1 in $(ls "${FQDIR}" | grep "_1_barcoded.fastq.gz$" | grep -P "${PAT}") ; do
    fq2=${fq1/%_1_barcoded.fastq.gz/_2_barcoded.fastq.gz}
    bam=${fq1/%_1_barcoded.fastq.gz/_12.bam}
    if [ $(echo "${PAT}" | grep run-bwa | wc -l) -gt 0 ]; then
        date && time -p "${BWA}" mem -t 32 "${HG19}" "${FQDIR}/${fq1}" "${FQDIR}/${fq2}" | samtools view -bh1 | samtools sort -@2 -o "${BAMDIR}/${bam}"
        samtools index "${BAMDIR}/${bam}"
    fi
done

BASELINEDIR="${EVALROOT}/${SOFTNAME}/data/${SRP}/smcounter2dir/"
mkdir -p -m uga+rwx "${BASELINEDIR}"
for fq1 in $(ls "${FQDIR}" | grep "_1_barcoded.fastq.gz$" | grep -P "${PAT}") ; do
    fq2=${fq1/%_1_barcoded.fastq.gz/_2_barcoded.fastq.gz}
    resultdir=${fq1/%_1_barcoded.fastq.gz/.resultdir}
    ## IMPORTANT: must run smcounter2 in docker here
done

IMBAperc=180 #200
IMBAfloat=1.8
IEXPR="(bDP >= 10000 && cDP >= 3000 & bFA >= 0.001 && cFA >= 0.001 && QUAL>=65 && MAX(bSB1) <= 250 && MAX(cSB1) <= 250 && MAX(aDB) <= 150 && MAX(bPBL) <= 250 && MAX(bPBR) <= 250 && MAX(cPBL) <= 250 && MAX(cPBR) <= 250 && MAX(cVQ3) >= 100 && MAX(bMMB) <= 250 && MAX(cMMB) <= 250)" 
IEXPR="(((1.0 - cFA - cFR)/ cFA) <= 0.75 && (cFA/bFA) <= 1.75 && MAX(bSB1) <= $IMBAperc && MAX(cSB1) <= $IMBAperc && MAX(bPBL) <= $IMBAperc && MAX(bPBR) <= $IMBAperc && MAX(cPBL) <= $IMBAperc && MAX(cPBR) <= $IMBAperc && MAX(bMMB) <= $IMBAperc && MAX(cMMB) <= $IMBAperc)"

VCFDIR="${EVALROOT}/${SOFTNAME}/data/${SRP}/vcf/${Tparamset}/"
otherVCFDIR="${EVALROOT}/${SOFTNAME}/data/${SRP}/vcf/"
mkdir -p -m uga+rwx "${VCFDIR}"
for bam in $(ls ${BAMDIR} | grep ".bam$" | grep -P "${PAT}"); do
    vcf=${bam/%.bam/_uvc1.vcf.gz}
    prepvcf=${bam/%.bam/_uvc1_prep.vcf.gz}
    if [ $(echo "${PAT}" | grep run-uvc1 | wc -l) -gt 0 ]; then
        date && time -p "${UVC}" -t ${ncpus} -f "${HG19}" -o "${VCFDIR}/${vcf}" -s "${bam}" "${BAMDIR}/${bam}" "${@:3}" --primerlen 23 2> "${VCFDIR}/${vcf}.stderr"
        date && time -p bcftools index -f -t "${VCFDIR}/${vcf}" --threads "${nbams}"
    fi
    ##
    isecdir="${VCFDIR}/${vcf/%.vcf.gz/.germline2call_isec.dir}"
    if [ $(echo "${PAT}" | grep -P "run-uvc1|eval-uvc1" | wc -l) -gt 0 ]; then
        date && time -p bcftools view "${VCFDIR}/${vcf}" -Oz -o "${VCFDIR}/${prepvcf}" -i 'QUAL >= 15 && ALT != "*"' --threads "${nbams}"
        date && time -p bcftools index -f -t "${VCFDIR}/${prepvcf}" --threads "${nbams}"
        date && time -p bcftools isec ${GROUNDTRUTH_ISEC_FLAGS} -Oz -p "${isecdir}" "${GERM_VARS}" "${VCFDIR}/${prepvcf}" -R "${CONF_REGIONS}" --threads "${nbams}"
    fi
    if [ $(echo "${PAT}" | grep eval-uvc1 | wc -l) -gt 0 ]; then
        #bcftools sort -Oz -o "${VCFDIR}/${sortedvcf}" "${VCFDIR}/${vcf}"
        #bcftools index -f -t "${VCFDIR}/${sortedvcf}"
        date && time -p bcftools index -f -t "${isecdir}/0001.vcf.gz" 
        # | bcftools filter -Ou -m+ -s dupedRawStrandBias  -e "(bSB1[0:0] * bDP1[0:0] + bSB1[0:1] * bDP1[0:1]) / (bDP1[0:0] + bDP1[0:1]) > $IMBAperc" \
        # | bcftools filter -Ou -m+ -s dedupRawStrandBias  -e "(cSB1[0:0] * cDP1[0:0] + cSB1[0:1] * cDP1[0:1]) / (cDP1[0:0] + cDP1[0:1]) > $IMBAperc" \
        # | bcftools filter -Ou -m+ -s dupedTriallelic     -e "1 / bFA - 1 - bFR / bFA > 0.33" \
        # | bcftools filter -Ou -m+ -s dedupTriallelic     -e "1 / cFA - 1 - cFR / cFA > 0.33" \
        # | bcftools filter -Ou -m+ -s lowNonrefAlleleDupedEfficiency -e "((1-cFR)/(1-bFR)) > 1.8" \
        # | bcftools filter -Ou -m+ -s highQT3 -e "MAX(cQT3) / MAX(cQT2) > 1.8" \
        echo NOT USED CODE: '''
        bcftools view -i "QUAL>=-40" "${isecdir}/0001.vcf.gz" \
            | bcftools filter -Ou -m+ -s lowThisAlleleDupedEfficiency   -e "(   cFA /   bFA ) > ${IMBAfloat}" \
            | bcftools filter -Ou -m+ -s lowRestAlleleDupedEfficiency   -e "((1-cFA)/(1-bFA)) > ${IMBAfloat} && (cFA < 0.8 || bFA < 0.8)" \
            | bcftools filter -Ou -m+ -s dupedQualThresBias  -e "MAX(bQT3) / MAX(bQT2) > ${IMBAfloat}" \
            | bcftools filter -Ou -m+ -s dedupQualThresBias  -e "MAX(cQT3) / MAX(cQT2) > ${IMBAfloat}" \
            | bcftools filter -Ou -m+ -s dupedDedupBias      -e "(aDB[0:0] * bDP1[0:0] + aDB[0:1] * bDP1[0:1]) / (bDP1[0:0] + bDP1[0:1]) > $IMBAperc" \
            | bcftools filter -Ou -m+ -s dedupDedupBias      -e "(aDB[0:0] * cDP1[0:0] + aDB[0:1] * cDP1[0:1]) / (cDP1[0:0] + cDP1[0:1]) > $IMBAperc" \
            | bcftools filter -Ou -m+ -s dupedAdjStrandBias  -e "(bSBR[0:0] * bDP1[0:0] + bSBR[0:1] * bDP1[0:1]) / (bDP1[0:0] + bDP1[0:1]) > $IMBAperc" \
            | bcftools filter -Ou -m+ -s dedupAdjStrandBias  -e "(cSBR[0:0] * cDP1[0:0] + cSBR[0:1] * cDP1[0:1]) / (cDP1[0:0] + cDP1[0:1]) > $IMBAperc" \
            | bcftools filter -Ou -m+ -s dupedPositionBias   -e "((bPBL[0:0] * bDP1[0:0] + bPBL[0:1] * bDP1[0:1]) / (bDP1[0:0] + bDP1[0:1]) > $IMBAperc) || ((bPBR[0:0] * bDP1[0:0] + bPBR[0:1] * bDP1[0:1]) / (bDP1[0:0] + bDP1[0:1]) > $IMBAperc)" \
            | bcftools filter -Ou -m+ -s dedupPositionBias   -e "((cPBL[0:0] * cDP1[0:0] + cPBL[0:1] * cDP1[0:1]) / (cDP1[0:0] + cDP1[0:1]) > $IMBAperc) || ((cPBR[0:0] * cDP1[0:0] + cPBR[0:1] * cDP1[0:1]) / (cDP1[0:0] + cDP1[0:1]) > $IMBAperc)" \
            | bcftools filter -Ou -m+ -s dupedMismatchBias   -e "(bMMB[0:0] * bDP1[0:0] + bMMB[0:1] * bDP1[0:1]) / (bDP1[0:0] + bDP1[0:1]) > ($IMBAperc)" \
            | bcftools filter -Ou -m+ -s dedupMismatchBias   -e "(cMMB[0:0] * cDP1[0:0] + cMMB[0:1] * cDP1[0:1]) / (cDP1[0:0] + cDP1[0:1]) > ($IMBAperc)" \
            | bcftools view   -Oz -o "${isecdir}/0001-flt2.vcf.gz" 
        bcftools index -f -t "${isecdir}/0001-flt2.vcf.gz"
        '''
        #bcftools view -i 'FORMAT/FTS !~ "&" && FORMAT/FTS !~ "\(&\|^\)QTD2\(&\|$\)$" && QUAL>33 && (ORAQs[:0] - VAQ < 0.0001) ' "${isecdir}/0001.vcf.gz" -Oz -o "${isecdir}/0001-flt.vcf.gz"
        #     ' cVQ2M[:0] - cVQ2[:1] == 0 && QUAL >= 40 && (30 * (aDPff[0:1] + aDPfr[0:1] + aDPrf[0:1] + aDPrr[0:1])) - aXMp1[0:1] >= 0 ' # 3388 / (0+0+31+39)
        
        bcftools view "${isecdir}/0001.vcf.gz" -Oz -o "${isecdir}/0001-flt1.vcf.gz" -i "${UMI_FILTER_EXPR_1}" --threads "${nbams}"
        bcftools view "${isecdir}/0001.vcf.gz" -Oz -o "${isecdir}/0001-flt2.vcf.gz" -i "${UMI_FILTER_EXPR_1} && ${UMI_FILTER_EXPR_2}" --threads "${nbams}"
        "${UVCNORM}" "${isecdir}/0001.vcf.gz" "${isecdir}/0001-uvcnorm.vcf.gz"
        bcftools view "${isecdir}/0001-uvcnorm.vcf.gz" -Oz -o "${isecdir}/0001-flt3.vcf.gz" -i "${UMI_FILTER_EXPR_2}" --threads "${nbams}"
        
        # -i 'FORMAT/FTS == "PASS" && QUAL > 33'
        bcftools index -f -t "${isecdir}/0001-flt1.vcf.gz" --threads "${nbams}"
        bcftools index -f -t "${isecdir}/0001-flt2.vcf.gz" --threads "${nbams}"
        bcftools index -f -t "${isecdir}/0001-flt3.vcf.gz" --threads "${nbams}"

        echo NOT USED CODE: '''
        evalout="${bam/%.bam/_uvc1.vcfeval-TLODQ1.outdir}"
        rm -r "${VCFDIR}/${evalout}" || true
        date && time -p "${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval ${evalflags} -f INFO.TLODQ1 \
            -b "${SOMA_VARS}" \
            -c "${isecdir}/0001-flt.vcf.gz" \
            -e "${CONF_REGIONS}" \
            -t "${EVALROOT}/datafiles/hg19_UCSC.fa.sdf" \
            -o "${VCFDIR}/${evalout}" || true
        logversion "${isecdir}/0001-flt.vcf.gz" "${VCFDIR}/${evalout}/uvc1"
        '''
        
        evalout="${bam/%.bam/_uvc1.vcfeval.outdir}"
        rm -r "${VCFDIR}/${evalout}" || true
        date && time -p ${vcfeval1soma} -f QUAL \
            -b "${SOMA_VARS}" \
            -c "${isecdir}/0001.vcf.gz" \
            -e "${CONF_REGIONS}" \
            -t "${EVALROOT}/datafiles/hg19_UCSC.fa.sdf" \
            -o "${VCFDIR}/${evalout}" || true
        logversion "${isecdir}/0001.vcf.gz" "${VCFDIR}/${evalout}/uvc1"        
        
        evalout="${bam/%.bam/_uvc1.vcfeval-flt1.outdir.disabled}"
        rm -r "${VCFDIR}/${evalout}" || true
        date && time -p ${vcfeval1soma} -f QUAL \
            -b "${SOMA_VARS}" \
            -c "${isecdir}/0001-flt1.vcf.gz" \
            -e "${CONF_REGIONS}" \
            -t "${EVALROOT}/datafiles/hg19_UCSC.fa.sdf" \
            -o "${VCFDIR}/${evalout}" || true
        logversion "${isecdir}/0001-flt1.vcf.gz" "${VCFDIR}/${evalout}/uvc1"
        
        evalout="${bam/%.bam/_uvc1.vcfeval-filtered.outdir}"
        rm -r "${VCFDIR}/${evalout}" || true
        date && time -p ${vcfeval1soma} -f QUAL \
            -b "${SOMA_VARS}" \
            -c "${isecdir}/0001-flt2.vcf.gz" \
            -e "${CONF_REGIONS}" \
            -t "${EVALROOT}/datafiles/hg19_UCSC.fa.sdf" \
            -o "${VCFDIR}/${evalout}" || true
        logversion "${isecdir}/0001-flt2.vcf.gz" "${VCFDIR}/${evalout}/uvc1"
        
        evalout="${bam/%.bam/_uvc1.any.vcfeval-filtered.outdir}"
        rm -r "${VCFDIR}/${evalout}" || true
        date && time -p ${vcfeval1soma} -f QUAL \
            -b "${SOMA_VARS}" \
            -c "${isecdir}/0001-flt3.vcf.gz" \
            -e "${CONF_REGIONS}" \
            -t "${EVALROOT}/datafiles/hg19_UCSC.fa.sdf" \
            -o "${VCFDIR}/${evalout}" || true
        logversion "${isecdir}/0001-flt3.vcf.gz" "${VCFDIR}/${evalout}/uvc1"
        
        echo Results are at "${VCFDIR}/${evalout}"

        echo '''
        evalout="${bam/%.bam/_uvc1.vcfeval-QUAL.outdir.disabled}"
        rm -r "${VCFDIR}/${evalout}" || true
        date && time -p "${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval ${evalflags} -f INFO.QUAL1 \
            -b "${SOMA_VARS}" \
            -c "${isecdir}/0001.vcf.gz" \
            -e "${CONF_REGIONS}" \
            -t "${EVALROOT}/datafiles/hg19_UCSC.fa.sdf" \
            -o "${VCFDIR}/${evalout}" || true
        logversion "${isecdir}/0001.vcf.gz" "${VCFDIR}/${evalout}/uvc1"        
        
        evalout="${bam/%.bam/_uvc1.vcfeval-filtered-QUAL.outdir.disabled}"
        rm -r "${VCFDIR}/${evalout}" || true
        date && time -p "${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval ${evalflags} -f INFO.QUAL1 \
            -b "${SOMA_VARS}" \
            -c "${isecdir}/0001-flt.vcf.gz" \
            -e "${CONF_REGIONS}" \
            -t "${EVALROOT}/datafiles/hg19_UCSC.fa.sdf" \
            -o "${VCFDIR}/${evalout}" || true
        logversion "${isecdir}/0001-flt.vcf.gz" "${VCFDIR}/${evalout}/uvc1"
        '''
    fi
done

sample="${SAMPLE}-${SRA}"
evalout=${sample}_smcounter2.vcfeval.outdir #"${bam/%.bam/_smcounter2.vcfeval.outdir}"
absvcf="${EVALROOT}/smcounter2/data/${SRP}/smcounter2dir/${SAMPLE}.smCounter.anno.vcf"
absvcf="${EVALROOT}/smcounter2/data/${SRP}/baseline/SRR/${SAMPLE}.smCounter.anno.vcf"
absvcfgz="${absvcf}.gz"
isecdir="${absvcfgz/%.vcf.gz/.germline2call_isec.dir}"

if [ $(echo "${PAT}" | grep eval-smcounter2 | wc -l) -gt 0 ]; then
    bcftools view "${absvcf}" -Oz -o "${absvcfgz}"
    bcftools index -f -t "${absvcfgz}"
    bcftools isec ${GROUNDTRUTH_ISEC_FLAGS} -Oz -p "${isecdir}" "${GERM_VARS}" "${absvcfgz}" -R "${CONF_REGIONS}"
    bcftools index -f -t "${isecdir}/0001.vcf.gz"
    rm -r "${otherVCFDIR}/${evalout}" || true 
    date && time -p ${vcfeval1soma} -f QUAL \
        -b "${SOMA_VARS}" \
        -c "${isecdir}/0001.vcf.gz" \
        -e "${CONF_REGIONS}" \
        -t "${EVALROOT}/datafiles/hg19_UCSC.fa.sdf" \
        -o "${otherVCFDIR}/${evalout}" || true
fi

for samplesra in N0261.SRR7526728 N13532.SRR7526729; do
    sample=$(echo $samplesra | awk -F"." '{print $1}')
    sra=$(   echo $samplesra | awk -F"." '{print $2}')
    if [ $(echo "${PAT}" | grep run-custom-prplot | grep -P "${sample}|${sra}" | wc -l) -gt 0 ]; then
        mkdir -p "${VCFDIR}/rocplots/" || true
        
        keysuffix=custom-prplot
        plothdr="VAF=0.005_and_0.010 sample=${sample} accession=${sra};notchinese-notbest"
        if [ $(echo "${PAT}" | grep patent-plot | wc -l) -gt 0 ]; then
            keysuffix=patent-plot
            plothdr="变异等为基因占比=0.005和0.010 样本名称=${sample} 编录号=${sra};ischinese-isbest"
        fi
        
        python "${EVALROOT}/common/custom-prplot.py" "${VCFDIR}/rocplots/${sample}-${sra}_all_methods_${keysuffix}_snp.pdf" "${plothdr}" 0,0 \
        $(ls "${VCFDIR}/"*${sra}*.any.vcfeval*".outdir/snp_roc.tsv.gz" | sort -r) $(ls "${otherVCFDIR}/"*${sra}*smcounter2.vcfeval*".outdir/snp_roc.tsv.gz" | sort -r)
        
        python "${EVALROOT}/common/custom-prplot.py" "${VCFDIR}/rocplots/${sample}-${sra}_all_methods_${keysuffix}_non_snp.pdf" "${plothdr}" 0,0 \
        $(ls "${VCFDIR}/"*${sra}*.any.vcfeval*".outdir/non_snp_roc.tsv.gz" | sort -r) $(ls "${otherVCFDIR}/"*${sra}*smcounter2.vcfeval*".outdir/non_snp_roc.tsv.gz" | sort -r)
    fi
done
