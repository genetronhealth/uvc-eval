#!/usr/bin/env bash
{
if [ $# -lt 2 ]; then
    printf "Usage: $0 <options> <output-subdirectory>, where <options> is a '/'-separated list of options.\n"
    printf "The <options> consists of a combination of the following terms.\n"
    printf "  use-tname-<tname> (<tname> is 1 22 or all) : use chromosome names, default is all.\n"
    printf "  use-refmat-HG00<x> (<x> is 1 2 or 5) : use the human genome reference material HG001 HG002 or HG005.\n"
    printf "  use-avgdep-<x>x (<x> is 30 or 60) : use 30x or 60x average depth.\n"
    printf "  run-bwa-bgi-pe100 run-bwa-mgi-pe150 : run BWA for BGI and MGI fastqs.\n"
    printf "  use-platform-bgi-pe100 use-platform-mgi-pe150 : "
    printf "    use BGI and MGI platform bam files instead of Illumina platform bam files (Illumina is the default platform).\n"
    printf "  run-gatk4rmdup, run-gatk4bqsr : run gatk4 removal of duplicated reads, run gatk4 BQSR.\n"
    printf "  run-<vc>-<stage> : run the variant caller <vc> at stage <stage>. \n"
    printf "    <vc> is either uvc1, gatk4hc, strelka2, freebayes, or bcftools. \n"
    printf "    <stage> is either call, norm, eval, or all. \n"
    exit 1
fi

set -evx

SCRIPTDIR=$(dirname $(which $0))
if [ -z "${EVALROOT}" ]; then
    source "${SCRIPTDIR}/../source-eval.sh"
else
    source "${EVALROOT}/source-eval.sh"
fi

#evalflags=" --ref-overlap --decompose --all-records --threads=${ncpus}" #" --all-records " #" --all-records --ref-overlap "
#evalflags=" --threads=${ncpus} --decompose --ref-overlap " #" --all-records " #" --all-records --ref-overlap "
VCREF="${HS37D5}"

BWA_BGI_GERM_PARAMS="" # " -A 2 -B 8 -O 20 -E 1 -L 20 "

chrnum=999
hg=HG001
avgdep=30
seqplat="NONE" # IL is Illumina, MS is MGISEQ, and BI is BGISEQ
seqplat_fullname="NONE"
aligner_fullname="NONE"

if isfound "${PAT}" "use-tname-all"; then
    chrnum=999
elif isfound "${PAT}" "use-tname-1"; then
    chrnum=1
elif isfound "${PAT}" "use-tname-22"; then
    chrnum=22
fi

if isfound "${PAT}" "use-refmat-hg001"; then
    hg=HG001
elif isfound "${PAT}" "use-refmat-hg002"; then
    hg=HG002
elif isfound "${PAT}" "use-refmat-hg005"; then
    hg=HG005
fi

if isfound "${PAT}" "use-avgdep-30x"; then
    avgdep=30
elif isfound "${PAT}" "use-avgdep-60x"; then
    avgdep=60
elif isfound "${PAT}" "use-avgdep-300x"; then
    avgdep=300
fi

suf=nochr${chrnum}
AVGDEP="${avgdep}x"

bam="${EVALROOT}/germline/bam/${hg}.hs37d5.${AVGDEP}.${suf}.leftaligned.bam"
#bam="${EVALROOT}/germline/bam/${hg}.hs37d5.${AVGDEP}.${suf}.bam"
#resrootdir="${EVALROOT}/germline/vcf/${PARAMSET}/"
#dir="${EVALROOT}/germline/vcf/${PARAMSET}/${hg}.hs37d5.${AVGDEP}.${suf}.uvc.resdir"

if isfound "${PAT}" "use-realign-abra2"; then
    bam="${EVALROOT}/germline/bam/${hg}.hs37d5.${AVGDEP}.${suf}.abra2.bam"
    #dir="${EVALROOT}/germline/vcf/${hg}.hs37d5.${AVGDEP}.${suf}.uvc-abra2.resdir"
fi

if [ "${hg}" == "HG001" ]; then
    TRUTH_BED="${HG001BED}"
    TRUTH_VCF="${HG001VCF}"
fi
if [ "${hg}" == "HG002" ]; then
    TRUTH_BED="${HG002BED}"
    TRUTH_VCF="${HG002VCF}"
fi
if [ "${hg}" == "HG005" ]; then
    TRUTH_BED="${HG005BED}"
    TRUTH_VCF="${HG005VCF}"
fi

${mkdir777} -p "${EVALROOT}/germline/bed/"
truth_bed="${EVALROOT}/germline/bed/${hg}-${suf}.bed"

if isfound "${PAT}" "run-prep-bed" || [ ! -f "${truth_bed}" ]; then
    if [ ${chrnum} -lt 25 ]; then
        cat "${TRUTH_BED}" | grep -P "^${chrnum}\t" > "${truth_bed}"
    else
        cp "${TRUTH_BED}" "${truth_bed}"
    fi
fi

lane1ord=L01
lane2ord=L02

if isfound "${PAT}" "use-platform-bgi-pe100"; then
    seqplat=BS
    seqplat_fullname="BGI/BGISEQ"
    aligner_fullname="BWA_MEM"
elif isfound "${PAT}" "use-platform-mgi-pe150"; then
    seqplat=MS
    seqplat_fullname="BGI/MGISEQ"
    aligner_fullname="BWA_MEM"
    if [ "${hg}" == "HG002" ]; then
        lane1ord=L03
        lane2ord=L04
    fi
else
    seqplat=IL
    seqplat_fullname="ILLUMINA/HiSeq"
    aligner_fullname="NovoAlign"
fi

if isfound "${PAT}" "run-bwa-bgi-pe100|run-bwa-mgi-pe150"; then
    bamfpref="${hg}${seqplat}30X"
    for lane in "${lane1ord}" "${lane2ord}"; do
        if isfound "${PAT}" "run-bwa-mgi-pe150"; then # MGI
            if [ "${hg}" == "HG001" ]; then
                fq1="${EVALROOT}/germline/fq/MGISEQ2000_PCR-free_NA12878_1_V100003043_${lane}_1.fq.gz"
                fq2="${EVALROOT}/germline/fq/MGISEQ2000_PCR-free_NA12878_1_V100003043_${lane}_2.fq.gz"
                chipUID=V100003043
            elif [ "${hg}" == "HG002" ]; then
                fq1="${EVALROOT}/germline/fq/MGISEQ2000_PCR-free_NA24385_V100002807_${lane}_1.fq.gz"
                fq2="${EVALROOT}/germline/fq/MGISEQ2000_PCR-free_NA24385_V100002807_${lane}_2.fq.gz"
                chipUID=V100002807
            elif [ "${hg}" == "HG005" ]; then
                fq1="${EVALROOT}/germline/fq/MGISEQ2000_PCR-free_NA24631_V100002812_${lane}_1.fq.gz"
                fq2="${EVALROOT}/germline/fq/MGISEQ2000_PCR-free_NA24631_V100002812_${lane}_2.fq.gz"
                chipUID=V100002812
            else
                echo "use-sample-hg00<x> must be specified, where x is either 1 2 or 5."
                exit -1
            fi
        elif isfound "${PAT}" "run-bwa-bgi-pe100"; then # MGI
            if [ "${hg}" == "HG001" ]; then
                fq1="${EVALROOT}/germline/fq/BGISEQ500_PCRfree_NA12878_CL100076243_${lane}_read_1.fq.gz"
                fq2="${EVALROOT}/germline/fq/BGISEQ500_PCRfree_NA12878_CL100076243_${lane}_read_2.fq.gz"
                chipUID=CL100076243
            elif [ "${hg}" == "HG002" ]; then
                fq1="${EVALROOT}/germline/fq/BGISEQ500_PCRfree_NA24385_CL100076190_${lane}_read_1.fq.gz"
                fq2="${EVALROOT}/germline/fq/BGISEQ500_PCRfree_NA24385_CL100076190_${lane}_read_2.fq.gz"
                chipUID=CL100076190
            elif [ "${hg}" == "HG005" ]; then
                fq1="${EVALROOT}/germline/fq/BGISEQ500_PCRfree_NA24631_CL100076244_${lane}_read_1.fq.gz"
                fq2="${EVALROOT}/germline/fq/BGISEQ500_PCRfree_NA24631_CL100076244_${lane}_read_2.fq.gz"
                chipUID=CL100076244
            else
                echo "use-sample-hg00<x> must be specified, where x is either 1 2 or 5."
                exit -1
            fi
        else
            echo "The DNA nanoball (DNB) sequencing type ${seqplat} is not valid (valid values are BS and MS (BGISEQ and MGISEQ))!"
            exit -1
        fi
        sample="${hg}"
        library="${hg}${seqplat}"
        platunit="${chipUID}.${lane}"
        bamfname="${hg}${seqplat}30X${lane}"
        # https://samtools.github.io/hts-specs/SAMv1.pdf : GATK accepts BGI instead of DNBSEQ, the correct PL according to the specs.
        # Hence, use PL:ILLUMINA and PM:MGISEQ as a compromise
        rgline="@RG\tID:${platunit}\tSM:${sample}\tLB:${library}\tPL:ILLUMINA\tPM:MGISEQ\tPU:${platunit}"
        "${BWA}" mem ${BWA_BGI_GERM_PARAMS} -R "${rgline}" -t ${ncpus} "${HS37D5}" "${fq1}" "${fq2}" \
                | samtools sort -@ ${nbams} -o "${EVALROOT}/germline/bam/${bamfname}.bam" -
        samtools index "${EVALROOT}/germline/bam/${bamfname}.bam"
    done
    bamfn60x="${hg}${seqplat}60X"
    samtools merge -@ ${nbams} "${EVALROOT}/germline/bam/${bamfn60x}.bam" \
        "${EVALROOT}/germline/bam/${bamfpref}${lane1ord}.bam" "${EVALROOT}/germline/bam/${bamfpref}${lane2ord}.bam"
    samtools index -@ ${nbams} "${EVALROOT}/germline/bam/${bamfn60x}.bam"
fi

if isfound "${PAT}" "use-platform-bgi-pe100|use-platform-mgi-pe150"; then
    if isfound "${PAT}" "use-avgdep-60x"; then
        bam="${EVALROOT}/germline/bam/${hg}${seqplat}60X.bam"
    elif isfound "${PAT}" "use-avgdep-30x"; then
        bam="${EVALROOT}/germline/bam/${hg}${seqplat}30X${lane1ord}.bam"
    else
        echo "BGI can only have 30x and 60x options"; exit -1
    fi
    if [ "${chrnum}" != "all" -a "${chrnum}" -le 22 ]; then
        oldbam="${bam}"
        bam=$(echo "${oldbam}" | sed "s/.bam$/-${chrnum}.bam/g")
        if isfound "${PAT}" "run-samtools-subsample"; then
            samtools view -@ "${nbams}" -bh "${oldbam}" "${chrnum}" > "${bam}"
            samtools index -@ "${nbams}" "${bam}"
        fi
    fi
else # Illumina
    if [ $(echo "${PAT}" | grep "run-samtools-subsample" | wc -l) -gt 0 ]; then
        if [ "${avgdep}" -lt 300 ]; then
            downfrac=$(python -c "print(float(${avgdep})/300)")
            #samtools view -s ${downfrac} -bh "${EVALROOT}/tn/silico-rawinput/${hg}.hs37d5.300x.bam" ${chrnum} 
            #   | "${EVALROOT}/tools/freebayes-1.3.2/bin/bamleftalign" -f "${VCREF}" > "${bam}"
            if [ ${chrnum} -lt 25 ]; then
                samtools view -@ "${nbams}" -s ${downfrac} -bh "${EVALROOT}/tn/silico-step2input/${hg}.hs37d5.300x.rg-leftaligned.bam" ${chrnum} > "${bam}"
            else
                samtools view -@ "${nbams}" -s ${downfrac} -bh "${EVALROOT}/tn/silico-step2input/${hg}.hs37d5.300x.rg-leftaligned.bam" > "${bam}"
            fi
            samtools index -@ "${nbams}" "${bam}"
        fi
    fi
fi

resultroot="${EVALROOT}/germline/vcf/$(echo ${bam} | awk -F"/" '{print $NF}' | sed 's/\.bam$//g')/"
${mkdir777} -p "${resultroot}"

source "${EVALROOT}/source-germline.sh"
recalbam="${bam/%.bam/.rmdup_recal.bam}"
vcfpref="${hg}${seqplat}${avgdep}X"
chrom="${chrnum}"
run_bamprep "${PAT}" "${HS37D5}" "${bam}" "${EVALROOT}/datafiles/hg19/dbsnp_138.b37.vcf" "${nbams}"

run_germline "${PAT}.enable-uvc1-all.enable-strelka2-all.enable-freebayes-all.enable-bcftools-all" \
        "${resultroot}/" "${HS37D5}" "${TRUTH_VCF}" "${truth_bed}" "${bam}" "${vcfpref},${hg}" "${ncpus}"
run_germline "${PAT}.enable-gatk4hc-all" \
        "${resultroot}/" "${HS37D5}" "${TRUTH_VCF}" "${truth_bed}" "${recalbam}" "${vcfpref},${hg}" "${ncpus}"

if isfound "${PAT}" "run-custom-prplot"; then
    headerline1="cell_line=${hg} sample=${vcfpref}"
    headerline2="average_sequencing_depth=${avgdep}x sequencing_platform=${seqplat_fullname}"
    headerline3="aligner=${aligner_fullname}"
    
    prplotHeader="${headerline1}\n${headerline2}\n${headerline3};disable-baseline-var-num-assertion"
    
    eval_uvc="${resultroot}/${vcfpref}_germline.${PARAMSET}_norm_vcfeval-all.dir"
    eval_gatk4hc="${resultroot}/${vcfpref}_germline.gatk4hc_norm_vcfeval-all.dir"
    eval_strelka2="${resultroot}/${vcfpref}_germline.strelka2_norm_vcfeval-all.dir"
    eval_freebayes="${resultroot}/${vcfpref}_germline.freebayes_norm_vcfeval-all.dir"
    eval_bcftools="${resultroot}/${vcfpref}_germline.bcftools_norm_vcfeval-all.dir"
    
    rocfile="/snp_roc.tsv.gz"
    python "${EVALROOT}/common/custom-prplot.py" "${resultroot}/${vcfpref}_germline_all-methods-custom-prplot_snp.pdf" "${prplotHeader}" 0.00,0.00 \
    "${eval_uvc}${rocfile}" \
    "${eval_gatk4hc}${rocfile}" \
    "${eval_strelka2}${rocfile}" \
    "${eval_freebayes}${rocfile}" \
    "${eval_bcftools}${rocfile}"
    
    rocfile="/non_snp_roc.tsv.gz"
    python "${EVALROOT}/common/custom-prplot.py" "${resultroot}/${vcfpref}_germline_all-methods-custom-prplot_non_snp.pdf" "${prplotHeader}" 0.00,0.00 \
    "${eval_uvc}${rocfile}" \
    "${eval_gatk4hc}${rocfile}" \
    "${eval_strelka2}${rocfile}" \
    "${eval_freebayes}${rocfile}" \
    "${eval_bcftools}${rocfile}"
fi
}
exit $?

