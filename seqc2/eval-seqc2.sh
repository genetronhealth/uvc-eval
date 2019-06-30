#!/usr/bin/env bash
{
set -evx

SCRIPTDIR=$(dirname $(which $0))

if [ -z "${EVALROOT}" ]; then
    source "${SCRIPTDIR}/../source-eval.sh"
else
    source "${EVALROOT}/source-eval.sh"
fi

TRUTH_BED="${EVALROOT}/seqc2/truth/High-Confidence_Regions.bed"
TRUTH_VCF="${EVALROOT}/seqc2/truth/high-confidence-sALL_in_HC-regions_v1.1.vcf.gz"
TRUTH_BOTHCONF_VCF="${EVALROOT}/seqc2/truth/high-confidence-sALL_in_HC-regions_v1.1.BothConf.vcf.gz"
TRUTH_HIGHCONF_VCF="${EVALROOT}/seqc2/truth/high-confidence-sALL_in_HC-regions_v1.1.HighConf.vcf.gz"

if isfound "${PAT}" "run-truth-prep" ; then
    TRUTH_VCF_SNV="${EVALROOT}/seqc2/truth/high-confidence-sSNV_in_HC-regions_v1.1.vcf.gz"
    TRUTH_VCF_INDEL="${EVALROOT}/seqc2/truth/high-confidence-sINDEL_in_HC-regions_v1.1.vcf.gz"
    bcftools index -ft "${TRUTH_VCF_SNV}" ; bcftools index -ft "${TRUTH_VCF_INDEL}"
    # Definitions for GT and GQ are from https://samtools.github.io/hts-specs/VCFv4.3.pdf
    bcftools concat -a "${TRUTH_VCF_SNV}" "${TRUTH_VCF_INDEL}" | awk -F"\t" \
            'OFS="\t" {if ($0 ~ "^##") {print $0;} else if ($0 ~ "^#CHROM") {print $0 "\tFORMAT\tHCC1395"} else {print $0 "\tGT:GQ\t0/1:50";} }' \
            | sed 's/##contig=<ID=chr1,length=248956422>/##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n##contig=<ID=chr1,length=248956422>/g' \
            |bcftools view -Oz -o "${TRUTH_VCF}"
    bcftools index -ft "${TRUTH_VCF}"
    
    setpass='OFS="\t" {if ($0 ~ "^#") {} else {$7 = "PASS"} ; print $0 }' 

    bcftools view               "${TRUTH_VCF}" | awk -F"\t" "${setpass}" | bcftools view -Oz -o "${TRUTH_BOTHCONF_VCF}" -
    bcftools index -ft "${TRUTH_BOTHCONF_VCF}"
    bcftools view -f "HighConf" "${TRUTH_VCF}" | awk -F"\t" "${setpass}" | bcftools view -Oz -o "${TRUTH_HIGHCONF_VCF}" -
    bcftools index -ft "${TRUTH_HIGHCONF_VCF}"
fi

if isfound "${PAT}" use-wes-data ; then
    fqdir="${EVALROOT}/seqc2/wes-fq"
    bamdir="${EVALROOT}/seqc2/wes-bam"
    vcfdir="${EVALROOT}/seqc2/wes-vcf"
elif isfound "${PAT}" use-wgs-data ; then
    fqdir="${EVALROOT}/seqc2/wgs-fq"
    bamdir="${EVALROOT}/seqc2/wgs-bam"
    vcfdir="${EVALROOT}/seqc2/wgs-vcf"
else
    echo "Either use-wes-data or use-wgs-data should be set"
    exit 1
fi

${mkdir777} -p "${fqdir}" "${bamdir}" "${vcfdir}"

SRRS=$(ls "${fqdir}/SRR"*.fastq.gz | awk -F"/|_" '{print $(NF-1)}' |sort -n | uniq | grep -P "${PAT}")
for SRR in $SRRS; do # start of big for-loop

bam="${bamdir}/${SRR}.bam"
bed="${bamdir}/${SRR}.bed"
depth="${bamdir}/${SRR}.depth"

tSraSample=$(cat "${EVALROOT}/seqc2/SraRunTable_SRP162370.txt" | grep ",breast carcinoma," | grep "${SRR}" | awk -F"," '{print $(NF-12)}')
nSraSample="${tSraSample/_T_/_N_}" 

TREFMAT="HCC1395"
NREFMAT="HCC1395BL"
if [ -z ${tSraSample} ]; then
    refmat="${NREFMAT}"
else
    refmat="${TREFMAT}"
fi

if isfound "${PAT}" "run-bwa" ; then
    r2="${fqdir}/${SRR}_2.fastq.gz"
    platform=ILLUMINA
    if [ ! -f "${r2}" ]; then
        r2=""
        platform=IONTORRENT
    fi
    rgline="@RG\tID:${SRR}\tSM:${refmat}\tLB:${SRR}\tPL:${platform}\tPU:${SRR}"
    "${BWA}" mem -t ${ncpus} -R "${rgline}" "${GRCH38}" "${fqdir}/${SRR}_1.fastq.gz" ${r2} \
        | samtools view -@ "${nbams}" -bh - | samtools sort -@ "${nbams}" -o "${bam}" -
    samtools index -@ ${ncpus} "${bam}"
fi

if isfound "${PAT}" "run-one-sample-bed-prep" ; then
    for depth_thres in 50 100 200; do
        samtools depth "${bam}" -b "${TRUTH_BED}" | awk "\$3 >= ${depth_thres}" | awk -F"\t" 'OFS="\t" {print $1, $2, $2+1}' \
            | bedtools merge -i - > "${bamdir}/${SRR}.${depth_thres}.bed" # "${bed}"
    done
fi

if [ -z ${tSraSample} ]; then
    echo "Skipping the tumor sample ${tSraSample} as it is not found in ${EVALROOT}/seqc2/SraRunTable_SRP162370.txt"
    continue
fi
nSRR=$(cat "${EVALROOT}/seqc2/SraRunTable_SRP162370.txt" | grep ",B lymphoblast," | grep "${nSraSample}" | awk -F"," '{print $(1)}')
if [ -z "${nSRR}" ]; then
    echo "Skipping the tumor sample ${tSraSample} as there is no matched normal sample ${nSraSample}"
    continue
fi

tbam="${bamdir}/${SRR}.bam" 
if [ ! -f "${bamdir}/${SRR}.bam" ]; then
    echo "Skipping the tumor sample ${tSraSample} as there is no tbam file ${bamdir}/${SRR}.bam"
    continue	
fi

nbam="${bamdir}/${nSRR}.bam" 
if [ ! -f "${bamdir}/${nSRR}.bam" ]; then
    echo "skipping the normal sample ${nSraSample} as there is no nbam file ${bamdir}/${nSRR}.bam"
    continue
fi

tExperiment=$(cat "${EVALROOT}/seqc2/SraRunTable_SRP162370.txt" | grep "^${SRR}," | awk -F"," '{print $(NF-10)}' | sed 's/ /-/g')
nExperiment=$(cat "${EVALROOT}/seqc2/SraRunTable_SRP162370.txt" | grep "^${nSRR}," | awk -F"," '{print $(NF-10)}' | sed 's/ /-/g')
tAssayType=$(cat "${EVALROOT}/seqc2/SraRunTable_SRP162370.txt" | grep "^${SRR}," | awk -F"," '{print $(3)}')
nAssayType=$(cat "${EVALROOT}/seqc2/SraRunTable_SRP162370.txt" | grep "^${nSRR}," | awk -F"," '{print $(3)}')
tPM=$(cat "${EVALROOT}/seqc2/SraRunTable_SRP162370.txt" | grep "^${SRR}," | awk -F"," '{print $(NF-14)}' | sed 's/ /-/g')
nPM=$(cat "${EVALROOT}/seqc2/SraRunTable_SRP162370.txt" | grep "^${nSRR}," | awk -F"," '{print $(NF-14)}' | sed 's/ /-/g')

if isfound "${PAT}" "run-tn-bed-prep" ; then
    for tdepth in 100 200; do
        for ndepth in 50 100; do
            evalbed="${bamdir}/${SRR}_${nSRR}.${tdepth}_${ndepth}.bed"
            bedtools intersect -a "${bamdir}/${SRR}.${tdepth}.bed" -b "${bamdir}/${nSRR}.${ndepth}.bed" | bedtools sort -i - | bedtools merge -i - > "${evalbed}"
        done
    done
fi

source "${EVALROOT}/source-tissue.sh"

recaltbam="${tbam/%.bam/.rmdup_recal.bam}"
recalnbam="${nbam/%.bam/.rmdup_recal.bam}"
dindeltbam="${tbam/%.bam/.rmdup_recal_dindel.bam}"
dindelnbam="${nbam/%.bam/.rmdup_recal_dindel.bam}"

run_bamprep "${PAT}" "${GRCH38}" "${tbam}" "${EVALROOT}/datafiles/GRCh38/GATK/All_20180418.vcf.gz" "${ncpus}"
run_bamprep "${PAT}" "${GRCH38}" "${nbam}" "${EVALROOT}/datafiles/GRCh38/GATK/All_20180418.vcf.gz" "${ncpus}"

if isfound "${PAT}" "enable-tn-pair" ; then
    
    is_call_done=0
    for tdepth in 100 200; do for ndepth in 50 100; do for conflevel in ".BothConf" ".HighConf"; do
    
    TRUTH_BED="${bamdir}/${SRR}_${nSRR}.${tdepth}_${ndepth}.bed"
    TRUTH_VCF_GZ="${EVALROOT}/seqc2/truth/high-confidence-sALL_in_HC-regions_v1.1${conflevel}.vcf.gz"
    
    if [ ${is_call_done} -gt 0 ]; then 
        enable_only_vcf_eval="enable-only-vcf-eval"
    else
        enable_only_vcf_eval="enable-not-only-vcf-eval"
    fi
    tnpair_resdir="${vcfdir}/${SRR}_${nSRR}.resdir/"
    tnpair_evaldir="${vcfdir}/${SRR}_${nSRR}.${tdepth}_${ndepth}${conflevel}.evaldir/"
    tnpair_evalsuffix="${tdepth}_${ndepth}${conflevel}"
    
    ${mkdir777} "${tnpair_resdir}"
    run_tnpair "${PAT}.enable-uvc1-all.enable-strelka2-all.enable-lolopicker-all.enable-varscan2-all.enable-somaticsniper-all.${enable_only_vcf_eval}" \
        "${tnpair_resdir},${tnpair_evalsuffix}" "${GRCH38}" "${TRUTH_VCF_GZ}" "${TRUTH_BED}" "${tbam}" "${SRR},${TREFMAT}" "${nbam}" "${nSRR},${NREFMAT}" "${ncpus}"
    run_tnpair "${PAT}.enable-gatk4mutect2-all.${enable_only_vcf_eval}" \
        "${tnpair_resdir},${tnpair_evalsuffix}" "${GRCH38}" "${TRUTH_VCF_GZ}" "${TRUTH_BED}" "${recaltbam}" "${SRR},${TREFMAT}" "${recalnbam}" "${nSRR},${NREFMAT}" "${nbams}"
    run_tnpair "${PAT}.enable-lofreq-all.${enable_only_vcf_eval}" \
        "${tnpair_resdir},${tnpair_evalsuffix}" "${GRCH38}" "${TRUTH_VCF_GZ}" "${TRUTH_BED}" "${dindeltbam}" "${SRR},${TREFMAT}" "${dindelnbam}" "${nSRR},${NREFMAT}" "${ncpus}"
    is_call_done=1
    
    if isfound "${PAT}" "run-custom-prplot"; then
        if [ ${conflevel} = ".BothConf" ]; then
            conftype="keep_both_medium_and_high_confidence_somatic_calls"
        else
            conftype="keep_only_high_confidence_somatic_calls"
        fi
        headerline1="tumor_cell_line=HCC1395 normal_cell_line=HCC1395BL"
        headerline2="tumor=${SRR}/${tSraSample} normal=${nSRR}/${nSraSample}"
        headerline3="tumor_assay=${tAssayType}/${tExperiment}/${tPM} normal_assay=${nAssayType}/${nExperiment}/${nPM}"
        headerline4="tumor_and_normal_regions_of_interest_min_sequencing_depths=${tdepth}x_and_${ndepth}x"
        headerline5="truth_generation=${conftype}"
        prplotHeader="${headerline1}\n${headerline2}\n${headerline3}\n${headerline4}\n${headerline5}"
        tsample="${SRR}"
        nsample="${nSRR}"
        resdir="${tnpair_resdir}"
        mkdir -p "${resdir}"
        eval_uvc="${resdir}/${tsample}_${nsample}.${PARAMSET}_norm_vcfeval-all${tnpair_evalsuffix}.dir"
        eval_mutect2raw="${resdir}/${tsample}_${nsample}.gatk4mutect2_norm_vcfeval-all${tnpair_evalsuffix}.dir"
        eval_mutect2filt="${resdir}/${tsample}_${nsample}.gatk4mutect2_norm_vcfeval-filt${tnpair_evalsuffix}.dir"
        eval_strelka2="${resdir}/${tsample}_${nsample}.strelka2_norm_vcfeval-all${tnpair_evalsuffix}.dir"
        eval_varscan2raw="${resdir}/${tsample}_${nsample}.varscan2_norm_vcfeval-all${tnpair_evalsuffix}.dir"
        eval_varscan2filt="${resdir}/${tsample}_${nsample}.varscan2_norm_vcfeval-filt${tnpair_evalsuffix}.dir"
        eval_somaticsniper="${resdir}/${tsample}_${nsample}.somaticsniper_norm_vcfeval-all${tnpair_evalsuffix}.dir"
        eval_lofreq="${resdir}/${tsample}_${nsample}.lofreq_norm_vcfeval-all${tnpair_evalsuffix}.dir"
        eval_lolopicker="${resdir}/${tsample}_${nsample}.lolopicker_norm_vcfeval-all${tnpair_evalsuffix}.dir"
        
        rocfile="/snp_roc.tsv.gz"
        python "${EVALROOT}/common/custom-prplot.py" "${resdir}/${tsample}_${nsample}_${tnpair_evalsuffix}_all-methods-custom-prplot_snp.pdf"     "${prplotHeader}" 0.00,0.00 \
            "${eval_uvc}${rocfile}" \
            "${eval_mutect2raw}${rocfile}" \
            "${eval_mutect2filt}${rocfile}" \
            "${eval_strelka2}${rocfile}" \
            "${eval_varscan2raw}${rocfile}" \
            "${eval_varscan2filt}${rocfile}" \
            "${eval_somaticsniper}${rocfile}" \
            "${eval_lofreq}${rocfile}" \
            
            #"${eval_lolopicker}${rocfile}" \

        rocfile="/non_snp_roc.tsv.gz"
        python "${EVALROOT}/common/custom-prplot.py" "${resdir}/${tsample}_${nsample}_${tnpair_evalsuffix}_all-methods-custom-prplot_non_snp.pdf" "${prplotHeader}" 0.00,0.00 \
            "${eval_uvc}${rocfile}" \
            "${eval_mutect2raw}${rocfile}" \
            "${eval_mutect2filt}${rocfile}"\
            "${eval_strelka2}${rocfile}" \
            "${eval_varscan2raw}${rocfile}" \
            "${eval_varscan2filt}${rocfile}" \
            "${eval_somaticsniper}${rocfile}" \
            "${eval_lofreq}${rocfile}" \
            
            #"${eval_lolopicker}${rocfile}" \
    
    fi
    done ; done ; done ;
fi

done # end of big for-loop
}
exit $?

