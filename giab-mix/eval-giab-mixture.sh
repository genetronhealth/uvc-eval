#!/usr/bin/env bash
{
if [ $# -lt 2 ]; then
    echo "Usage: $0 optionList outputPrefix"
    echo "  <optionList> is a string consisting of '.'-separated (period-separated or dot-separated) substrings that correspond to options controlling the run. "
    echo "  <outputPrefix> is the sub-directory containing variant-call results for UVC. "
    
    echo "  The options are listed as follows in the order that they should be run. "
    echo "  The user can run a subset of the options given that their pre-requisite options were already successfully run. "
    echo "  In general, the order of commands is: setup, prep (prepare), and run. "
    
    echo "  enable-tumor-only : let variant caller(s) run in tumor-only mode (so far only uvc1, gatk4mutect2, and lofreq can do so). "
    echo "  enable-tn-pair : let variant caller(s) run in tumor-normal-pair mode (all variant callers can do so). "
    
    echo "  use-<type>-mixture : create virtual tumor-normal pairs using <type> original bam data, where <type> can be one of the following three words. "
    echo "    physical : use the standard physical mixture from GIAB. "
    echo "    pa12: use HG001 300x as tumor and HG002 300x as the matched normal. "
    echo "    pa21: use HG002 300x as tumor and HG001 300x as the matched normal. "
    echo "  use-scenario-<sc> : only use the scenario named <sc> for in-silico mixtures. "
    echo "  use-tname-<tname> : only use reads aligned to the chromosome <tname>. "
    
    echo "  step 1 : setup-addreplacerg-<hg> : add read group to the hg (hg is HG001, HG002 or HG005) original novoalign bam file. "
    echo "  step 2 : setup-bamleftalign-<hg> : generate left-aligned bam from the previous bam with rg. "
    echo "  step 3 : prep-truth: prepare truth vcf and bed files for in-silico mixture data depending on the use-scenario-<sc> option. "
    echo "  step 4 : prep-tname-<tname> : use samtools view to retrieve reads aligned to the chromosome <tname>. "
    echo "  step 5 : run-silico-mix: simulate silico tumor normal bams with the given scenario and tnpair. This run option is not applicable to physical mixtures. "
    echo "  step 6 : run-bwa: run bwa mem to map fastq reads in the bam format. "
    echo "      This step is the initial step for physical-mix, and this step is skipped for silico-*-mix. "
    echo "  step 7 : run-gatk4rmdup run-gatk4bqsr run-dindel: three sequential stages for processing bam files. "
    echo "      First, gatk4rmdup removes duplicated reads. "
    echo "      Then, gatk4bqsr performs base-quality-score-recalibration. "
    echo "      Finally, dindel inserts InDel basecall-like qualilites. "
    echo "      These are the preprocessing steps for their variant callers to make variant calls. "
    echo "  step 8 : run-<variantcaller>-<stage> "
    echo "      <variantcaller> can be uvc1, gatk4mutect2, strelka2, lofreq, lolopicker, octopus, lancet, somaticsniper, and varscan2. "
    echo "      <stage> can be call, norm, eval, or all (run all three stages), The order of the three stages is call-norm-eval. "
    echo "      The call stage performs variant calling, the norm stage normalizes variants for evaluation, and the eval stage evaluates calling performance. "
    echo "  step 9 : run-custom-prplot: generate precision-recall curves. "
    exit 1
fi

set -evx

# PAT
SCRIPTDIR=$(dirname $(which $0))
if [ -z "${EVALROOT}" ]; then
    source "${SCRIPTDIR}/../source-eval.sh"
else
    source "${EVALROOT}/source-eval.sh"
fi

MIXRDIR="${EVALROOT}/giab-mix/"

TNAMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
PUPL="PU:L001\tPL:ILLUMINA"
PROCSUFFIX="rg-leftaligned"

if isfound "${PAT}" "test-exit-with-code-1" ; then
    echo "Will exit wtih code 1 to complete the test"
    exit 1
fi

if issubstr "${PAT}" "run-bamformat" ; then
    ${mkdir777} -p "${MIXRDIR}/silico-step2input/"
    for hg in HG001 HG002 HG005; do
        bam1="${MIXRDIR}/silico-rawinput/${hg}.hs37d5.300x.bam"
        bam2="${MIXRDIR}/silico-step2input/${hg}.hs37d5.300x.rg.bam"
        bam3="${MIXRDIR}/silico-step2input/${hg}.hs37d5.300x.rg-leftaligned.bam"
        bam4="${MIXRDIR}/silico-step2input/${hg}.hs37d5.300x.rg-fail.bam"
        
        RGID="${hg}HS300X"
        RGLINE="@RG\tID:${RGID}\tSM:${hg}\tLB:${RGID}\t${PUPL}" 
        if isfound "run-addreplacerg-${hg}" ; then
            time -p $(samtools addreplacerg -@ ${ncpus} -r "${RGLINE}" "${bam1}" | samtools view -@ ${ncpus} -bh > "${bam2}")
            time -p samtools index -@ ${ncpus} "${bam2}"
        fi
        if isfound "run-bamleftalign-${hg}" ; then
            time -p $(samtools view -bh -@ ${ncpus} "${bam2}" | "${EVALROOT}/tools/freebayes-1.3.2/bin/bamleftalign" -f "${HS37D5}" 2> "${bam2/.bam/.stderr}" | \
                   "${EVALROOT}/common/bamfilter.out" - - "${bam4}" 2> "${bam3/.bam/.stderr}" | samtools view -bh -@ ${ncpus}  > "${bam3}")
            
            # LeftAlignIndels seems to run out of memory and crash when started, tried multiple times
            #bam2="${MIXRDIR}/silico-step2input/${hg}.hs37d5.300x.gatk4leftaligned.bam"
            #time ${gatk4lowmem} LeftAlignIndels --input "${bam1}" --reference "${HS37D5}" 2>"${bam2/.bam/.stderr}" | \
            #        samtools addreplacerg -@ ${ncpus} -r "${RGLINE}" - | samtools view -@ ${ncpus} -bh > "${bam2}"
            time -p samtools index -@ ${ncpus} "${bam3}"
        fi        
    done
    exit 0
fi

platform_full=Illumina
aligner_full=BWA
scenario=00
chrom=00
mutect2chromparams=""
if isfound "${PAT}" "use-pa12-mixture|use-pa21-mixture|use-pb12-mixture|use-pb21-mixture" ; then
    for nochr in $(seq 1 22); do
        if isfound "${PAT}" "use-tname-${nochr}" ; then
            chrom=${nochr}
        fi
    done
    if [ $chrom == "00" ]; then echo "chrom cannot be NONE!"; exit -1; fi
    for sc in $(seq -w 1 40) 99; do
        if isfound "${PAT}" use-scenario-${sc} ; then
            scenario=${sc}
        fi
    done
    if [ $scenario == "00" ]; then echo "scenario cannot be NONE!"; exit -1; fi
    mutect2chromparams=" -L ${chrom} "
    
    scline=$(cat "${MIXRDIR}/silico/silico-design.tsv" | awk "\$1 == ${scenario}")
    fracTinT=$(echo $scline | awk '{print $2}') 
    fracTinN=$(echo $scline | awk '{print $3}')
    fracNinT=$(echo $scline | awk '{print $4}')    
    fracNinN=$(echo $scline | awk '{print $5}')
    tpurity=$(echo $scline | awk '{print $8}') 
    npurity=$(echo $scline | awk '{print $9}')
    tdepth=$(echo $scline | awk '{print $10}')    
    ndepth=$(echo $scline | awk '{print $11}')
    if isfound "${PAT}" "use-pa21-mixture|use-pb21-mixture"; then
        mixtype="silico-pa21-mix" 
        silicoTname=HG002
        silicoNname=HG001
        #pn="pl21_${silicoTname}_${silicoNname}_${scenario}"
        pn=HGM21PA${scenario}
        HGnormalBED="${HG001BED}"
        HGnormalVCF="${HG001VCF}"
        HGtumorBED="${HG002BED}"
        HGtumorVCF="${HG002VCF}"
    elif isfound "${PAT}" "use-pa12-mixture|use-pb12-mixture" ; then # let-tn-be-hg002-hg001
        mixtype="silico-pa12-mix"
        silicoTname=HG001
        silicoNname=HG002
        #pn="pl12_${silicoTname}_${silicoNname}_${scenario}"
        pn=HGM12PA${scenario}
        HGnormalBED="${HG002BED}"
        HGnormalVCF="${HG002VCF}"
        HGtumorBED="${HG001BED}"
        HGtumorVCF="${HG001VCF}"
    else
        echo "Runtime error: use-<pa12|pa21|pb12|pb21>-mixture must be found!"
        exit -1
    fi
    if isfound "${PAT}" "use-pb12-mixture|use-pb21-mixture"; then
        mixtype=$(echo "${mixtype}" | sed 's/silico-pa/silico-pb/g')
        pn=$(echo "${pn}" | sed 's/HGM12PA/HGM12PB/g' | sed 's/HGM21PA/HGM21PB/g')
        silicoTbam="${EVALROOT}/germline/bam/${silicoTname}MS60X.bam"
        silicoNbam="${EVALROOT}/germline/bam/${silicoNname}MS60X.bam"
        tdepth=$(python -c "print($tdepth / 5)")
        ndepth=$(python -c "print($ndepth / 5)")
        platform_full=BGI/MGISEQ
        aligner_full=BWA_MEM
    else
        silicoTbam="${MIXRDIR}/silico-step2input/${silicoTname}.hs37d5.300x.${PROCSUFFIX}.bam" 
        silicoNbam="${MIXRDIR}/silico-step2input/${silicoNname}.hs37d5.300x.${PROCSUFFIX}.bam"    
        platform_full=Illumina/HiSeq
        aligner_full=NovoAlign
    fi
    
    nsample=${pn}N
    tsample=${pn}T
    tfname=${pn}T
    nfname=${pn}N
    headerline1="tumor=${tsample}/${silicoTname} normal=${nsample}/${silicoNname}"
    headerline2="tumor_purity=${tpurity} normal_purity=${npurity}"
    headerline3="tumor_average_depth=${tdepth}x normal_average_depth=${ndepth}x"
    headerline4="platform=${platform_full}  mix_type=in_silico"
    headerline5="aligner=${aligner_full} chromosome=${chrom}"
    prplotHeader="${headerline1}\n${headerline2}\n${headerline3}\n${headerline4}\n${headerline5};notchinese-notbest"
    
    if isfound "${PAT}" enable-tumor-only; then
        headerline1="tumor=${tsample}/${silicoTname} normal=NotAvailable"
        headerline2="tumor_purity=${tpurity}"
        headerline3="tumor_average_depth=${tdepth}x"
        headerline4="platform=${platform_full} mix_type=in_silico"
        headerline5="aligner=${aligner_full} chromosome=${chrom}"
        prplotHeader="${headerline1}\n${headerline2}\n${headerline3}\n${headerline4}\n${headerline5};notchinese-notbest"
    fi
    
    NinT="${MIXRDIR}/${mixtype}/bam/tmp/${pn}-${chrom}.NinT.bam"
    NinN="${MIXRDIR}/${mixtype}/bam/tmp/${pn}-${chrom}.NinN.bam"
    TinT="${MIXRDIR}/${mixtype}/bam/tmp/${pn}-${chrom}.TinT.bam"
    TinN="${MIXRDIR}/${mixtype}/bam/tmp/${pn}-${chrom}.TinN.bam"
    NTinN="${MIXRDIR}/${mixtype}/bam/subset/${nfname}-${chrom}.bam"
    NTinT="${MIXRDIR}/${mixtype}/bam/subset/${tfname}-${chrom}.bam"
    if isfound "${PAT}" run-silico-mix ; then
        divbam="${EVALROOT}/common/divbam.out"
        ${mkdir777} -p "${MIXRDIR}/${mixtype}/bam/subset/"
        ${mkdir777} -p "${MIXRDIR}/${mixtype}/bam/tmp/"
        samtools view "${silicoTbam}" -h "${chrom}" -@ ${ncpus} | "${divbam}" - "${fracTinT}" "${TinT}" "${fracTinN}" "${TinN}"
        samtools view "${silicoNbam}" -h "${chrom}" -@ ${ncpus} | "${divbam}" - "${fracNinT}" "${NinT}" "${fracNinN}" "${NinN}"
        samtools merge -u -@ ${ncpus} - "${NinN}" "${TinN}" | samtools reheader -c 'grep -v ^@RG' - | samtools addreplacerg -u -@ ${ncpus} -r "@RG\tID:${nsample}\tSM:${nsample}\tLB:${nsample}:\t${PUPL}" - | samtools view -@ ${ncpus} -bh > "${NTinN}" && samtools index -@ ${ncpus} "${NTinN}"
        samtools merge -u -@ ${ncpus} - "${NinT}" "${TinT}" | samtools reheader -c 'grep -v ^@RG' - | samtools addreplacerg -u -@ ${ncpus} -r "@RG\tID:${tsample}\tSM:${tsample}\tLB:${tsample}:\t${PUPL}" - | samtools view -@ ${ncpus} -bh > "${NTinT}" && samtools index -@ ${ncpus} "${NTinT}"
        exit 0
    fi
elif isfound "${PAT}" "use-physical-mixture"; then
    headerline1="tumor=HG001 normal=HG002"
    headerline2="tumor_purity=0.3 normal_purity=1.0"
    headerline3="tumor_depth=90 normal_depth=30"
    headerline4="platform=Illumina mix_type=physical"
    headerline5="aligner=BWA_MEM"
    prplotHeader="${headerline1}\n${headerline2}\n${headerline3}\n${headerline4}\n${headerline5};notchinese-notbest"
    mixtype=physical-mix
    tsample=HGM12IL27X63X #NA12878-${nsample}-mixture
    nsample=HG002IL30X #NA24385
    tfname=24385-12878-30-200_RX_001
    nfname=24385-200_AH5G7WCCXX_S4_L004_RX_001
else
    echo "The options must contain either use-physical-mixture, use-pa12-mixture, or use-pa21-mixture"
    exit -1
fi

FQDIR="${MIXRDIR}/${mixtype}/fq/"
BAMDIR="${MIXRDIR}/${mixtype}/bam/"
VCFDIR="${MIXRDIR}/${mixtype}/vcf/"
${mkdir777} "${FQDIR}"
${mkdir777} "${BAMDIR}"
${mkdir777} "${VCFDIR}"

if isfound "${PAT}" "use-physical-mix|use-physical-mixture"; then
    TRUTH_VCF_GZ="${FQDIR}/na12878-na24385-somatic-truth.vcf.gz"
    TRUTH_BED_ALL="${FQDIR}/na12878-na24385-somatic-truth-regions.bed"
else
    TRUTH_VCF_GZ="${FQDIR}/${silicoTname}-${silicoNname}-somatic-truth.vcf.gz"
    TRUTH_BED_ALL="${FQDIR}/${silicoTname}-${silicoNname}-somatic-truth-regions.bed"
    if isfound "${PAT}" "prep-truth" ; then
        bedtools intersect -a "${HGtumorBED}" -b "${HGnormalBED}" > "${TRUTH_BED_ALL}"
        #bcftools merge "${HGtumorVCF}" "${HGnormalVCF}" | bcftools view -f "$TRUTH_VIEW_EXPR" -Oz -o "${TRUTH_VCF_GZ}"
        bcftools isec ${GROUNDTRUTH_ISEC_FLAGS} --threads 8 "${HGtumorVCF}" "${HGnormalVCF}" -Oz -p "${FQDIR}/${silicoTname}-${silicoNname}-somatic-truth.isecdir"
        bcftools norm ${GROUNDTRUTH_NORM_FLAGS} --threads 8 -Oz -o "${TRUTH_VCF_GZ}" "${FQDIR}/${silicoTname}-${silicoNname}-somatic-truth.isecdir/0000.vcf.gz"
        bcftools index --threads 8 -ft "${TRUTH_VCF_GZ}"
    fi
fi

alltbam="${BAMDIR}/${tsample}.bam"
allnbam="${BAMDIR}/${nsample}.bam"

if isfound "${PAT}" "run-bwa" ; then
    ulimit -n 4096

    fq1="${FQDIR}/24385-12878-30-200_R1_001.fastq.gz"
    fq2="${FQDIR}/24385-12878-30-200_R2_001.fastq.gz"
    echo "Running bwa mem on the GIAB tumor-like sample"
    date && time -p "${BWA}" mem -R "@RG\tID:${tsample}\tSM:${tsample}\tLB:${tsample}\tPU:L001\tPL:ILLUMINA" \
            -t "${ncpus}" "${HS37D5}" "${fq1}" "${fq2}" | samtools view -bh1 -@ "${nbams}" | samtools sort -@ "${nbams}" -o "${alltbam}"
    time -p samtools index -@ "${ncpus}" "${alltbam}"
    fq1="${FQDIR}/24385-200_AH5G7WCCXX_S4_L004_R1_001.fastq.gz"
    fq2="${FQDIR}/24385-200_AH5G7WCCXX_S4_L004_R2_001.fastq.gz"
    echo "Running bwa mem on the GIAB normal-like sample"
    date && time -p "${BWA}" mem -R "@RG\tID:${nsample}\tSM:${nsample}\tLB:${nsample}\tPU:L001\tPL:ILLUMINA" \
            -t "${ncpus}" "${HS37D5}" "${fq1}" "${fq2}" | samtools view -bh1 -@ "${nbams}" | samtools sort -@ "${nbams}" -o "${allnbam}"
    time -p samtools index -@ "${ncpus}" "${allnbam}"
fi

TRUTH_BED="${TRUTH_BED_ALL}"
tbam="${alltbam}"
nbam="${allnbam}"
resultRootDir="${VCFDIR}/"
outname=${tsample}-${nsample}-nochrall

chmod uga+rwxs "${resultRootDir}" || true
#setfacl -d -m g::rwx "${resultRootDir}"
#setfacl -d -m o::rwx "${resultRootDir}"

for tname in $TNAMES; do
if isfound "${PAT}" use-tname-${tname} ; then
    #${mkdir777} -p "${BAMDIR}/subset/${tname}"
    ${mkdir777} -p "${VCFDIR}/subset/${tname}"
    
    TRUTH_BED="${VCFDIR}/subset/${tname}/somatic-truth-regions.bed"
    cat "${TRUTH_BED_ALL}" | grep -P "^${tname}\t" > "${TRUTH_BED}" || true
    tbam="${BAMDIR}/subset/${tfname}-${tname}.bam"
    nbam="${BAMDIR}/subset/${nfname}-${tname}.bam"
    if isfound "${PAT}" "prep-tname-${tname}" ; then
        samtools view -bh -@ "${ncpus}" "${alltbam}" ${tname} > "${tbam}"
        samtools index -@ "${ncpus}" "${tbam}"
        samtools view -bh -@ "${ncpus}" "${allnbam}" ${tname} > "${nbam}"
        samtools index -@ "${ncpus}" "${nbam}"
    fi
    resultRootDir="${VCFDIR}/subset/${tname}/" #${scenario}/"
    outname=${tsample}-${nsample}-nochr${tname}
    for headnum in 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536; do
        headxend="head-${headnum}-bedlines"
        if isfound "${PAT}" ${headxend} ; then
            
            prevTRUTH_BED="${TRUTH_BED}"
            TRUTH_BED="${TRUTH_BED/%.bed/}.${tname}.${headxend}.bed"
            head -n ${headnum} "${prevTRUTH_BED}" > "${TRUTH_BED}"
            
            tname_begpos=$(head -n1 "${TRUTH_BED}" | awk '{print $2}')
            tname_endpos=$(tail -n1 "${TRUTH_BED}" | awk '{print $3}')
            samtools_view_region="${tname}:${tname_begpos}-${tname_endpos}"
            tsbamdir="${BAMDIR}/subset/${tfname}-${tname}.bam-region.dir/"
            nsbamdir="${BAMDIR}/subset/${nfname}-${tname}.bam-region.dir/"
            
            ${mkdir777} -p "${tsbamdir}" "${nsbamdir}"
            prevtbam="${tbam}"
            prevnbam="${nbam}"
            tbam="${tsbamdir}/${tfname}-${tname}.${headxend}.bam"
            nbam="${nsbamdir}/${nfname}-${tname}.${headxend}.bam"
            if isfound "${PAT}" run-bedlines ; then
                samtools view "${prevtbam}" "${samtools_view_region}" -bh > "${tbam}"
                samtools index "${tbam}"
                samtools view "${prevnbam}" "${samtools_view_region}" -bh > "${nbam}"
                samtools index "${nbam}"
            fi
            resultRootDir="${VCFDIR}/subset/${tname}-${headxend}/"
        fi
    done
fi
done

source "${EVALROOT}/source-tissue.sh"

recaltbam="${tbam/%.bam/.rmdup_recal.bam}"
recalnbam="${nbam/%.bam/.rmdup_recal.bam}"
dindeltbam="${tbam/%.bam/.rmdup_recal_dindel.bam}"
dindelnbam="${nbam/%.bam/.rmdup_recal_dindel.bam}"

run_bamprep "${PAT}" "${HS37D5}" "${tbam}" "${EVALROOT}/datafiles/hg19/dbsnp_138.b37.vcf" "${ncpus}"
run_bamprep "${PAT}" "${HS37D5}" "${nbam}" "${EVALROOT}/datafiles/hg19/dbsnp_138.b37.vcf" "${ncpus}"

if isfound "${PAT}" enable-tn-pair ; then
    run_tnpair "${PAT}.enable-uvc1-all.enable-strelka2-all.enable-lolopicker-all.enable-octopus-all.enable-varscan2-all.enable-somaticsniper-all" \
        "${resultRootDir}/${scenario}/" "${HS37D5}" "${TRUTH_VCF_GZ}" "${TRUTH_BED}" "${tbam}" "${tsample}" "${nbam}" "${nsample}" "${ncpus}"
    run_tnpair "${PAT}.enable-gatk4mutect2-all.enable-lancet-all" \
        "${resultRootDir}/${scenario}/" "${HS37D5}" "${TRUTH_VCF_GZ}" "${TRUTH_BED}" "${recaltbam}" "${tsample}" "${recalnbam}" "${nsample}" "${ncpus}"
    run_tnpair "${PAT}.enable-lofreq-all" \
        "${resultRootDir}/${scenario}/" "${HS37D5}" "${TRUTH_VCF_GZ}" "${TRUTH_BED}" "${dindeltbam}" "${tsample}" "${dindelnbam}" "${nsample}" "${ncpus}"
elif isfound "${PAT}" enable-tumor-only; then
    run_tonly "${PAT}.enable-uvc1-all.enable-strelka2-all.enable-varscan2-all.enable-outlyzer-all" \
        "${resultRootDir}/${scenario}/" "${HS37D5}" "${TRUTH_VCF_GZ}" "${TRUTH_BED}" "${tbam}" "${tsample}" "${HGnormalVCF}" "${ncpus}"
    run_tonly "${PAT}.enable-gatk4mutect2-all" \
        "${resultRootDir}/${scenario}/" "${HS37D5}" "${TRUTH_VCF_GZ}" "${TRUTH_BED}" "${recaltbam}" "${tsample}" "${HGnormalVCF}" "${ncpus}"
    run_tonly "${PAT}.enable-lofreq-all" \
        "${resultRootDir}/${scenario}/" "${HS37D5}" "${TRUTH_VCF_GZ}" "${TRUTH_BED}" "${dindeltbam}" "${tsample}" "${HGnormalVCF}" "${ncpus}"
else
    echo "The options should contain either enable-tn-pair or enable-tumor-only."
    exit -2
fi

if isfound "${PAT}" "run-powerlaw-plot"; then
    resdir="${resultRootDir}/${scenario}/"
    python "${EVALROOT}/common/custom-qualplot.py" "${prplotHeader}" \
            "${resdir}/${tsample}_tonly.${PARAMSET}_norm_vcfeval-all.dir/fp.vcf.gz" \
            "${TRUTH_BED}" \
            "${resdir}/${tsample}_${nsample}_uvc-powerlaw-plot.pdf" > "${resdir}/${tsample}_${nsample}_uvc-powerlaw-plot.tsv.xls"
fi

if isfound "${PAT}" "run-custom-prplot"; then
    
    resdir="${resultRootDir}/${scenario}/"
    mkdir -p "${resdir}"

    if isfound "${PAT}" enable-tn-pair ; then
        eval_uvc="${resdir}/${tsample}_${nsample}.${PARAMSET}_norm_vcfeval-all.dir"
        eval_mutect2raw="${resdir}/${tsample}_${nsample}.gatk4mutect2_norm_vcfeval-all.dir"
        eval_mutect2filt="${resdir}/${tsample}_${nsample}.gatk4mutect2_norm_vcfeval-filt.dir"
        eval_strelka2="${resdir}/${tsample}_${nsample}.strelka2_norm_vcfeval-all.dir"
        eval_varscan2raw="${resdir}/${tsample}_${nsample}.varscan2_norm_vcfeval-all.dir"
        eval_varscan2filt="${resdir}/${tsample}_${nsample}.varscan2_norm_vcfeval-filt.dir"
        eval_somaticsniper="${resdir}/${tsample}_${nsample}.somaticsniper_norm_vcfeval-all.dir"
        eval_lofreq="${resdir}/${tsample}_${nsample}.lofreq_norm_vcfeval-all.dir"
        eval_lolopicker="${resdir}/${tsample}_${nsample}.lolopicker_norm_vcfeval-all.dir"
        eval_octopus="${resdir}/${tsample}_${nsample}.octopus_norm_vcfeval-all.dir"
        eval_lancet="${resdir}/${tsample}_${nsample}.lancet_norm_vcfeval-all.dir"

        rocfile="/snp_roc.tsv.gz"
        python "${EVALROOT}/common/custom-prplot.py" "${resdir}/${tsample}_${nsample}_all-methods-custom-prplot_snp.pdf"     "${prplotHeader}" 0.00,0.00 \
            "${eval_uvc}${rocfile}" \
            "${eval_mutect2raw}${rocfile}" \
            "${eval_mutect2filt}${rocfile}" \
            "${eval_strelka2}${rocfile}" \
            "${eval_varscan2raw}${rocfile}" \
            "${eval_varscan2filt}${rocfile}" \
            "${eval_somaticsniper}${rocfile}" \
            "${eval_lofreq}${rocfile}" \
            "${eval_lolopicker}${rocfile}" \
            "${eval_octopus}${rocfile}" \
            "${eval_lancet}${rocfile}" \
 
        rocfile="/non_snp_roc.tsv.gz"
        python "${EVALROOT}/common/custom-prplot.py" "${resdir}/${tsample}_${nsample}_all-methods-custom-prplot_non_snp.pdf" "${prplotHeader}" 0.00,0.00 \
            "${eval_uvc}${rocfile}" \
            "${eval_mutect2raw}${rocfile}" \
            "${eval_mutect2filt}${rocfile}"\
            "${eval_strelka2}${rocfile}" \
            "${eval_varscan2raw}${rocfile}" \
            "${eval_varscan2filt}${rocfile}" \
            "${eval_somaticsniper}${rocfile}" \
            "${eval_lofreq}${rocfile}" \
            "${eval_lolopicker}${rocfile}" \
            "${eval_octopus}${rocfile}" \
            "${eval_lancet}${rocfile}" \
 
    elif isfound "${PAT}" enable-tumor-only; then
        eval_uvc="${resdir}/${tsample}_tonly.${PARAMSET}_norm_vcfeval-all.dir"
        eval_mutect2raw="${resdir}/${tsample}_tonly.gatk4mutect2_norm_vcfeval-all.dir"
        eval_mutect2filt="${resdir}/${tsample}_tonly.gatk4mutect2_norm_vcfeval-filt.dir"
        eval_lofreq="${resdir}/${tsample}_tonly.lofreq_norm_vcfeval-all.dir"
        
        rocfile="/snp_roc.tsv.gz"
        python "${EVALROOT}/common/custom-prplot.py" "${resdir}/${tsample}_tonly_all-methods-custom-prplot_snp.pdf"     "${prplotHeader}" 0.00,0.00 \
            "${eval_uvc}${rocfile}" \
            "${eval_mutect2raw}${rocfile}" \
            "${eval_mutect2filt}${rocfile}" \
            "${eval_lofreq}${rocfile}" \
        
        rocfile="/non_snp_roc.tsv.gz"
        python "${EVALROOT}/common/custom-prplot.py" "${resdir}/${tsample}_tonly_all-methods-custom-prplot_non_snp.pdf" "${prplotHeader}" 0.00,0.00 \
            "${eval_uvc}${rocfile}" \
            "${eval_mutect2raw}${rocfile}" \
            "${eval_mutect2filt}${rocfile}"\
            "${eval_lofreq}${rocfile}" \
        
    fi
fi

}
exit $?

