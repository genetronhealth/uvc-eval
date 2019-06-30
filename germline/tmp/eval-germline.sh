#!/usr/bin/env bash
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

chrnum=1
hg=HG001
avgdep=30

if [ $(echo "${PAT}" | grep -P "use-1-nochr|samtools-view-1-end" | wc -l) -gt 0 ]; then
    chrnum=1
fi
if [ $(echo "${PAT}" | grep -P "use-22-nochr|samtools-view-22-end" | wc -l) -gt 0 ]; then
    chrnum=22
fi
if [ $(echo "${PAT}" | grep -P "use-all-nochr|use-999-nochr|samtools-view-all-end|samtools-view-999-end" | wc -l) -gt 0 ]; then
    chrnum=999
fi

if [ $(echo "${PAT}" | grep "use-hg001-sample" | wc -l) -gt 0 ]; then
    hg=HG001
fi
if [ $(echo "${PAT}" | grep "use-hg002-sample" | wc -l) -gt 0 ]; then
    hg=HG002
fi
if [ $(echo "${PAT}" | grep "use-hg005-sample" | wc -l) -gt 0 ]; then
    hg=HG005
fi

if [ $(echo "${PAT}" | grep "use-30x-avgdep" | wc -l) -gt 0 ]; then
    avgdep=30
fi
if [ $(echo "${PAT}" | grep "use-60x-avgdep" | wc -l) -gt 0 ]; then
    avgdep=60
fi
if [ $(echo "${PAT}" | grep "use-300x-avgdep" | wc -l) -gt 0 ]; then
    avgdep=300
fi

suf=nochr${chrnum}
AVGDEP="${avgdep}x"

bam="${EVALROOT}/germline/bam/${hg}.hs37d5.${AVGDEP}.${suf}.leftaligned.bam"
#bam="${EVALROOT}/germline/bam/${hg}.hs37d5.${AVGDEP}.${suf}.bam"
resrootdir="${EVALROOT}/germline/vcf/${PARAMSET}/"
dir="${EVALROOT}/germline/vcf/${PARAMSET}/${hg}.hs37d5.${AVGDEP}.${suf}.uvc.resdir"

if [ $(echo "${PAT}" | grep "use-abra2" | wc -l) -gt 0 ]; then
    bam="${EVALROOT}/germline/bam/${hg}.hs37d5.${AVGDEP}.${suf}.abra2.bam"
    dir="${EVALROOT}/germline/vcf/${hg}.hs37d5.${AVGDEP}.${suf}.uvc-abra2.resdir"
fi

resvcf="${dir}/${hg}-${suf}.uvc1.vcf.gz"
${mkdir777} -p "${dir}"

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

if [ ${chrnum} -lt 25 ]; then
    cat "${TRUTH_BED}" | grep -P "^${chrnum}\t" > "${truth_bed}"
else
    cp "${TRUTH_BED}" "${truth_bed}"
fi

function generate_disconcord_summary {
    germvcf=$1
    germout=$2
    bcftools merge --force-samples "${TRUTH_VCF}" "${germvcf}" -R ${truth_bed} | bcftools stats -s - > "${germvcf}.mstats"
    bcftools isec -Oz -c none "${germout}/fn.vcf.gz" "${germout}/fp.vcf.gz" -p "${germout}/f.isecdir/"
    bcftools view "${germout}/tp.vcf.gz" -i "GQ<10" -Oz -o "${germout}/f.isecdir/lowGQ.vcf.gz"
    bcftools index -ft "${germout}/f.isecdir/lowGQ.vcf.gz"
    bcftools isec -Oz -c none "${germout}/f.isecdir/0000.vcf.gz" "${germvcf}" -p "${germout}/f.isecdir/0000/"
    chmod -R 777 "${germout}"
}

if [ $(echo "${PAT}" | grep "run-bwa-BGI" | wc -l) -gt 0 ]; then
    for lane in L01 L02; do
        if [ "${hg}" == "HG001" ]; then
            fq1="${EVALROOT}/germline/fq/HG001/BGISEQ500_PCRfree_NA12878_CL100076243_${lane}_read_1.fq.gz"
            fq2="${EVALROOT}/germline/fq/HG001/BGISEQ500_PCRfree_NA12878_CL100076243_${lane}_read_2.fq.gz"
        fi
        if [ "${hg}" == "HG002" ]; then
            fq1="${EVALROOT}/germline/fq/HG002/BGISEQ500_PCRfree_NA24385_CL100076190_${lane}_read_1.fq.gz"
            fq2="${EVALROOT}/germline/fq/HG002/BGISEQ500_PCRfree_NA24385_CL100076190_${lane}_read_2.fq.gz"
        fi
        if [ "${hg}" == "HG005" ]; then
            fq1="${EVALROOT}/germline/fq/HG005/BGISEQ500_PCRfree_NA24631_CL100076244_${lane}_read_1.fq.gz"
            fq2="${EVALROOT}/germline/fq/HG005/BGISEQ500_PCRfree_NA24631_CL100076244_${lane}_read_2.fq.gz"
        fi
        sample="${hg}"
        library="${hg}BS30X"
        platunit="${hg}BS30X${lane}"
        "${BWA}" mem -t ${ncpus} "${HS37D5}" "${fq1}" "${fq2}" | samtools addreplacerg -r "@RG\tID:${platunit}\tSM:${sample}\tLB:${library}\tPU:${platunit}\tPL:BGI" - \
                | samtools sort -@ ${nbams} -o "${EVALROOT}/germline/bam/${platunit}.bam" -
        samtools index "${EVALROOT}/germline/bam/${platunit}.bam"
    done
    lib60x="${hg}BS60X"
    samtools merge -@ ${ncpus} "${EVALROOT}/germline/bam/${lib60x}.bam" "${EVALROOT}/germline/bam/${library}L01.bam" "${EVALROOT}/germline/bam/${library}L02.bam"
    samtools index -@ ${ncpus} "${EVALROOT}/germline/bam/${lib60x}.bam"
fi

if [ $(echo "${PAT}" | grep "run-samtools-view" | wc -l) -gt 0 ]; then
    if [ "${avgdep}" -ge 300 ]; then
        samtools view -bh "${EVALROOT}/tn/insilico-rawinput/${hg}.hs37d5.300x.bam" ${chrnum} | "${EVALROOT}/tools/freebayes-1.3.2/bin/bamleftalign" -f "${VCREF}" > "${bam}"
    else
        downfrac=$(python -c "print(float(${avgdep})/300)")
        #samtools view -s ${downfrac} -bh "${EVALROOT}/tn/insilico-rawinput/${hg}.hs37d5.300x.bam" ${chrnum} | "${EVALROOT}/tools/freebayes-1.3.2/bin/bamleftalign" -f "${VCREF}" > "${bam}"
        if [ ${chrnum} -lt 25 ]; then
            samtools view -@ ${ncpus} -s ${downfrac} -bh "${EVALROOT}/tn/insilico-step2input/${hg}.hs37d5.300x.rg.leftaligned.bam" ${chrnum} > "${bam}"
        else
            samtools view -@ ${ncpus} -s ${downfrac} -bh "${EVALROOT}/tn/insilico-step2input/${hg}.hs37d5.300x.rg.leftaligned.bam" > "${bam}"
        fi
    fi
    samtools index -@ ${ncpus} "${bam}"
fi

date
if [ $(echo "${PAT}" | grep "run-uvc1" | wc -l) -gt 0 ]; then
    ${mkdir777} -p "${dir}" # leftaligned
    date && "${UVC}" -q5 --somaticGT 0 --outvar-flag 63 -f "${VCREF}" -o "${resvcf}" -s ${hg} -t "${ncpus}" "${bam}" 1> "${resvcf/.vcf.gz/.stdout}" 2> "${resvcf/.vcf.gz/.stderr}"
    date && bcftools index -ft "${resvcf}"
fi

# bcftools view ${resvcf} -i "GT != '0/0' && ANY_VAR=1 && ((type=='snp' && TLODQ>=31) || (TYPE!='snp' && (GT != '0/.' && GT != './0' && GT != './.')))" \
date
if [ $(echo "${PAT}" | grep "eval-uvc-somagerm1" | wc -l) -gt 0 ]; then
    germvcf="${dir}/${hg}-${suf}-germ-uvc.nothomref-withtriallelic.vcf.gz" # && GT != '0/.' && GT != './0' && GT != './.' 
    bcftools view ${resvcf} -i "GT != '0/0' && ANY_VAR=1 && ((type=='snp' && TLODQ>=31) || (TYPE!='snp' && TLODQ>=31))" \
        | python "${EVALROOT}/scripts/bcfanno.py" \
        | bcftools norm -m +both - \
        | awk 'OFS="\t" {if ($5 ~ ",") { gsub("^(0|1|2|.)/(0|1|2|.):", "1/2:", $10);} print; }' \
        | bcftools view -Oz -o "${germvcf}"
    bcftools index -ft "${germvcf}"
    #rm -r "${dir}/uvc1.vcfeal.outdir" || true 
    #"${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval -b "${TRUTH_VCF}" -c "${germvcf}" -e ${truth_bed} \
    #    -t "${EVALROOT}/datafiles/Homo_sapiens_assembly19.sdf" -o "${dir}/uvc1.vcfeal.outdir" ${evalflags} || true
    
    FTvcf="${dir}/${hg}-${suf}-germ-uvc.nothomref-withtriallelic-FT.vcf.gz"
    bcftools view "${germvcf}" -i "FORMAT/FT == 'PASS' || FORMAT/FT == 'GPBLR'" -Oz -o "${FTvcf}"
    bcftools index -ft "${FTvcf}"
    
    #rm -r  "${dir}/uvcFT.vcfeal.TLODQ-resvcf.outdir" "${dir}/uvcFT.vcfeal.TLODQ.outdir" "${dir}/uvcFT.vcfeal.squash-ploidy.outdir" "${dir}/uvcFT.vcfeal.outdir" || true
    
    #"${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval -b "${TRUTH_VCF}" -c "${resvcf}" -e ${truth_bed} \
    #    -t "${EVALROOT}/datafiles/Homo_sapiens_assembly19.sdf" -o "${dir}/uvcFT.vcfeal.TLODQ-resvcf.outdir" ${evalflags} -f INFO.TLODQ --squash-ploidy || true
    #"${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval -b "${TRUTH_VCF}" -c "${FTvcf}" -e ${truth_bed} \
    #    -t "${EVALROOT}/datafiles/Homo_sapiens_assembly19.sdf" -o "${dir}/uvcFT.vcfeal.TLODQ.outdir" ${evalflags} -f INFO.TLODQ --squash-ploidy || true
    #"${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval -b "${TRUTH_VCF}" -c "${FTvcf}" -e ${truth_bed} \
    #    -t "${EVALROOT}/datafiles/Homo_sapiens_assembly19.sdf" -o "${dir}/uvcFT.vcfeal.squash-ploidy.outdir" ${evalflags} --squash-ploidy || true
    
    #"${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval -b "${TRUTH_VCF}" -c "${FTvcf}" -e ${truth_bed} \
    #    -t "${EVALROOT}/datafiles/Homo_sapiens_assembly19.sdf" -o "${dir}/uvcFT.vcfeal.GLa-squash-ploidy.outdir" ${evalflags} -f GLA --squash-ploidy || true
    
    #rm -r "${dir}/uvc1.vcfeal.outdir" || true
    #"${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval -b "${TRUTH_VCF}" -c "${germvcf}" -e ${truth_bed} \
    #    -t "${EVALROOT}/datafiles/Homo_sapiens_assembly19.sdf" -o "${dir}/uvc1.vcfeal.outdir" ${evalflags} || true
    
    rm -r "${dir}/uvcFT.vcfeal.outdir" || true
    "${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval -b "${TRUTH_VCF}" -c "${FTvcf}" -e ${truth_bed} \
        -t "${EVALROOT}/datafiles/Homo_sapiens_assembly19.sdf" -o "${dir}/uvcFT.vcfeal.outdir" ${evalflags} || true
    
    rm -r "${dir}/uvcFT.vcfeal.NonHomrefQ.outdir" || true
    "${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval -b "${TRUTH_VCF}" -c "${FTvcf}" -e ${truth_bed} \
        -t "${EVALROOT}/datafiles/Homo_sapiens_assembly19.sdf" -o "${dir}/uvcFT.vcfeal.NonHomrefQ.outdir" ${evalflags} -f NonHomrefQ --squash-ploidy || true
    
    rm -r "${dir}/uvcFT.vcfeal.INFO-TLODQ.outdir" || true
    "${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval -b "${TRUTH_VCF}" -c "${FTvcf}" -e ${truth_bed} \
        -t "${EVALROOT}/datafiles/Homo_sapiens_assembly19.sdf" -o "${dir}/uvcFT.vcfeal.INFO-TLODQ.outdir" ${evalflags} -f INFO.TLODQ --squash-ploidy || true
fi

if [ $(echo "${PAT}" | grep -P "eval-uvc1|eval-uvc-germline1" | wc -l) -gt 0 ]; then
    germvcf="${dir}/${hg}-${suf}.uvc1-germline.vcf.gz" 
    # bcftools view "${resvcf}" -i 'GERMLINE = 1 && GT != "ref"' -Oz -o "${germvcf}"
    bcftools view "${resvcf}" -i 'GERMLINE = 1' -Oz -o "${germvcf}"
    bcftools index -ft "${germvcf}"
    
    rm -r "${dir}/uvcgerm.vcfeal.outdir" || true
    "${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval -b "${TRUTH_VCF}" -c "${germvcf}" -e ${truth_bed} \
        -t "${EVALROOT}/datafiles/Homo_sapiens_assembly19.sdf" -o "${dir}/uvcgerm.vcfeal.outdir" ${evalflags} || true
    generate_disconcord_summary "${germvcf}" "${dir}/uvcgerm.vcfeal.outdir"
    #bcftools merge --force-samples "${TRUTH_VCF}" "${germvcf}" -R ${truth_bed} | bcftools stats -s - > "${germvcf}.mstats"
    #bcftools isec -Oz -c none "${dir}/uvcgerm.vcfeal.outdir/fn.vcf.gz" "${dir}/uvcgerm.vcfeal.outdir/fp.vcf.gz" -p "${dir}/uvcgerm.vcfeal.outdir/f.isecdir/"
    #bcftools view "${dir}/uvcgerm.vcfeal.outdir/tp.vcf.gz" -i "GQ<10" -Oz -o "${dir}/uvcgerm.vcfeal.outdir/f.isecdir/lowGQ.vcf.gz"
    #bcftools index -ft "${dir}/uvcgerm.vcfeal.outdir/f.isecdir/lowGQ.vcf.gz"
    #bcftools isec -Oz -c none "${dir}/uvcgerm.vcfeal.outdir/f.isecdir/0000.vcf.gz" "${germvcf}" -p "${dir}/uvcgerm.vcfeal.outdir/f.isecdir/0000/"
    chmod -R 777 "${dir}/uvcgerm.vcfeal.outdir/"
fi

dir="${EVALROOT}/germline/vcf/gatk4hc/${hg}.hs37d5.${AVGDEP}.${suf}.gatk4hc.resdir"
resvcftxt="${dir}/${hg}-${suf}.gatk4hc.vcf"
resvcf="${dir}/${hg}-${suf}.gatk4hc.vcf.gz"
date
${mkdir777} -p "${dir}"
bam2="${dir}/${hg}-${suf}_rmdup.bam"
bam2a="${dir}/${hg}-${suf}_rmdup_rg.bam"
bam3="${dir}/${hg}-${suf}_rmdup_recal.bam"
rmdup_txt="${dir}/${hg}-${suf}_rmdup.txt"
recal_txt="${dir}/${hg}-${suf}_rmdup_recal.txt"

if [ $(echo "${PAT}" | grep -P "run-gatk4all|run-gatk4prep|run-gatk4rmdup" | wc -l) -gt 0 ]; then
    date && ${gatk4lowmem} MarkDuplicates --ASSUME_SORT_ORDER coordinate --REMOVE_DUPLICATES true -I "${bam}" -M "${rmdup_txt}" -O "${bam2}"
fi
if [ $(echo "${PAT}" | grep -P "run-gatk4all|run-gatk4prep|run-gatk4bqsr" | wc -l) -gt 0 ]; then
    echo NOT USED'''
    date &&  ${gatk4} AddOrReplaceReadGroups \
      -I "${bam2}" \
      -O "${bam2a}" \
      -RGID 4 \
      -RGLB lib1 \
      -RGPL illumina \
      -RGPU unit1 \
      -RGSM 20
    date && ${gatk4} BaseRecalibrator -I "${bam2a}" --known-sites "${EVALROOT}/datafiles/hg19/dbsnp_138.b37.vcf" -O "${recal_txt}" -R "${VCREF}"
    date && ${gatk4} ApplyBQSR -bqsr "${recal_txt}" -I "${bam2a}" -O "${bam3}"
    '''
    date && ${gatk4} BaseRecalibrator -I "${bam2}" --known-sites "${EVALROOT}/datafiles/hg19/dbsnp_138.b37.vcf" -O "${recal_txt}" -R "${VCREF}"
    date && ${gatk4} ApplyBQSR -bqsr "${recal_txt}" -I "${bam2}" -O "${bam3}"
fi
if [ $(echo "${PAT}" | grep -P "run-gatk4all|run-gatk4hc|run-haplotypercaller" | wc -l) -gt 0 ]; then
    date && ${gatk4} HaplotypeCaller -I "${bam3}" -O "${resvcftxt}" -R "${VCREF}"
    date && bcftools view -Oz -o "${resvcf}" "${resvcftxt}"
    bcftools index -ft "${resvcf}"
fi

date
if [ $(echo "${PAT}" | grep -P "eval-gatk4hc|eval-haplotypecaller" | wc -l) -gt 0 ]; then
    rm -r "${dir}/gatk4hc.vcfeal.outdir" || true
    "${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval -b "${TRUTH_VCF}" -c "${resvcf}" -e ${truth_bed} \
        -t "${EVALROOT}/datafiles/Homo_sapiens_assembly19.sdf" -o "${dir}/gatk4hc.vcfeal.outdir" || true
    generate_disconcord_summary "${resvcf}" "${dir}/gatk4hc.vcfeal.outdir"
fi

dir="${EVALROOT}/germline/${hg}.hs37d5.${AVGDEP}.${suf}.freebayes.resdir"
resvcf="${dir}/${hg}-${suf}-germ-freebayes.vcf.gz"
date
if [ $(echo "${PAT}" | grep "run-freebayes" | wc -l) -gt 0 ]; then
    ${mkdir777} -p "${dir}"
    "${EVALROOT}/tools/freebayes-1.3.2/bin/freebayes" -f "${VCREF}" "${bam}" | bcftools view -Oz -o "${resvcf}"
    bcftools index -ft "${resvcf}"
fi

date
if [ $(echo "${PAT}" | grep "eval-freebayes" | wc -l) -gt 0 ]; then
    rm -r "${dir}/freebayes.vcfeal.outdir" || true
    "${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval -b "${TRUTH_VCF}" -c "${resvcf}" -e ${truth_bed} \
        -t "${EVALROOT}/datafiles/Homo_sapiens_assembly19.sdf" -o "${dir}/freebayes.vcfeal.outdir" || true
    generate_disconcord_summary "${resvcf}" "${dir}/freebayes.vcfeal.outdir"
fi

dir="${EVALROOT}/germline/vcf/strelka2/${hg}.hs37d5.${AVGDEP}.${suf}.strelka2.resdir"
rundir="${dir}/${hg}-${suf}-germ-strelka2.rundir"
date
if [ $(echo "${PAT}" | grep "run-strelka2" | wc -l) -gt 0 ]; then
    ${mkdir777} -p "${dir}"
    rm -r "${rundir}" || true
    "${EVALROOT}/tools/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py" \
        --referenceFasta="${VCREF}" --runDir="${rundir}" --bam="${bam}"
    "${rundir}/runWorkflow.py" -j "${ncpus}" -m local 1> "${rundir}/stdout-disableEVS.log" 2> "${rundir}/stderr-disableEVS.log"
fi

resvcf="${rundir}/results/variants/variants.vcf.gz"

date
if [ $(echo "${PAT}" | grep "eval-strelka2" | wc -l) -gt 0 ]; then
    rm -r "${dir}/strelka2.vcfeal.outdir" || true
    "${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval -b "${TRUTH_VCF}" -c "${resvcf}" -e ${truth_bed} \
        -t "${EVALROOT}/datafiles/Homo_sapiens_assembly19.sdf" -o "${dir}/strelka2.vcfeal.outdir" || true
    generate_disconcord_summary "${resvcf}" "${dir}/strelka2.vcfeal.outdir"
    #bcftools merge --force-samples "${TRUTH_VCF}" "${resvcf}" -R ${truth_bed} | bcftools stats -s - > "${resvcf}.mstats"
    #bcftools isec -Oz -c none "${dir}/strelka2.vcfeal.outdir/fn.vcf.gz" "${dir}/strelka2.vcfeal.outdir/fp.vcf.gz" -p "${dir}/strelka2.vcfeal.outdir/f.isecdir/"
    #bcftools view "${dir}/strelka2.vcfeal.outdir/tp.vcf.gz" -i "GQ<10" -Oz -o "${dir}/strelka2.vcfeal.outdir/f.isecdir/lowGQ.vcf.gz"
    #bcftools index -ft "${dir}/strelka2.vcfeal.outdir/f.isecdir/lowGQ.vcf.gz"
fi

java8="${EVALROOT}/tools/jre/bin/java"

dir="${EVALROOT}/germline/${hg}.hs37d5.${AVGDEP}.${suf}.abra2.dir"
resvcf="${EVALROOT}/germline/HG001.hs37d5.${AVGDEP}.nochr22.abra2.resdir/HG001.hs37d5.${AVGDEP}.nochr22.cadabra2.vcf.gz"
date
if [ $(echo "${PAT}" | grep "run-cadabra2" | wc -l) -gt 0 ]; then
    ${mkdir777} -p "${dir}"
    "${java8}" -cp "${EVALROOT}/tools/abra2-2.22.jar" abra.cadabra.Cadabra --threads=${ncpus} --tumor "${EVALROOT}/germline/HG001.hs37d5.${AVGDEP}.nochr22.abra2.bam" --ref "${EVALROOT}/datafiles/Homo_sapiens_assembly19.fasta" | bcftools view -Oz -o "${resvcf}" -
    bcftools index -ft "${resvcf}"
fi

date
if [ $(echo "${PAT}" | grep "eval-cadabra2" | wc -l) -gt 0 ]; then
    rm -r "${dir}/cadabra2.vcfeal.outdir" || true
    "${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval -b "${TRUTH_VCF}" -c "${resvcf}" -e ${truth_bed} \
        -t "${EVALROOT}/datafiles/Homo_sapiens_assembly19.sdf" -o "${dir}/cadabra2.vcfeal.outdir" || true
fi

