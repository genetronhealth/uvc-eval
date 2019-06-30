#!/usr/bin/env bash

# https://dev.biologists.org/content/145/6/dev153049#T1 
# https://www.ebi.ac.uk/ena/browser/view/PRJEB38047  ERP121418
{
set -evx

patients="236 245 282 300 469 575 610 664 708 726 735 772 781 815 861 GI_21"

SCRIPTDIR=$(dirname $(which $0))

if [ -z "${EVALROOT}" ]; then
    source "${SCRIPTDIR}/../source-eval.sh"
else
    source "${EVALROOT}/source-eval.sh"
fi

fqdir="${EVALROOT}/ERP015684/fq/"
bamdir="${EVALROOT}/ERP015684/bam/"
vcfdir="${EVALROOT}/ERP015684/vcf/"

${mkdir777} -p "${fqdir}" "${bamdir}" "${vcfdir}"

for patient in ${patients}; do

if [ $(echo "${PAT}" | grep -cP "use-patient-${patient}-end|use-patient-all-end") -eq 0 ]; then
    continue
fi

tbam="${bamdir}/${patient}_TUMOR_R1R2_trimPaired.bam"
tsample=$(echo "T${patient}" | sed 's/GI_//g')
nbam="${bamdir}/${patient}_NORMAL_R1R2_trimPaired.bam"
nsample=$(echo "N${patient}" | sed 's/GI_//g')

resultRootDir="${vcfdir}/${tsample}"
${mkdir777} -p "${resultRootDir}"

if isfound "${PAT}" "run-bwa" ; then
    r1="${fqdir}/${patient}_TUMOR_R1_trimPaired.fastq.gz"
    r2="${fqdir}/${patient}_TUMOR_R2_trimPaired.fastq.gz"
    platform=ILLUMINA
    rgline="@RG\tID:${tsample}\tSM:${tsample}\tLB:${tsample}\tPL:${platform}\tPU:${tsample}"
    "${BWA}" mem -t ${ncpus} -R "${rgline}" "${HS37D5}" "${r1}" "${r2}" \
        | samtools view -@ "${nbams}" -bh - | samtools sort -@ "${nbams}" -o "${tbam}" -
    samtools index -@ ${nbams} "${tbam}"
    
    r1="${fqdir}/${patient}_NORMAL_R1_trimPaired.fastq.gz"
    r2="${fqdir}/${patient}_NORMAL_R2_trimPaired.fastq.gz"
    platform=ILLUMINA
    rgline="@RG\tID:${nsample}\tSM:${nsample}\tLB:${nsample}\tPL:${platform}\tPU:${nsample}"
    "${BWA}" mem -t ${ncpus} -R "${rgline}" "${HS37D5}" "${r1}" "${r2}" \
        | samtools view -@ "${nbams}" -bh - | samtools sort -@ "${nbams}" -o "${nbam}" -
    samtools index -@ ${nbams} "${nbam}"
fi

#recalbam="${tbam/%.bam/.rmdup_recal.bam}"
#dindeltbam="${tbam/%.bam/.rmdup_recal_dindel.bam}"
#run_bamprep "${PAT}" "${HS37D5}" "${tbam}" "${EVALROOT}/datafiles/hg19/dbsnp_138.b37.vcf" "${ncpus}"

if isfound "${PAT}" enable-tn-pair; then
    source "${EVALROOT}/source-tissue.sh"
    run_tnpair "${PAT}.enable-uvc1-all.enable-strelka2-all.enable-lolopicker-all.enable-varscan2-all.enable-somaticsniper-all" \
        "${resultRootDir}/" "${HS37D5}" "{TRUTH_VCF_GZ}" "{TRUTH_BED}" "${tbam}" "${tsample}" "${nbam}" "${nsample}" "${ncpus}"
else
    echo "enable-tn-pair is not set in options"
fi
done
}
exit $?

