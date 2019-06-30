#!/usr/bin/env bash
# need java8, fastq-dump (with aspera for faster download), bwa, samtools, bcftools, and vcfeval (from RTG)
# need GNU parallel but can be substituted

SCRIPTDIR=$(dirname $(which $0))
source "${SCRIPTDIR}/../source-eval.sh"

MER_VAR_OPT_TRUTH="-m+any" #all
MER_VAR_OPT="-m+any" #"-m- all" #all

function nop() {
    echo no-operation-for "${@}"
}

function sortvcf() {
    vcf="${1}"
    cat <(bcftools view "${vcf}" | grep -B10000 -m1 ^#CHROM) <(bcftools view -H "${vcf}" | sort -k1V,1 -k2n,2 -k4V,4 -k5V,5) #| bcftools view -Oz -o "${2}"
}

SRADIR="${EVALROOT}/mageri/data/SRP064450/sra/" 
mkdir -p "${SRADIR}"

FQDIR="${EVALROOT}/mageri/data/SRP064450/fq/" 
mkdir -p "${FQDIR}"
if [ $(echo "${PAT}" | grep "run-fastq-dump" | wc -l ) -gt 0 ]; then
    for sra in $(ls "${SRADIR}" | grep ".sra$" | grep -P "${PAT}") ; do
        date && time -p ${fastqdump} --gzip --split-files -O "${FQDIR}" "${SRADIR}/${sra}"
    done
fi

if [ $(echo "${PAT}" | grep "run-extract-barcode" | wc -l ) -gt 0 ]; then
    for fq1 in $(ls "${FQDIR}" | grep "_1.fastq.gz$" | grep -P "${PAT}") ; do
        fq2=${fq1/%_1.fastq.gz/_2.fastq.gz}
        out1=${fq1/%_1.fastq.gz/_1_barcoded.fastq.gz}
        out2=${fq2/%_2.fastq.gz/_2_barcoded.fastq.gz}
        echo time -p "${ROOTDIR}/scripts/extract-barcodes.py" "${FQDIR}/${fq1}" "${FQDIR}/${fq2}" "${FQDIR}/${out1}" "${FQDIR}/${out2}" 0 11
    done | parallel -j $ncpus
fi

BAMDIR="${EVALROOT}/mageri/data/SRP064450/bam/"
mkdir -p -m uga+rwx "${BAMDIR}"
if [ $(echo "${PAT}" | grep "run-bwa" | wc -l ) -gt 0 ]; then
    for fq1 in $(ls "${FQDIR}" | grep "_1_barcoded.fastq.gz$" | grep -P "${PAT}") ; do
        fq2=${fq1/%_1_barcoded.fastq.gz/_2_barcoded.fastq.gz}
        bam=${fq1/%_1_barcoded.fastq.gz/_12.bam}
        if [ $(echo "${PAT}" | grep run-bwa | wc -l) -gt 0 ]; then
            rm "${BAMDIR}/${bam}"* || true
            date
            time -p "${BWA}" mem -t $ncpus "${HG19}" "${FQDIR}/${fq1}" "${FQDIR}/${fq2}" | samtools view -bh1 | samtools sort -@8 -o "${BAMDIR}/${bam}"
            date
            time -p samtools index -@ $ncpus "${BAMDIR}/${bam}"
        fi
    done
fi

TRUTHVCF="${EVALROOT}/datafiles/hd734.hg19.vcf"
TRUTHVCFGZ="${EVALROOT}/datafiles/hd734.hg19.vcf.gz"
TRUTH_MA_VCFGZ="${EVALROOT}/datafiles/hd734.hg19.ma-p-all.vcf.gz"
MAGERI="${EVALROOT}/mageri/mageri.jar"
PRESET1="${EVALROOT}/mageri/mageri-paper/processing/hd734_and_donors/preset.xml"
PRESET2="${EVALROOT}/mageri/mageri-paper/processing/hd734_and_donors/preset_1_1_1.xml"
META="${EVALROOT}/mageri/mageri-paper/processing/hd734_and_donors/meta/"
BASELINEDIR="${EVALROOT}/mageri/data/SRP064450/baseline/"

if [ $(echo "${PAT}" | grep run-dataprep | wc -l) -gt 0 ]; then
    bcftools view -Oz -o "${TRUTHVCFGZ}" "${TRUTHVCF}"
    bcftools index -ft "${TRUTHVCFGZ}"
    sortvcf "${TRUTHVCF}" | bcftools norm - "${MER_VAR_OPT_TRUTH}" -Oz -o "${TRUTH_MA_VCFGZ}" 
    bcftools index -ft "${TRUTH_MA_VCFGZ}"
    cat "${PRESET1}" | sed 's/1.1.1-SNAPSHOT/1.1.1/g' > "${PRESET2}"
fi
mkdir -p -m uga+rwx "${BASELINEDIR}"
VCFDIR="${EVALROOT}/mageri/data/SRP064450/vcf/${Tparamset}/"
mkdir -p -m uga+rwx "${VCFDIR}"

for fq1 in $(ls "${FQDIR}" | grep "_1_barcoded.fastq.gz$" | grep -P "${PAT}") ; do
    fq2=${fq1/%_1_barcoded.fastq.gz/_2_barcoded.fastq.gz}
    resultdir=${fq1/%_1_barcoded.fastq.gz/.resultdir}
    
    if [ $(echo "${PAT}" | grep run-mageri | wc -l) -gt 0 ]; then
        date && time -p ${java8} -Xmx32G -jar ${MAGERI} \
            --import-preset ${PRESET2} \
            --references    ${META}/refs.fa \
            --contigs       ${META}/contigs.txt \
            -M2             ${META}/primers.txt \
            --bed           ${META}/refs.bed \
            -R1             "${FQDIR}/${fq1}" \
            -R2             "${FQDIR}/${fq2}" \
            --threads       $ncpus \
            -O              "${BASELINEDIR}/${resultdir}" ;
    fi
    vcf="${BASELINEDIR}/${resultdir}/my_project.my_sample.vcf"
    vcfgz="${BASELINEDIR}/${resultdir}/my_project.my_sample.vcf.gz"
    bam="${fq1/%_1_barcoded.fastq.gz/.bam}"
    if [ $(echo "${PAT}" | grep "eval-mageri" | wc -l) -gt 0 ]; then
        sortvcf "${vcf}" | sed "s/my_project.my_sample/${bam}/g" | bcftools view -Oz -o "${vcfgz}"
        date && time -p bcftools index -ft "${vcfgz}"
        isecdir="${vcf/%.vcf/_baseline2call_isec.dir}"
        date && time -p bcftools isec -c none -Oz -p "${isecdir}" "${TRUTHVCFGZ}" "${vcfgz}"
        date && time -p bcftools index "${isecdir}/0003.vcf.gz"
    fi
done

for bam in $(ls ${BAMDIR} | grep ".bam$" | grep -P "${PAT}"); do
    vcf=${bam/%.bam/_uvc1.vcf.gz}
    if [ $(echo "${PAT}" | grep "run-uvc1" | wc -l) -gt 0 -o $(echo "${PAT}" | grep "run-bwa" | wc -l) -gt 0 ]; then
        echo time -p "${UVC}" -f "${HG19}" -o "${VCFDIR}/${vcf}" -t $ncpus -s "${bam}" "${BAMDIR}/${bam}" ${uvcparams} "2> ${VCFDIR}/${vcf/%.vcf.gz/.stderr.log}" \
                "&&" time -p bcftools index -ft "${VCFDIR}/${vcf}"
    fi
done | parallel -j ${nbams} || true ; 

for bam in $(ls ${BAMDIR} | grep ".bam$" | grep -P "${PAT}"); do
    vcf=${bam/%.bam/_uvc1.vcf.gz}
    isecdir="${vcf/%_uvc1.vcf.gz/_uvc1_baseline2call_isec.dir}"
    if [ $(echo "${PAT}" | grep "eval-uvc1" | wc -l) -gt 0 ]; then
        date && time -p bcftools isec -c none -Oz -p "${VCFDIR}/${isecdir}" "${TRUTHVCFGZ}" "${VCFDIR}/${vcf}"
        date && time -p bcftools annotate --collapse none -a "${TRUTHVCFGZ}" -c ID -Oz -o "${VCFDIR}/${isecdir}/0003_anno.vcf.gz" "${VCFDIR}/${isecdir}/0003.vcf.gz"
        #date && time -p bcftools view -i "FORMAT/DP>=1000" -Oz -o "${VCFDIR}/${vcf/%_uvc1.vcf.gz/_uvc1_hd734_atleast1000dp.vcf.gz}" "${VCFDIR}/${isecdir}/0003_anno.vcf.gz"
        date && time -p bcftools view -i "FORMAT/DP>=1000" "${VCFDIR}/${isecdir}/0003_anno.vcf.gz" \
            | bcftools annotate -x "FORMAT/VTD,FORMAT/gapSa" -Oz -o "${VCFDIR}/${vcf/%_uvc1.vcf.gz/_uvc1_hd734_atleast1000dp.vcf.gz}"         
        date && time -p bcftools index -ft "${VCFDIR}/${vcf/%_uvc1.vcf.gz/_uvc1_hd734_atleast1000dp.vcf.gz}"
    fi
done

function bcfmerge() {
    if [ $# -gt 2 -a "$2" != "-" -a "$3" != "-" ]; then
        date && time -p bcftools merge -m none --force-samples -Oz -o $@
    elif [ $# -eq 2 -a "$2" != "-" ]; then
        date && time -p cp "${2}" "${1}"
    else
        echo "${1} will be an emtpy gzipped vcf file."
        printf "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n" | bcftools view -Oz -o "${1}" -
    fi
    bcftools index -ft "${1}"
}

allSRAregex=="SRR2556933|SRR2556934|SRR2556935|SRR2556936|SRR2556937|SRR2556938|SRR2556939|SRR2556940|SRR2556941|SRR2556942|SRR2556943|SRR2556944|SRR2556945|SRR2556946|SRR2556947|SRR2556948"
healthySRAregex="SRR2556933|SRR2556934|SRR2556935|SRR2556936|SRR2556937|SRR2556938|SRR2556939|SRR2556940"
hd734SRAregexs=("SRR2556941|SRR2556942|SRR2556943|SRR2556944" "SRR2556945|SRR2556946|SRR2556947|SRR2556948")
hd734SRAregexAll="${hd734SRAregexs[0]}|${hd734SRAregexs[1]}"

if [ $(echo "${PAT}" | grep "eval-uvc1" | wc -l) -gt 0 ]; then
    
    bcfmerge "${VCFDIR}/SRP064450_uvc1_healthy_8.vcf.gz" $(ls "${VCFDIR}/"*"_uvc1_hd734_atleast1000dp.vcf.gz" | grep -P "${healthySRAregex}")
    bcftools view "${VCFDIR}/SRP064450_uvc1_healthy_8.vcf.gz" | sed 's;\t[0-4.][|/][0-4.];\t0/0;g' | sed 's/chr//g' | sed 's/my_project.my_sample/my_project.healthyCtrl/g' \
        | bcftools view -i "${UMI_FILTER_EXPR_2}" -Oz -o  "${VCFDIR}/SRP064450_uvc1_healthy_8_healthychr.vcf.gz" -
    bcftools index -ft "${VCFDIR}/SRP064450_uvc1_healthy_8_healthychr.vcf.gz"

    for i in 0 1; do
        bcfmerge "${VCFDIR}/SRP064450_uvc1_hd734_8.vcf.gz" $(ls "${VCFDIR}/"*"_uvc1_hd734_atleast1000dp.vcf.gz" | grep -P "${hd734SRAregexs[$i]}")
        bcftools merge --force-samples -m none -Oz -o "${VCFDIR}/SRP064450_uvc1_hd734_healthychr_8x2.vcf.gz" "${VCFDIR}/SRP064450_uvc1_hd734_8.vcf.gz" "${VCFDIR}/SRP064450_uvc1_healthy_8_healthychr.vcf.gz" 
        firstSample=$(bcftools view -h "${VCFDIR}/SRP064450_uvc1_hd734_healthychr_8x2.vcf.gz" | tail -n1 | awk '{print $(10)}')
        sortvcf "${VCFDIR}/SRP064450_uvc1_hd734_healthychr_8x2.vcf.gz" | sed 's;\t[0-4.][|/][0-4.];\t0/1;g' | bcftools norm - ${MER_VAR_OPT} -Ou -o - | bcftools view -s ${firstSample} -i "${UMI_FILTER_EXPR_2}" -Oz -o "${VCFDIR}/SRP064450_uvc1_hd734_healthychr_8x2_oneSample${i}.vcf.gz" -
        bcftools index -ft "${VCFDIR}/SRP064450_uvc1_hd734_healthychr_8x2_oneSample${i}.vcf.gz"
    done

    python "${EVALROOT}/scripts/evaluate-hd734-vs-healthy.py" <(ls "${VCFDIR}/"*"_uvc1_hd734_atleast1000dp.vcf.gz" | grep -P "${healthySRAregex}") <(ls "${VCFDIR}/"*"_uvc1_hd734_atleast1000dp.vcf.gz" | grep -P "${hd734SRAregexAll}") "all,pos-all" "DP>=1000" > "${VCFDIR}/evaluation-results-uvc1.tsv"
    python "${EVALROOT}/scripts/evaluate-hd734-vs-healthy.py" <(ls "${VCFDIR}/"*"_uvc1_hd734_atleast1000dp.vcf.gz" | grep -P "${healthySRAregex}") <(ls "${VCFDIR}/"*"_uvc1_hd734_atleast1000dp.vcf.gz" | grep -P "${hd734SRAregexAll}") "all,pos-all,FIELD=TLODQ" "DP>=1000" > "${VCFDIR}/evaluation-results-uvc1.TLODQ.tsv"
    echo Results are at "${VCFDIR}/evaluation-results-uvc1.tsv"
fi


if [ $(echo "${PAT}" | grep "eval-mageri" | wc -l) -gt 0 ]; then

    bcfmerge "${BASELINEDIR}/SRP064450_mageri_healthy_8.vcf.gz" $(ls "${BASELINEDIR}/"*"/my_project.my_sample_baseline2call_isec.dir/0003.vcf.gz" | grep -P "${healthySRAregex}")
    bcftools view "${BASELINEDIR}/SRP064450_mageri_healthy_8.vcf.gz" | sed 's;\t[0-4.][|/][0-4.];\t0/0;g' | sed 's/chr//g' | sed 's/my_project.my_sample/my_project.healthyCtrl/g' | bcftools view -Oz -o  "${BASELINEDIR}/SRP064450_mageri_healthy_8_healthychr.vcf.gz" -
    bcftools index -ft "${BASELINEDIR}/SRP064450_mageri_healthy_8_healthychr.vcf.gz"

    for i in 0 1; do
        bcfmerge "${BASELINEDIR}/SRP064450_mageri_hd734_8.vcf.gz" $(ls "${BASELINEDIR}/"*"/my_project.my_sample_baseline2call_isec.dir/0003.vcf.gz" | grep -P "${hd734SRAregexs[$i]}")
        bcftools merge -m none --force-samples -Oz -o "${BASELINEDIR}/SRP064450_mageri_hd734_healthychr_8x2.vcf.gz" "${BASELINEDIR}/SRP064450_mageri_hd734_8.vcf.gz" "${BASELINEDIR}/SRP064450_mageri_healthy_8_healthychr.vcf.gz" 
        firstSample=$(bcftools view -h "${BASELINEDIR}/SRP064450_mageri_hd734_healthychr_8x2.vcf.gz" | tail -n1 | awk '{print $(10)}')
        sortvcf "${BASELINEDIR}/SRP064450_mageri_hd734_healthychr_8x2.vcf.gz" | sed 's;\t[0-4.][|/][0-4.];\t0/1;g' | bcftools norm - ${MER_VAR_OPT} -Ou -o - | bcftools view -s ${firstSample} -Oz -o "${BASELINEDIR}/SRP064450_mageri_hd734_healthychr_8x2_oneSample${i}.vcf.gz" - 
        bcftools index -ft "${BASELINEDIR}/SRP064450_mageri_hd734_healthychr_8x2_oneSample${i}.vcf.gz"
    done

    python "${EVALROOT}/scripts/evaluate-hd734-vs-healthy.py" <(ls "${BASELINEDIR}/"*"/my_project.my_sample_baseline2call_isec.dir/0003.vcf.gz" | grep -P "${healthySRAregex}") <(ls "${BASELINEDIR}/"*"/my_project.my_sample_baseline2call_isec.dir/0003.vcf.gz" | grep -P "${hd734SRAregexAll}") "all,pos-all" "DP>=1000" > "${VCFDIR}/evaluation-results-mageri.tsv"

fi

if false; then
    FLAGS="--all-records" #" --squash-ploidy --ref-overlap --decompose" #" --ref-overlap --decompose "
    for i in 0 1; do
        rm -r "${VCFDIR}/SRP064450_uvc1_hd734_healthychr_8x2rep${i}.vcfeval.outdir" || true
        date && time -p "${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval --threads=16 \
            -f QUAL $FLAGS \
            -b "${TRUTH_MA_VCFGZ}" \
            -c "${VCFDIR}/SRP064450_uvc1_hd734_healthychr_8x2_oneSample${i}.vcf.gz" \
            -t "${EVALROOT}/datafiles/b37_plus_hg19_UCSC.sdf" \
            -o "${VCFDIR}/SRP064450_uvc1_hd734_healthychr_8x2rep${i}.vcfeval.outdir"

        rm -r "${VCFDIR}/SRP064450_mageri_hd734_healthychr_8x2rep${i}.vcfeval.outdir" || true
        date && time -p "${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" vcfeval --threads=16 \
            -f QUAL $FLAGS \
            -b "${TRUTH_MA_VCFGZ}" \
            -c "${BASELINEDIR}/SRP064450_mageri_hd734_healthychr_8x2_oneSample${i}.vcf.gz" \
            -t "${EVALROOT}/datafiles/b37_plus_hg19_UCSC.sdf" \
            -o "${VCFDIR}/SRP064450_mageri_hd734_healthychr_8x2rep${i}.vcfeval.outdir"

        date && time -p "${java8}" -jar "${EVALROOT}/tools/rtg-tools-3.10.1/RTG.jar" rocplot \
            "${VCFDIR}/SRP064450_"*"_hd734_healthychr_8x2rep${i}.vcfeval.outdir/"*"snp_roc.tsv.gz" --scores \
            --svg "${VCFDIR}/SRP064450_mageri_hd734_healthychr_8x2rep${i}.vcfeval.outdir/SRP064450.all_methods_all_vars_rocplot.svg"
    done
fi


