
strelka2="${EVALROOT}/tools/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py"
freebayes="${EVALROOT}/tools/freebayes-1.3.2/bin/freebayes"

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

function setup_germline {
    vcname="${1}"
    evaladir="${calldir}/${gsample}_germline.${vcname}_norm_vcfeval-all.dir"
    evalfdir="${calldir}/${gsample}_germline.${vcname}_norm_vcfeval-filt.dir"
    normvcftxt="${calldir}/${gsample}_germline.${vcname}_norm.vcf"
    normvcf="${normvcftxt}.gz"
}

function setup_hcparam {
    if [ $(samtools view -H "${1}" | grep -cP "^@SQ\tSN:${2}\t") -gt 0 ]; then
        hctnameparam=" -L ${chrom} "
    else
        hctnameparam=""
    fi
    hctnames=$(samtools view -H "${1}" | grep -P "^@SQ\tSN:\t" | awk '{print $2}' | awk -F":" '{print $2}')
}

function run_germline {
    if [ $1 == "help" ]; then
        echo "
        The function run_germline is used to run germline vars
        arg1: variant-caller name (run-uvc1-all, run-uvc1-call, run-uvc1-eval, ... (uvc1 can be replaced by gatk4hc, strelka2))
        arg2: output directory
        arg3: reference FASTA file indexed with .fai
        arg4: germline groundtruth vcf indexed with .tbi
        arg5: germline groundtruth bed
        arg6: germline bam indexed with .bai
        arg7: germline sample name
        arg8: number of threads to use
        "
        return 0
    fi
    calldir=$(echo "${2}" | awk -F"," '{print $1}')
    evalsuffix=$(echo "${2}" | awk -F"," '{print $2}')
    
    fastaref="${3}"
    truthvcf="${4}"
    truthbed="${5}"
    gbam="${6}"
    gsample=$(echo "${7}" | awk -F"," '{print $1}')
    grefmat=$(echo "${7}" | awk -F"," '{print $2}')
    ncpus="${8}"
    
    sdfref="${fastaref}.sdf"
    refsample=$(bcftools view -h "${truthvcf}" | tail -n1 | awk '{print $NF}')
    ${mkdir777} "${calldir}"
    
    setup_germline "${PARAMSET}"
    vcname=uvc1
    if isfound "${1}" "enable-${vcname}-all" ; then
        callvcf="${calldir}/${gsample}_germline.${PARAMSET}.vcf.gz"
        prenormvcf="${calldir}/${gsample}_germline.${PARAMSET}_prenorm.vcf.gz"
        surrogatevcf="${calldir}/${gsample}_germline.${PARAMSET}_surrogate.vcf.gz"
        surrnormvcf="${calldir}/${gsample}_germline.${PARAMSET}_surrogate_norm.vcf.gz"
        
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-call" ; then
            rm -r "${callvcf}" || true
            "${UVC}" --outvar-flag 0x1F -t "${ncpus}" -f "${fastaref}" "${gbam}" -s "${grefmat}" -o "${callvcf}"
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-norm" ; then
            bcftools view  --threads ${ncpus} -i "GERMLINE=1 || ADDITIONAL_INDEL_CANDIDATE=1" "${callvcf}" -Oz -o "${prenormvcf}"
            bcftools index --threads ${ncpus} -ft "${prenormvcf}"
            bcftools view  --threads ${ncpus} -i "GERMLINE=1" "${prenormvcf}" -Oz -o "${normvcf}"
            bcftools index --threads ${ncpus} -ft "${normvcf}"
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-norm|run-${vcname}-surrogate-align" ; then
            bash "${UVCROOT}/uvcSurrogateAlign.sh" "${surrogatevcf}" "${prenormvcf}" "${fastaref}" "${gbam}"
            bcftools view --threads ${ncpus} -i "GERMLINE=1" -Oz -o "${surrnormvcf}" "${surrogatevcf}"
            bcftools index --threads ${ncpus} -ft "${surrnormvcf}"
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-eval" ; then 
            rm -r "${evaladir}" "${evalfdir}" || true
            date ; ${monitortime} ${vcfeval1germ} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" -f GQ
            generate_disconcord_summary "${normvcf}" "${evaladir}"
            date ; ${monitortime} ${vcfeval1germ} -o "${evalfdir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${surrnormvcf}" -f GQ
            generate_disconcord_summary "${normvcf}" "${evalfdir}"
        fi
    fi
    
    setup_germline gatk4hc
    if isfound "${1}" "enable-${vcname}-all" ; then
        callvcftxt="${calldir}/${gsample}.${vcname}.vcf.gz"
        setup_hcparam "${gbam}" "${chrom}"
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-call" ; then
            date ; ${monitortime} ${gatk4lowmem} HaplotypeCaller -I "${gbam}" -O "${callvcftxt}" -R "${fastaref}" ${hctnameparam}
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-norm" ; then
            bcftools view -Oz -o "${normvcf}" "${callvcftxt}"
            bcftools index -ft "${normvcf}"
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-eval" ; then
            rm -r "${evaladir}" "${evalfdir}" || true
            date ; ${monitortime} ${vcfeval1germ} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" -f GQ
            generate_disconcord_summary "${normvcf}" "${evaladir}"
            date ; ${monitortime} ${vcfeval2germ} -o "${evalfdir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" -f GQ
            generate_disconcord_summary "${normvcf}" "${evalfdir}"
        fi
    fi
    
    setup_germline strelka2
    if isfound "${1}" "enable-${vcname}-all" ; then
        callvcfdir="${calldir}/${vcname}/"
        
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-call" ; then
            rm -r "${callvcfdir}" || true
            date ; ${monitortime} "${strelka2}" --referenceFasta="${fastaref}" --runDir="${callvcfdir}" --bam="${gbam}"
            date ; ${monitortime} "${callvcfdir}/runWorkflow.py" -j "${ncpus}" -m local 1> "${callvcfdir}/strelka2-stdout.log" 2> "${callvcfdir}/strelka2-stderr.log"
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-norm" ; then
            bcftools view -Oz -o "${normvcf}" "${callvcfdir}/results/variants/variants.vcf.gz" 
            bcftools index -ft "${normvcf}"
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-eval" ; then
            rm -r "${evaladir}" "${evalfdir}" || true
            date ; ${monitortime} ${vcfeval1germ} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" -f GQ
            generate_disconcord_summary "${normvcf}" "${evaladir}"
            date ; ${monitortime} ${vcfeval2germ} -o "${evalfdir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" -f GQ
            generate_disconcord_summary "${normvcf}" "${evalfdir}"
        fi
    fi
    
    setup_germline freebayes
    if isfound "${1}" "enable-${vcname}-all" ; then
        callvcf="${calldir}/${gsample}.${vcname}.vcf.gz"
        
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-call" ; then
            date ; ${monitortime} "${freebayes}" -f "${fastaref}" --bam="${gbam}" | bcftools view -Oz -o "${callvcf}"
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-norm" ; then
            cp "${callvcf}" "${normvcf}"
            bcftools index -ft "${normvcf}"
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-eval" ; then
            rm -r "${evaladir}" "${evalfdir}" || true
            date ; ${monitortime} ${vcfeval1germ} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" -f GQ
            generate_disconcord_summary "${normvcf}" "${evaladir}"
            date ; ${monitortime} ${vcfeval2germ} -o "${evalfdir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" -f GQ
            generate_disconcord_summary "${normvcf}" "${evalfdir}"
        fi
    fi
    
    setup_germline bcftools
    if isfound "${1}" "enable-${vcname}-all" ; then
        callvcf="${calldir}/${gsample}.${vcname}.vcf.gz"
        
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-call" ; then
            date ; ${monitortime} bcftools mpileup --threads "${ncpus}" -Ou -f "${fastaref}" "${gbam}" \
                    | bcftools call --threads "${ncpus}" -f GQ,GP -m -v -Oz -o "${callvcf}" -
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-norm" ; then
            cp "${callvcf}" "${normvcf}"
            bcftools index -ft "${normvcf}"
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-eval" ; then
            rm -r "${evaladir}" "${evalfdir}" || true
            date ; ${monitortime} ${vcfeval1germ} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" -f GQ
            date ; ${monitortime} ${vcfeval2germ} -o "${evalfdir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" -f GQ
        fi
    fi

}

