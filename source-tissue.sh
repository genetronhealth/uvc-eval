
strelka2="${EVALROOT}/tools/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py"
lofreq="${EVALROOT}/tools/lofreq-2.1.5/dist/lofreq_star-2.1.5/src/lofreq/lofreq"
somaticsniper="${EVALROOT}/tools/somatic-sniper-1.0.5.0/build/bin/bam-somaticsniper"
lolopicker="${EVALROOT}/tools/LoLoPicker/build/scripts-2.7/LoLoPicker_somatic.py"
varscan2="${java8} -jar ${EVALROOT}/tools/VarScan.v2.4.2.jar "
outlyzer="${EVALROOT}/tools/outLyzer-3.0/outLyzer.py"
octopus="${EVALROOT}/tools/venv/bin/octopus"
lancet="${EVALROOT}/tools/lancet-1.1.0/lancet"

function ensure_valid_header_for_vcfeval {
    if [ $(bcftools view -H "${1}" | wc -l) -eq 0 ]; then
        mv "${1}" "${1}.bak"
        #${UVC} "/only-print-vcf-header/" | bcftools view -Oz -o "${1}" -
        bcftools reheader -f "${2}.fai" "${1}.bak" | bcftools view -Oz -o "${1}"
        bcftools index -ft "${1}"
    fi
}

function setup_tnpair {
    vcname="${1}"
    evaladir="${calldir}/${tsample}_${nsample}.${vcname}_norm_vcfeval-all${evalsuffix}.dir"
    evalfdir="${calldir}/${tsample}_${nsample}.${vcname}_norm_vcfeval-filt${evalsuffix}.dir"
    normvcftxt="${calldir}/${tsample}_${nsample}.${vcname}_norm.vcf"
    normvcf="${normvcftxt}.gz"
    callvcfdir="${calldir}/${tsample}_${nsample}.${vcname}.vcfdir"
}

function setup_tonly {
    vcname="${1}"
    evaladir="${calldir}/${tsample}_tonly.${vcname}_norm_vcfeval-all${evalsuffix}.dir"
    evalfdir="${calldir}/${tsample}_tonly.${vcname}_norm_vcfeval-filt${evalsuffix}.dir"
    normvcftxt="${calldir}/${tsample}_tonly.${vcname}_norm.vcf"
    normvcf="${normvcftxt}.gz"
    germisecdir="${calldir}/${tsample}_tonly.${vcname}-norm.isecdir"
    nogermvcf="${germisecdir}/0001.vcf.gz"
    callvcfdir="${calldir}/${tsample}_tonly.${vcname}.vcfdir"
}

function setup_mutect2param {
    if [ $(samtools view -H "${1}" | grep -cP "^@SQ\tSN:${2}\t") -gt 0 ]; then
        mutect2tnameparam=" -L ${chrom} "
    else
        mutect2tnameparam=""
    fi
    mutect2tnames=$(samtools view -H "${1}" | grep -P "^@SQ\tSN" | awk '{print $2}' | awk -F":" '{print $2}')
}

function run_tnpair {
    if [ $1 == "help" ]; then
        echo "
        The function run_tnpair is used for evaluating a somatic variant caller with tumor-normal pair of bams
        arg1: variant-caller name (run-uvc1-all, run-uvc1-call, run-uvc1-eval, ... (uvc1 can be replaced by gatk4mutect2, strelka2, lofreq, lolopicker, varscan2, somaticsniper))
        arg2: output directory
        arg3: reference FASTA file indexed with .fai
        arg4: somatic groundtruth vcf indexed with .tbi
        arg5: somatic groundtruth bed
        arg6: tumor  bam indexed with .bai
        arg7: tumor  sample name
        arg8: normal bam indexed with .bai
        arg9: normal sample name
        arg10: number of threads to use
        "
        return 0
    fi
    calldir=$(echo "${2}" | awk -F"," '{print $1}')
    evalsuffix=$(echo "${2}" | awk -F"," '{print $2}')
    
    fastaref="${3}"
    truthvcf="${4}"
    truthbed="${5}"
    tbam="${6}"
    tsample=$(echo "${7}" | awk -F"," '{print $1}')
    trefmat=$(echo "${7},${7}" | awk -F"," '{print $2}')
    nbam="${8}"
    nsample=$(echo "${9}" | awk -F"," '{print $1}')
    nrefmat=$(echo "${9},${9}" | awk -F"," '{print $2}')
    ncpus="${10}"
    
    sdfref="${fastaref}.sdf"
    refsample=$(bcftools view -h "${truthvcf}" | tail -n1 | awk '{print $10}')
    ${mkdir777} "${calldir}"
    
    setup_tnpair "${PARAMSET}"
    vcname=uvc1
    if isfound "${1}" enable-${vcname}-all ; then
        if isnotfound "${1}" "enable-only-vcf-eval" ; then
            prevcallvcf="${calldir}/${tsample}_${nsample}.${Tparamset}.vcf.gz"
            callvcf="${calldir}/${tsample}_${nsample}.${PARAMSET}.vcf.gz"
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-call" ; then
                
                # use tumor vcf from previous run if available
                tumorvcf1="${prevcallvcf}.byproduct/${tsample}_uvc1.vcf.gz"
                tumorvcf2="${callvcf}.byproduct/${tsample}_uvc1.vcf.gz"
                ${mkdir777} "${callvcf}.byproduct"
                cp -s "${tumorvcf1}" "${tumorvcf2}" || true
                cp -s "${tumorvcf1}.tbi" "${tumorvcf2}.tbi" || true
                
                date ; ${monitortime} "${UVCTN}" "${fastaref}" "${tbam}" "${nbam}" "${callvcf}" "${trefmat},${nrefmat}" -t ${ncpus} ${uvcparams}
            fi
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-norm" ; then
                date ; ${monitortime} "${UVCNORM}" "${callvcf}" "${normvcf}"
            fi
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-eval" ; then
            rm -r "${evaladir}" "${evalfdir}" || true
            date ; ${monitortime} ${vcfeval1soma} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},${trefmat}" -f QUAL
            chmod uga+rwxs "${evaladir}" "${evalfdir}" || true
        fi
    fi
    
    setup_tnpair gatk4mutect2
    if  isfound "${1}" enable-${vcname}-all ; then
        callvcf="${calldir}/${tsample}_${nsample}.${vcname}.vcf"
        #callvcfdir="${calldir}/${tsample}_${nsample}.${vcname}.vcfdir"
        if isnotfound "${1}" "enable-only-vcf-eval" ; then
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-call" ; then
                setup_mutect2param "${tbam}" "${chrom}"
                if [ ${nchroms} -gt 1 -a -z "${mutect2tnameparam}" ]; then
                    ${mkdir777} "${callvcfdir}"
                    for tname in ${mutect2tnames}; do
                        printf "${monitortime} ${gatk4lowmem} Mutect2 -R \"${fastaref}\" -I \"${tbam}\" -I \"${nbam}\" -O \"${callvcfdir}/${tname}.vcf\" "
                        printf " --tumor \"${trefmat}\" --normal \"${nrefmat}\" ${mutect2tnameparam} "
                        printf " && bcftools view -Oz -o \"${callvcfdir}/${tname}.vcf.gz\" \"${callvcfdir}/${tname}.vcf\" "
                        printf " && bcftools index -ft \"${callvcfdir}/${tname}.vcf\" \n"
                    done > "${callvcfdir}/${tsample}_${nsample}.${vcname}.sh"
                    cat "${callvcfdir}/${tsample}_${nsample}.${vcname}.sh" | parallel -j ${nchroms} 
                    bcftools merge -Ov -o "${callvcf}" "${callvcfdir}/"*".vcf.gz"
                else
                    date ; ${monitortime} ${gatk4lowmem} Mutect2 -R "${fastaref}" -I "${tbam}" -I "${nbam}"  -O "${callvcf}" \
                        --tumor "${trefmat}" --normal "${nrefmat}" ${mutect2tnameparam}
                fi
            fi
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-norm" ; then
                date ; ${monitortime} ${gatk4lowmem} FilterMutectCalls --reference "${fastaref}" -V "${callvcf}" -O "${normvcftxt}"
                txtvcf2idxvcf "${normvcftxt}" "${normvcf}"
            fi
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-eval" ; then
            rm -r "${evaladir}" "${evalfdir}" || true
            date ; ${monitortime} ${vcfeval1soma} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},${trefmat}" -f INFO.TLOD
            date ; ${monitortime} ${vcfeval2soma} -o "${evalfdir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},${trefmat}" -f INFO.TLOD
            chmod uga+rwxs "${evaladir}" "${evalfdir}" || true
        fi
    fi
    
    setup_tnpair strelka2
    if isfound "${1}" enable-${vcname}-all ; then
        #callvcfdir="${calldir}/${tsample}_${nsample}.${vcname}.vcfdir"
        if isnotfound "${1}" "enable-only-vcf-eval" ; then
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-call" ; then
                rm -r "${callvcfdir}" || true
                date ; ${monitortime} "${strelka2}" --referenceFasta="${fastaref}" --runDir="${callvcfdir}" --tumorBam="${tbam}" --normalBam="${nbam}"
                date ; ${monitortime} "${callvcfdir}/runWorkflow.py" -j "${ncpus}" -m local 1> "${callvcfdir}/stdout-default.log" 2> "${callvcfdir}/stderr-default.log"
            fi
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-norm" ; then
                bcftools concat -a "${callvcfdir}/results/variants/somatic.snvs.vcf.gz" "${callvcfdir}/results/variants/somatic.indels.vcf.gz" \
                    | awk 'OFS="\t" { if ($0 !~ "^#") {$9="GT:"$9; $10="0/1:"$10; $11="0/1:"$11; print; } else { print;}}' \
                    | bcftools view -Oz -o "${normvcf}"
                bcftools index -ft "${normvcf}"
            fi
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-eval" ; then
            rm -r "${evaladir}" "${evalfdir}" || true
            date ; ${monitortime} ${vcfeval1soma} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},TUMOR" -f INFO.SomaticEVS
            date ; ${monitortime} ${vcfeval2soma} -o "${evalfdir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},TUMOR" -f INFO.SomaticEVS
            chmod uga+rwxs "${evaladir}" "${evalfdir}" || true
        fi
    fi

    setup_tnpair lofreq
    if isfound "${1}" enable-${vcname}-all ; then
        #callvcfdir="${calldir}/${vcname}/"
        if isnotfound "${1}" "enable-only-vcf-eval" ; then
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-call" ; then
                rm -r "${callvcfdir}" || true
                ${mkdir777} "${callvcfdir}" 
                ${monitortime} "${lofreq}" somatic -t "${tbam}" -n "${nbam}" -o "${callvcfdir}/" -f "${fastaref}" --threads ${ncpus} --call-indels
            fi
            if  isfound "${1}" "run-${vcname}-all|run-${vcname}-norm" ; then
                bcftools concat -a "${callvcfdir}/somatic_final.snvs.vcf.gz" "${callvcfdir}/somatic_final.indels.vcf.gz" \
                    | awk '{if ($1 ~ "^##") { print $0; } else {if ($1 ~ "^#") { print $0 "\tFORMAT\tTUMOR\tNORMAL" } else { print $0 "\tGT\t0/1\t0/1" } }}' \
                    | bcftools view -Oz -o "${normvcf}"
                ensure_valid_header_for_vcfeval "${normvcf}" "${fastaref}"
                bcftools index -ft "${normvcf}"
            fi
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-eval" ; then 
            rm -r "${evaladir}" "${evalfdir}" || true
            date ; ${monitortime} ${vcfeval1soma} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},TUMOR" -f QUAL
            date ; ${monitortime} ${vcfeval2soma} -o "${evalfdir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},TUMOR" -f QUAL
            chmod uga+rwxs "${evaladir}" "${evalfdir}" || true
        fi
    fi

    setup_tnpair lolopicker
    if isfound "${1}" enable-${vcname}-all ; then
        #callvcfdir="${calldir}/${vcname}/"
        callvcf="${callvcfdir}/raw_somatic_variants.vcf"
        if isnotfound "${1}" "enable-only-vcf-eval" ; then
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-call" ; then
                ${mkdir777} "${callvcfdir}"
                date ; ${monitortime} python "${lolopicker}" -t "${tbam}" -n "${nbam}" -r "${fastaref}" -b "${truthbed}" -o "${callvcfdir}"
            fi
            if  isfound "${1}" "run-${vcname}-all|run-${vcname}-norm" ; then
                printf "##fileformat=VCFv4.2\n" > "${callvcf}"
                printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${trefmat}\n" >> "${callvcf}"
                cat "${callvcfdir}/raw_somatic_var"*"ants.txt" | grep pass_to_test \
                    | awk 'OFS="\t" {print $1, $2, ".", $3, $4, "200", "PASS", "t_ref="$5";t_alt="$6";n_ref="$7";n_alt="$8";judge="$9, "GT", "0/1"}' \
                    | sort -k1,2n >> "${callvcf}"
                bcftools view -Oz -o "${normvcf}" "${callvcf}"
                ensure_valid_header_for_vcfeval "${normvcf}" "${fastaref}"
                bcftools index -ft "${normvcf}"
            fi
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-eval" ; then 
            rm -r "${evaladir}" "${evalfdir}" || true
            date ; ${monitortime} ${vcfeval1soma} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},${trefmat}" -f QUAL
            date ; ${monitortime} ${vcfeval2soma} -o "${evalfdir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},${trefmat}" -f QUAL
            chmod uga+rwxs "${evaladir}" "${evalfdir}" || true
        fi
    fi
    
    setup_tnpair octopus
    if isfound "${1}" enable-${vcname}-all ; then
        #callvcfdir="${calldir}/${vcname}/"
        callvcf="${callvcfdir}/octopus.vcf"
        if isnotfound "${1}" "enable-only-vcf-eval" ; then
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-call" ; then
                ${mkdir777} "${callvcfdir}"
                abam="${callvcfdir}/${tsample}_${nsample}_merged.bam"
                samtools merge -f -@ ${ncpus} "${abam}" "${tbam}" "${nbam}"
                samtools index -@ ${ncpus} "${abam}"
                date ; ${monitortime} "${octopus}" --target-read-buffer-footprint 1GB -I "${abam}" -N "${nrefmat}" -R "${fastaref}.fa" -o "${callvcf}" --threads ${ncpus} --somatics-only
            fi
            if  isfound "${1}" "run-${vcname}-all|run-${vcname}-norm" ; then
                bcftools view "${callvcf}" | awk 'OFS="\t" { if ($0 !~ "^#") {$9="GT:"$9; $10="0/1:"$10; $11="0/1:"$11; print; } else { print;}}' \
                    | sed 's/\tGT:GT:/\tGT:GTX:/g' | bcftools view -Oz -o "${normvcf}" - 
                ensure_valid_header_for_vcfeval "${normvcf}" "${fastaref}"
                bcftools index -ft "${normvcf}"
            fi
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-eval" ; then 
            rm -r "${evaladir}" "${evalfdir}" || true
            date ; ${monitortime} ${vcfeval1soma} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},${trefmat}" -f QUAL
            date ; ${monitortime} ${vcfeval2soma} -o "${evalfdir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},${trefmat}" -f QUAL
            chmod uga+rwxs "${evaladir}" "${evalfdir}" || true
        fi
    fi
    
    setup_tnpair lancet
    if isfound "${1}" enable-${vcname}-all ; then
        callvcf="${callvcfdir}/lancet.vcf"
        if isnotfound "${1}" "enable-only-vcf-eval" ; then
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-call" ; then
                ${mkdir777} "${callvcfdir}"
                date ; ${monitortime} "${lancet}" -t "${tbam}" -n "${nbam}" -r "${fastaref}" --num-threads ${ncpus} -B "${truthbed}" > "${callvcf}"
            fi
            if  isfound "${1}" "run-${vcname}-all|run-${vcname}-norm" ; then
                bcftools view "${callvcf}" -Oz -o "${normvcf}"
                ensure_valid_header_for_vcfeval "${normvcf}" "${fastaref}"
                bcftools index -ft "${normvcf}"
            fi
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-eval" ; then 
            rm -r "${evaladir}" "${evalfdir}" || true
            date ; ${monitortime} ${vcfeval1soma} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},${trefmat}" -f QUAL
            date ; ${monitortime} ${vcfeval2soma} -o "${evalfdir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},${trefmat}" -f QUAL
            chmod uga+rwxs "${evaladir}" "${evalfdir}" || true
        fi
    fi
    
    setup_tnpair varscan2
    if isfound "${1}" enable-${vcname}-all ; then
        callprefix="${calldir}/${tsample}_${nsample}.${vcname}"
        if isnotfound "${1}" "enable-only-vcf-eval" ; then
            if  isfound "${1}" "run-${vcname}-all|run-${vcname}-call" ; then
                date ; ${monitortime} ${varscan2} somatic <(samtools mpileup -f "${fastaref}" "${nbam}") <(samtools mpileup -f "${fastaref}" "${tbam}") "${callprefix}" --output-vcf
            fi
            filtvcftxt="${callprefix}.fpfilter.vcf"
            filtvcf="${filtvcftxt}.gz"
            if  isfound "${1}" "run-${vcname}-all|run-${vcname}-norm" ; then
                txtvcf2idxvcf "${callprefix}.snp.vcf" "${callprefix}.snp.vcf.gz" 
                txtvcf2idxvcf "${callprefix}.indel.vcf" "${callprefix}.indel.vcf.gz" 
                date ; ${monitortime} bcftools concat -a "${callprefix}.snp.vcf.gz" "${callprefix}.indel.vcf.gz" | bcftools view -i "SS='2'" -Oz -o "${normvcf}"
                date ; ${monitortime} bcftools index -ft "${normvcf}"
                ensure_valid_header_for_vcfeval "${normvcf}" "${fastaref}"
                bcftools view "${normvcf}" -Ov -o "${normvcftxt}" 
                
                # apply additional false-positive filter # commented out due to lack of documentation from varscan2
                # ${varscan2} readcounts <(samtools mpileup -f "${fastaref}" "${tbam}") --variants-file "${normvcftxt}" --output-file "${normvcftxt}.readcounts"
                # ${varscan2} fpfilter "${normvcftxt}" "${normvcftxt}.readcounts" --output-file "${filtvcftxt}"
                # txtvcf2idxvcf "${filtvcftxt}" "${filtvcf}"
            fi
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-eval" ; then 
            rm -r "${evaladir}" "${evalfdir}" || true
            date ; ${monitortime} ${vcfeval1soma} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},TUMOR" -f INFO.SSC
            date ; ${monitortime} ${vcfeval2soma} -o "${evalfdir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},TUMOR" -f INFO.SSC
            chmod uga+rwxs "${evaladir}" "${evalfdir}" || true
            
            # evaluate additional false-positive filter
            # varscan2evaldir="${evaladir/varscan2/varscan2fpfilter}"
            # normvcf="${filtvcf}"
            # rm -r "${evaladir}" "${evalfdir}" || true
            # date ; ${monitortime} ${vcfeval1soma} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},TUMOR" -f QUAL
            # date ; ${monitortime} ${vcfeval2soma} -o "${evalfdir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},TUMOR" -f QUAL
            # chmod uga+rwxs "${evaladir}" "${evalfdir}" || true
        fi
    fi
    
    setup_tnpair somaticsniper
    if isfound "${1}" enable-${vcname}-all ; then
        callvcf="${calldir}/${tsample}_${nsample}.${vcname}.vcf"
        if isnotfound "${1}" "enable-only-vcf-eval" ; then
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-call" ; then
                date ; ${monitortime} "${somaticsniper}" -F vcf -f "${fastaref}" "${tbam}" "${nbam}" "${callvcf}"
            fi
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-norm" ; then
                date ; ${monitortime} bcftools view "${callvcf}" -i "SS=2" -Oz -o "${normvcf}"
                date ; ${monitortime} bcftools index -ft "${normvcf}"
            fi
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-eval" ; then
            rm -r "${evaladir}" "${evalfdir}" || true
            date ; ${monitortime} ${vcfeval1soma} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},TUMOR" -f FORMAT.SSC
            date ; ${monitortime} ${vcfeval2soma} -o "${evalfdir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${normvcf}" --sample "${refsample},TUMOR" -f FORMAT.SSC
            chmod uga+rwxs "${evaladir}" "${evalfdir}" || true
        fi
    fi
}

function run_tonly {
    if [ $1 == "help" ]; then
        echo "
        The function run_tonly is used for evaluating a variant caller with tumor-only bam without distinguishing germline-vs-somatic origin
        arg1: variant-caller name (run-uvc1-all, run-uvc1-call, run-uvc1-eval, ... (uvc1 can be replaced by gatk4mutect2, lofreq, varscan2))
        arg2: output directory
        arg3: reference FASTA file indexed with .fai
        arg4: somatic groundtruth vcf indexed with .tbi
        arg5: somatic groundtruth bed
        arg6: tumor  bam indexed with .bai
        arg7: tumor  sample name
        arg8: germline groundtruth vcf indexed with .tbi
        arg9: number of threads to use
        "
        return 0
    fi
    
    calldir=$(echo "${2}" | awk -F"," '{print $1}')
    evalsuffix=$(echo "${2}" | awk -F"," '{print $2}')
    
    fastaref="${3}"
    truthvcf="${4}"
    truthbed="${5}"
    tbam="${6}"
    tsample=$(echo "${7}" | awk -F "," '{print $1}')
    trefmat=$(echo "${7},${7}" | awk -F "," '{print $2}')
    germvcf="${8}"
    ncpus="${9}"
    sdfref="${fastaref}.sdf"
    refsample=$(bcftools view -h "${truthvcf}" | tail -n1 | awk '{print $10}')
    ${mkdir777} "${calldir}"
    
    setup_tonly "${PARAMSET}"
    vcname=uvc1
    if isfound "${1}" "enable-${vcname}-all" ; then
        callvcf="${calldir}/${tsample}_tonly.${PARAMSET}.vcf.gz"
        # callvcf="${calldir}/${tsample}_tonly.${vcname}.vcf.gz"
        if isnotfound "${1}" "enable-only-vcf-eval" ; then
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-call" ; then
                date ; ${monitortime} "${UVC}" -f "${fastaref}" "${tbam}" -o "${callvcf}" -s "${trefmat}"
            fi
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-norm" ; then
                bcftools index -ft "${callvcf}"
                date ; ${monitortime} "${UVCNORM}" "${callvcf}" "${normvcf}"
                date ; ${monitortime} bcftools isec -Oz ${GROUNDTRUTH_ISEC_FLAGS} -p "${germisecdir}" "${germvcf}" "${normvcf}"
                date ; ${monitortime} bcftools index -ft "${nogermvcf}"
            fi
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-eval" ; then
            rm -r "${evaladir}" "${evalfdir}" || true
            date ; ${monitortime} ${vcfeval1soma} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${nogermvcf}" -f QUAL
            chmod uga+rwxs "${evaladir}" "${evalfdir}" || true
        fi
    fi

    setup_tonly gatk4mutect2
    if isfound "${1}" "enable-${vcname}-all" ; then
        callvcf="${calldir}/${tsample}_tonly.${vcname}.vcf"
        if isnotfound "${1}" "enable-only-vcf-eval" ; then
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-call" ; then
                setup_mutect2param "${tbam}" "${chrom}"
                date ; ${monitortime} ${gatk4lowmem} Mutect2 -R "${fastaref}" -I "${tbam}" -O "${callvcf}" --tumor "${trefmat}" ${mutect2tnameparam}
            fi
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-norm" ; then
                date ; ${monitortime} ${gatk4lowmem} FilterMutectCalls --reference "${fastaref}" -V "${callvcf}" -O "${normvcftxt}"
                txtvcf2idxvcf "${normvcftxt}" "${normvcf}"
                date ; ${monitortime} bcftools isec ${GROUNDTRUTH_ISEC_FLAGS} -Oz -p "${germisecdir}" "${germvcf}" "${normvcf}"
                date ; ${monitortime} bcftools index -ft "${nogermvcf}"
            fi
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-eval" ; then
            rm -r "${evaladir}" "${evalfdir}" || true
            date ; ${monitortime} ${vcfeval1soma} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${nogermvcf}" --sample "${refsample},${trefmat}" -f INFO.TLOD
            date ; ${monitortime} ${vcfeval2soma} -o "${evalfdir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${nogermvcf}" --sample "${refsample},${trefmat}" -f INFO.TLOD
        fi
    fi
    
    setup_tonly lofreq
    if isfound "${1}" "enable-${vcname}-all" ; then
        callvcf="${calldir}/${tsample}_tonly.${vcname}.vcf.gz"
        if isnotfound "${1}" "enable-only-vcf-eval" ; then
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-call" ; then
                rm -r "${callvcf}" || true 
                "${lofreq}" call-parallel --pp-threads ${ncpus} "${tbam}" -o "${callvcf}" -f "${fastaref}" --call-indels
            fi
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-norm" ; then
                bcftools view "${callvcf}" \
                    | awk '{if ($1 ~ "^##") { print $0; } else {if ($1 ~ "^#") { print $0 "\tFORMAT\tTUMOR" } else { print $0 "\tGT\t0/1" } }}' \
                    | bcftools view -Oz -o "${normvcf}"
                bcftools index -ft "${normvcf}"
                date ; ${monitortime} bcftools isec ${GROUNDTRUTH_ISEC_FLAGS} -Oz -p "${germisecdir}" "${germvcf}" "${normvcf}"
                date ; ${monitortime} bcftools index -ft "${nogermvcf}"
                ensure_valid_header_for_vcfeval "${nogermvcf}" "${fastaref}"
            fi
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-eval" ; then 
            rm -r "${evaladir}" "${evalfdir}" || true
            date ; ${monitortime} ${vcfeval1soma} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${nogermvcf}" -f QUAL
            date ; ${monitortime} ${vcfeval2soma} -o "${evalfdir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${nogermvcf}" -f QUAL
        fi
    fi
    
    setup_tonly outlyzer
    if isfound "${1}" "enable-${vcname}-all" ; then
        #callvcfdir="${calldir}/${tsample}_tumoroly.${vcname}.vcfdir"
        mkdir -p "${callvcfdir}/"
        callvcf="${callvcfdir}/resultsExample.vcf"
        if isnotfound "${1}" "enable-only-vcf-eval" ; then
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-call" ; then
                rm -r "${callvcf}" || true 
                python3 "${outlyzer}" calling -core "${ncpus}"  -bed "${truthbed}" -bam "${tbam}" -output "${callvcfdir}/" -ref "${fastaref}"
            fi
            if isfound "${1}" "run-${vcname}-all|run-${vcname}-norm" ; then
                bcftools view "${callvcf}" -Oz -o "${normvcf}"
                bcftools index -ft "${normvcf}"
                date ; ${monitortime} bcftools isec ${GROUNDTRUTH_ISEC_FLAGS} -Oz -p "${germisecdir}" "${germvcf}" "${normvcf}"
                date ; ${monitortime} bcftools index -ft "${nogermvcf}"
                ensure_valid_header_for_vcfeval "${nogermvcf}" "${fastaref}"
            fi
        fi
        if isfound "${1}" "run-${vcname}-all|run-${vcname}-eval" ; then 
            rm -r "${evaladir}" "${evalfdir}" || true
            date ; ${monitortime} ${vcfeval1soma} -o "${evaladir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${nogermvcf}" -f QUAL
            date ; ${monitortime} ${vcfeval2soma} -o "${evalfdir}" -t "${sdfref}" -b "${truthvcf}" -e "${truthbed}" -c "${nogermvcf}" -f QUAL
        fi
    fi

}

