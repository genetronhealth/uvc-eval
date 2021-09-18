set -evx

# isfound text pattern
function isfound {
    [ $(echo "${1}" | grep -cP "(\.|^)${2}(\.|$)") -gt 0 ]
}

function isnotfound {
    [ $(echo "${1}" | grep -cP "(\.|^)${2}(\.|$)") -eq 0 ]
}

function issubstr {
    [ $(echo "${1}" | grep -cP "${2}") -gt 0 ]
}

function logversion {
    bcftools view -h "${1}" | grep -F "##variantCallerVersion" | awk '{print $1}' | awk -F"." '{print substr($(NF), 0, 7)}' > "${2}.version.info"
}

function vcfsubtract {
    bcftools isec -Oz -p "${1}.isecdir" "${1}" "${2}" -T "${3}"
    ${mkdir777} -p "${1}.isecdir/backup/"
    mv "${1}" "${1}.isecdir/backup/"
    mv "${1}.tbi" "${1}.isecdir/backup/"
    cp "${1}.isecdir/0000.vcf.gz" "${1}"
    bcftools index -ft "${1}"
}

function txtvcf2idxvcf {
    bcftools view "${1}" -Oz -o "${1}.gz"
    bcftools index -ft "${1}.gz"
}

function run_index_fastaref {
    if isfound "${1}" run-faidx-index ; then
        date ; time -p samtools faidx "${2}"
    fi
    if isfound "${1}" run-gatk4-index ; then
        date ; time -p ${gatk4lowmem} CreateSequenceDictionary --REFERENCE "${2}"
    fi
    if isfound "${1}" run-bwa-index ; then
        date ; time -p bwa index "${2}" 
    fi
    if isfound "${1}" run-rtg-index ; then
        rm -r "${1}.sdf"
        date ; time -p "${java8}" -jar "${RTGJAR}" format -o "${1}.sdf" "${2}"
    fi
}

SCRIPTDIR=$(dirname $(which $0))
if [ -z "${EVALROOT}" ]; then
    EVALROOT=$(cd "${SCRIPTDIR}/../" && pwd)
fi

monitortime="${EVALROOT}/tools/time"

if [ -e "${monitortime}" ]; then
    echo "The executable ${monitortime} is found. "
    monitortime="${monitortime} -v "
else
    echo "The executable ${monitortime} is not found. Falling back to shell default. "
    monitortime="time -p "
fi

mkdir777="mkdir -m uga+rwxs -p "
java8="${EVALROOT}/tools/jre/bin/java"
parallel="${EVALROOT}/tools/parallel"
fastqdump="${SOFT}/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump"

BWA="${EVALROOT}/tools/bwa-0.7.17/bwa"
if [ -z "${UVCROOT}" ]; then
    UVCROOT="${EVALROOT}/tools/uvc/"
fi
UVC="${UVCROOT}/uvc1"
UVCTN="${UVCROOT}/uvcTN.sh"
UVCNORM="${UVCROOT}/uvcnorm.sh"

export PATH="${UVCROOT}:${EVALROOT}/tools/jre/bin:${EVALROOT}/bin:${EVALROOT}/tools:${PATH}"
export TMPDIR="${EVALROOT}/systmp/"

RTGJAR="${EVALROOT}/tools/rtg-tools-3.11/RTG.jar"
vcfeval1germ=" ${java8} -jar ${RTGJAR} vcfeval --ref-overlap --decompose --all-records "
vcfeval2germ=" ${java8} -jar ${RTGJAR} vcfeval --ref-overlap --decompose "
vcfeval1soma=" ${java8} -jar ${RTGJAR} vcfeval --ref-overlap --decompose --all-records --squash-ploidy "
vcfeval2soma=" ${java8} -jar ${RTGJAR} vcfeval --ref-overlap --decompose               --squash-ploidy "
GATK_VERSION="4.1.9.0"
gatk4lowmem=" ${java8} -Xmx4g -Djava.io.tmpdir=${EVALROOT}/systmp -jar ${EVALROOT}/tools/gatk-${GATK_VERSION}/gatk-package-${GATK_VERSION}-local.jar "
#gatkmain=" ${EVALROOT}/tools/gatk-4.1.9.0/gatk --java-options -Djava.io.tmpdir=${EVALROOT}/systmp "


HG001BED="${EVALROOT}/HG001_HG002_groundtruth/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed"
HG002BED="${EVALROOT}/HG001_HG002_groundtruth/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed"
HG005BED="${EVALROOT}/HG001_HG002_groundtruth/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.2_highconf_noMetaSV.bed"

HG001VCF="${EVALROOT}/HG001_HG002_groundtruth/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz"
HG002VCF="${EVALROOT}/HG001_HG002_groundtruth/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz"
HG005VCF="${EVALROOT}/HG001_HG002_groundtruth/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.2_highconf.vcf.gz"

# GIAB_TRUTH_VIEW_EXPR=' ([GT != "ref"][0] && [GT != "mis"][1]) && (FORMAT/GT[0] != FORMAT/GT[1]) && ([GT = "ref"][1] || [GT = "mis"][1]) '
# https://www.ncbi.nlm.nih.gov/books/NBK44417/#Content.what_is_a_reference_snp_or__rs_i
GROUNDTRUTH_ISEC_FLAGS=" -c both "
GROUNDTRUTH_NORM_FLAGS=" -m-both "

UMI_FILTER_EXPR_1=" cVQ2M[:0] - cVQ2[:1] == 0 && QUAL >= 40 "
# use QUAL thres to prevent vcfeval from skipping regions that are too complex # 844283 / 20586 
UMI_FILTER_EXPR_2="((APXM[:0] - (30 * APDP[:0])) <= 0 || TYPE != \"snp\") && (TYPE = \"snp\" || cMmQ[:1] >= 20)"

PAT="${1}"
Tparamset=$(echo "${2}" | awk -F"," '{print $1}')
PARAMSET=$(echo "${2}" | awk -F"," '{print $(2)}')
EVALOUT=$(echo "${2}" | awk -F"," '{print $(3)}')
if [ -z "${Tparamset}" -o "${Tparamset}" == '.' ]; then
    Tparamset=$(${UVC} -v | grep -i uvc | sed 's/\./-/g')
fi
if [ -z "${PARAMSET}" -o "${PARAMSET}" == '.' ]; then
    PARAMSET=${Tparamset}
fi
if [ -z "${EVALOUT}" -o "${EVALOUT}" == '.' ]; then
    EVALOUT=${PARAMSET}
fi
uvcparams="${@:3}"

ncpus=32
for i in $(seq 0 48); do
    if isfound "${PAT}" use-${i}-cpus ; then
        ncpus=${i}
    fi
done

nbams=8
for i in $(seq 0 48); do
    if isfound "${PAT}" use-${i}-bams ; then
        nbams=${i}
    fi
done

nchroms=0
for i in $(seq 0 48); do
    if isfound "${PAT}" use-${i}-chroms ; then
        nchroms=${i}
    fi
done

#REF="${EVALROOT}/datafiles/Homo_sapiens_assembly19.fasta" # This ref should not be used
HG19="${EVALROOT}/datafiles/hg19_UCSC.fa"
HS37D5="${EVALROOT}/datafiles/hs37d5.fa" # the 1000genome b37 human genome reference without chr prefix
GRCH38="${EVALROOT}/datafiles/GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.fna"

mkdir -p "${TMPDIR}" || true

function run_bamprep {
    if [ $1 == "help" ]; then
        echo "
        The function runbamprep is used for preparing bam files
        arg1: program name (run-gatk4prep (including run-gatk4rmdup and run-gatk4bqsr), run-dindel)
        arg2: reference FASTA file
        arg3: bam with index
        arg4: germline vcf file indexed with .idx
        arg5: number of threads to use, use 6 by default
        "
        return 0
    fi
    fastaref="${2}"
    rawbam="${3}"
    germvcf="${4}"
    ncpus2="${5}"
    if [ -z "${ncpus2}" ]; then
        ncpus2=6
    fi

    rmdupbam="${rawbam/%.bam/.rmdup.bam}"
    recalbam="${rawbam/%.bam/.rmdup_recal.bam}"
    dindelbam="${rawbam/%.bam/.rmdup_recal_dindel.bam}"
    rmdupmetrics="${rawbam/%.bam/.rmdup.metrics}"
    recaltable="${rawbam/%.bam/.rmdup_recal.table}"
    abra2bam="${rawbam/%.bam/.abra2.bam}"
    
    if isfound "${1}" "run-gatk4prep|run-gatk4rmdup" ; then
        date ; time -p ${gatk4lowmem} MarkDuplicates --ASSUME_SORT_ORDER coordinate --REMOVE_DUPLICATES true -I "${rawbam}" -M ${rmdupmetrics} -O "${rmdupbam}"
        date ; time -p samtools index -@ "${ncpus2}" "${rmdupbam}"
    fi
    if isfound "${1}" "run-gatk4prep|run-gatk4bqsr" ; then
        date ; time -p ${gatk4lowmem} BaseRecalibrator -I "${rmdupbam}" --known-sites "${germvcf}" -O "${recaltable}" -R "${fastaref}"
        date ; time -p ${gatk4lowmem} ApplyBQSR -bqsr "${recaltable}" -I "${rmdupbam}" -O "${recalbam}"
        date ; time -p samtools index -@ "${ncpus2}" "${recalbam}"
    fi
    if isfound "${1}" "run-dindel" ; then
        rm "${dindelbam}" || true
        date ; time -p "${lofreq}" indelqual --dindel -f "${fastaref}" -o "${dindelbam}" "${recalbam}"
        date ; time -p samtools index -@ "${ncpus2}" "${dindelbam}"
    fi
    if isfound "${1}" "run-abra2"; then
        date ; time -p "${java8}" -jar "${EVALROOT}/tools/abra2-2.22.jar" --in "${rawbam}" --out "${abra2bam}" --ref "${fastaref}" --threads "${ncpus}"
    fi
}

