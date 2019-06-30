set -evx

SCRIPTDIR=$(dirname $(which $0))
if [ -z "${EVALROOT}" ]; then
    source "${SCRIPTDIR}/../source-eval.sh"
else
    source "${EVALROOT}/source-eval.sh"
fi

truthbed1="${EVALROOT}/HG001_HG002_groundtruth/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed"
truthvcf1="${EVALROOT}/HG001_HG002_groundtruth/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz" 

bamfile1="${EVALROOT}/giab-mix/silico-rawinput/HG001.hs37d5.300x.bam"
vcfpref1="${EVALROOT}/powerlaw/data/HG001.hs37d5.300x" #${nochr}-q${mq}.vcf.gz"

truthbed12intersected="${EVALROOT}/HG001_HG002_groundtruth/intersect_HG001_HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_nosomaticdel_noinconsistent.bed"
truthvcf12merged="${EVALROOT}/HG001_HG002_groundtruth/intersect_HG001_HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_nosomaticdel_noinconsistent.union.vcf.gz"

#truthbed12intersected2="${EVALROOT}/HG001_HG002_groundtruth/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed"
#truthvcf12merged2="${EVALROOT}/HG001_HG002_groundtruth/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz"

bamfile2="${EVALROOT}/powerlaw/data/ERX1966272/BNA12878-1.read12.bam" #"${EVALROOT}/powerlaw/data/bam/SRR5683134.bam" #"${EVALROOT}/powerlaw/data/bam/SRR10134980.bam"
vcfpref2="${EVALROOT}/powerlaw/data/ERX1966272/BNA12878-1.read12.var" #"${EVALROOT}/powerlaw/data/vcf/SRR10134980" #${nochr}-q${mq}.vcf.gz"
bamfile3="${EVALROOT}/powerlaw/data/SRR7526729/SRR7526729.bam"
vcfpref3="${EVALROOT}/powerlaw/data/SRR7526729/SRR7526729.var"

truthbed=${truthbed1}
truthvcf=${truthvcf1}

if [ "$1" -eq 1 ]; then
    #truthbed=${truthbed1}
    #truthvcf=${truthvcf1}
    bamfile="${bamfile1}"
    vcfpref="${vcfpref1}"
    nbins=100
    mindp=200
    plotheader="Illumina/HiSeq/whole-genome-sequencing cell_line=HG001 min_depth_required=${mindp}x" # average_sequencing_depth=300x 
    univ_offset=1
elif [ "$1" -eq 3 ]; then
    #truthbed=${truthbed12intersected}
    #truthvcf=${truthvcf12merged}
    bamfile="${bamfile2}"
    vcfpref="${vcfpref2}"
    nbins=300
    mindp=600
    plotheader="BGI/BGISEQ/whole-exome-sequencing cell_line=HG001 min_depth_required=${mindp}x"
    univ_offset=1
elif [ "$1" -eq 2 ]; then
    truthbed=${truthbed12intersected}
    truthvcf=${truthvcf12merged}
    bamfile="${bamfile3}"
    vcfpref="${vcfpref3}"
    nbins=1000
    mindp=2000
    plotheader="Illumina/NextSeq/amplicon-sequencing cell_line=HG001_HG002_mixed_at_1_to_99 min_depth_required=${mindp}x"
    univ_offset=2
    if [ ! -f "${truthbed}" ]; then
        bedtools intersect -a "${truthbed1}" -b "${truthbed12intersected2}" > "${truthbed}"
    fi
    if [ ! -f "${truthvcf}" ]; then
        bcftools merge -Ou "${truthvcf1}" "${truthvcf12merged2}" | bcftools view -T "${truthbed}" -Oz -o "${truthvcf}" && bcftools index -ft "${truthvcf}"
    fi
else
    echo "The flag ($1) is neither 1 nor 2"
    exit -1
fi

$mkdir777 "${EVALROOT}/powerlaw/data/vcf/"

MAXDP=$((1000*1000))
for mq in 0 60; do
    
    for nochr in $(seq 1 22); do
         echo bcftools mpileup --max-depth ${MAXDP} --max-idepth ${MAXDP} --min-BQ 20 -f "${EVALROOT}/datafiles/hs37d5.fa" "${bamfile}" -Oz -o "${vcfpref}.${nochr}-q${mq}.vcf.gz"  -a AD,ADF,ADR,DP,SP -r $nochr -q ${mq}
    done
    
    for nochr in $(seq 1 22); do
        #truthvcf="${EVALROOT}/HG001_HG002_groundtruth/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz" 
        #truthbed="${EVALROOT}/HG001_HG002_groundtruth/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed"
        echo "bcftools view ${vcfpref}.${nochr}-q${mq}.vcf.gz -T ${truthbed} | python ${EVALROOT}/powerlaw/powerlaw-table.py ${nochr} ${truthvcf} $mq ${nbins} ${mindp} > ${vcfpref}.${nochr}-q${mq}.tsv 2> ${vcfpref}.${nochr}-q${mq}.stderr"
    done
done

if isfound run-reduce-chroms "${2}"; then
    python "${EVALROOT}/powerlaw/powerlaw-plot.py" \
            "${plotheader}" \
            $univ_offset \
            "${vcfpref}.powerlaw.pdf" \
            ${vcfpref}.*.tsv \
            > ${vcfpref}.tsv.xls
fi

