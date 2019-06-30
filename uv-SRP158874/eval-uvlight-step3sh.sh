SCRIPTDIR=$(dirname $(which $0))
if [ -z "${EVALROOT}" ]; then
    EVALROOT=$(cd "${SCRIPTDIR}/../" && pwd)
fi

REF="${EVALROOT}/datafiles/hs37d5.fa"

if [ $(echo $1 | grep -c "disable-1-end") -eq 0 ]; then
    for d in SRR7757437_SRR7757438_NT SRR7757439_SRR7757440_NT SRR7757441_SRR7757442_NT SRR7757443_SRR7757444_NT; do 
        bcftools view ${EVALROOT}/uv-SRP158874/${d}/${d}*_N_uvc1.vcf.gz -i "QUAL>0" 3 \
            | python ${EVALROOT}/uv-SRP158874/eval-uvlight-step3plot.py ${EVALROOT}/uv-SRP158874/${d}_plot.pdf 1>/dev/null ; 
    done
fi

if [ $(echo $1 | grep -c "disable-2-end") -eq 0 ]; then
    for i in 0 1; do 
        for sra in SRR7757437 SRR7757439 SRR7757441 SRR7757443; do 
            for t in  N; do 
                f=$(ls ${EVALROOT}/uv-SRP158874/${sra}*_NT/${sra}*_NT_${t}_uvc1.vcf.gz); 
                s=$(bcftools view $f -h | tail -n1 |awk '{print $(NF-1)}')
                bcftools view -h "${f}" | grep "^##variantCallerVersion="
                line1=$(echo $f:RPLA $(bcftools view $f -v snps -i "QUAL>0 && ALT=='A'" 19:49990694-49990694 \
                    | bcftools query -s ${s} -f " %tADCR{1} %tADCR{0} [%cDP2f{1}] [%cDP2r{1}] [%cDP2f{0}] [%cDP2r{0}]") | awk '{print $2,"/",$2+$3, " & ", $4+$5,"/",$4+$5+$6+$7}');
                line2=$(echo $f:DPH3 $(bcftools view $f -v snps -i "QUAL>0 && ALT=='T'"  3:16306504-16306504 \
                    | bcftools query -s ${s} -f " %tADCR{1} %tADCR{0} [%cDP2f{1}] [%cDP2r{1}] [%cDP2f{0}] [%cDP2r{0}]") | awk '{print $2,"/",$2+$3, " & ", $4+$5,"/",$4+$5+$6+$7}');
                echo "RPLA//DPH3 ${line1} & ${line2}"
                #echo $f:DPH3 $(
                #    bcftools view $f -v snps -i "QUAL>0 && ALT=='T'" 3:16306504-16306504 | 
                #    bcftools query -s ${s} -f "[%cRDTT{${i}}] [%cRDT1{${i}}] [%cRDTN{${i}}] [%cADTT{${i}}] [%cADT1{$i}] [%cADTN{$i}] " 
                #) | awk '{print $1, " & ", $5-$6-$7 " / " $2-$3-$4}' ; 
                #echo $f:RPLA $(
                #    bcftools view $f -v snps -i "QUAL>0 && ALT=='A'" 19:49990694-49990694 | 
                #    bcftools query -s ${s} -f "[%cRDTT{${i}}] [%cRDT1{${i}}] [%cRDTN{${i}}] [%cADTT{${i}}] [%cADT1{$i}] [%cADTN{$i}] " 
                #) | awk '{print $1, " & ", $5-$6-$7 " / " $2-$3-$4}' ; 
            done ; 
        done ; 
    done
fi

if [ $(echo $1 | grep -c "disable-3-end") -eq 0 ]; then
    for srr in SRR7757437 SRR7757439 SRR7757441 SRR7757443; do
        echo $srr with original T-N; 
        bcftools view ${EVALROOT}/uv-SRP158874/${srr}_*_NT/${srr}_*_NT_N_uvc1.vcf.gz -i "SomaticQ>25" -v snps -H ; 
        echo $srr with reversed T-N; 
        bcftools view ${EVALROOT}/uv-SRP158874/${srr}_*_TN/${srr}_*_TN_N_uvc1.vcf.gz -i "SomaticQ>25" -v snps -H ; 
    done
fi

