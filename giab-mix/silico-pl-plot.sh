#!/usr/bin/env sh

SCRIPTDIR=$(dirname $(which $0))
source "${SCRIPTDIR}/../source-eval.sh"

evalflag=" --squash-ploidy --ref-overlap "

if [ $(echo "${PAT}" | grep -c run-HG001) -gt 0 ] ; then
    hg00x=HG001
    hg00xVCF="${HG001VCF}"
    hg00xBED="${HG001BED}"
fi
if [ $(echo "${PAT}" | grep -c run-HG002) -gt 0  ] ; then
    hg00x=HG002
    hg00xVCF="${HG002VCF}"
    hg00xBED="${HG002BED}"
fi

for nochr in 22; do
    #nochrlen=$(cat ${REF}.dict | grep -P "\tSN:${nochr}\t" | awk '{print $3}' | awk -F":" '{print $2}')
    
    resDIR="${ROOTDIR}/eval/tn/insilico-pl/${Tparamset}/"
    mkdir -p "${resDIR}"
    nochrBED="${resDIR}/${hg00x}.hs37d5.300x.${nochr}.bed"
    nochrISEC="${resDIR}/${hg00x}.hs37d5.300x.${nochr}.isecdir" 
    truthVCF="${resDIR}/${hg00x}.hs37d5.300x.${nochr}.truth.vcf.gz"
    uvcresVCF="${resDIR}/${hg00x}.hs37d5.300x.${nochr}.uvc1.vcf.gz"
    vcallVCF="${resDIR}/${hg00x}.hs37d5.300x.${nochr}.vcall.vcf.gz"
    fpVCF="${resDIR}/${hg00x}.hs37d5.300x.${nochr}.isecdir/0001.vcf.gz"
    mindepth=200
    cat "${hg00xBED}" | grep -P "^${nochr}\t" > "${nochrBED}"
    regsize=$(cat "${nochrBED}" | awk '{regsize += $3-$2} END {print regsize}')
    
if [ $(echo ${PAT} | grep -c "skip-result-generation") -eq 0 ]; then
    
    bcftools view "${hg00xVCF}" -r ${nochr} | bcftools norm ${GROUNDTRUTH_NORM_FLAGS} -Oz -o "${truthVCF}" -
    if [ $(echo $1 | grep -c "disable-1-end") -eq 0 ]; then
        samtools view -bh "${ROOTDIR}/eval/tn/insilico-rawinput/${hg00x}.hs37d5.300x.bam" 22 > "${ROOTDIR}/eval/tn/insilico-pl/${hg00x}.hs37d5.300x.22.bam"
        samtools index -@ ${ncpus} "${ROOTDIR}/eval/tn/insilico-pl/${hg00x}.hs37d5.300x.${nochr}.bam"
    fi
    if [ $(echo $1 | grep -c "disable-2-end") -eq 0 ]; then
        "${UVC}" -t ${ncpus} -f "${REF}" -o "${uvcresVCF}" "${ROOTDIR}/eval/tn/insilico-pl/${hg00x}.hs37d5.300x.${nochr}.bam"
        bcftools index -ft "${uvcresVCF}"
    fi
    
    bcftools view "${uvcresVCF}" -r ${nochr} | bcftools norm ${GROUNDTRUTH_NORM_FLAGS} -Oz -o "${vcallVCF}" -
    bcftools index -ft "${vcallVCF}"
    ## ${GROUNDTRUTH_ISEC_FLAGS} is not applied here 
fi
if [ $(echo ${PAT} | grep -c "skip-result-isec") -eq 0 ]; then
    bcftools index -ft "${truthVCF}"
    bcftools isec --threads 4 -c all -Oz -T "${nochrBED}" -p "${nochrISEC}" "${truthVCF}" "${vcallVCF}"
    bcftools index -ft "${fpVCF}"
    
    rm -r "${nochrISEC}/uvc1.vcfeval.outdir/" || true
    ${vcfeval1} -f QUAL \
            -b "${truthVCF}" \
            -c "${vcallVCF}" \
            -e "${nochrBED}" \
            -t "${HS37D5}" \
            -o "${nochrISEC}/uvc1.vcfeval.outdir/" || true
fi
    version=$(bcftools view -h ${fpVCF} | grep "^##variantCallerVersion=" | awk -F"=" '{print $2}' | awk '{print $1}' | awk -F"." '{print $NF}')
    for vartype in snps indels; do
        if [ $(echo ${vartype} | grep indels | wc -l) -gt 0 ]; then
            vartype_string=InDels
        else
            vartype_string=SNVs
        fi
        
        bcftools view "${fpVCF}" -v ${vartype} -i "bDP>=${mindepth}" | bcftools query -f " [ %bADR ] [ %bDP ] \n" | sed 's/,/ /g' | awk '{print int($2*100 / $3)}' | sort -n | uniq -c | awk '{print $2, $1}' > "${ROOTDIR}/eval/tn/insilico-pl/${hg00x}.hs37d5.300x.${nochr}.atleast${mindepth}DP_100xFA_${vartype}.tsv"
        # bcftools view "${fpVCF}" -v ${vartype} -i "bDP>=${mindepth}" | bcftools query -f " [%bFA] \n" | awk '{print int($1*100)}' | sort -n | uniq -c | awk '{print $2, $1}' > "${ROOTDIR}/eval/tn/insilico-pl/${hg00x}.hs37d5.300x.${nochr}.atleast${mindepth}DP_100xFA_${vartype}.tsv"
            c=$(cat <<EOF
import numpy as np
import scipy
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def powlaw(x, a, b) :
    return a * np.power(x, b)
def linlaw(x, a, b) :
    return a + x * b

def curve_fit_log(xdata, ydata) :
    """Fit data to a power law with weights according to a log scale"""
    # Weights according to a log scale
    # Apply fscalex
    xdata_log = np.log10(xdata)
    # Apply fscaley
    ydata_log = np.log10(ydata)
    # Fit linear
    popt_log, pcov_log = curve_fit(linlaw, xdata_log, ydata_log)
    #print(popt_log, pcov_log)
    # Apply fscaley^-1 to fitted data
    ydatafit_log = np.power(10, linlaw(xdata_log, *popt_log))
    # There is no need to apply fscalex^-1 as original data is already available
    return (popt_log, pcov_log, ydatafit_log)

def func_powerlaw(x, m, c, c0):
    return 0 + x**m * c

target_func = powlaw
vafs1 = []
cnts1 = []
with open("${ROOTDIR}/eval/tn/insilico-pl/${hg00x}.hs37d5.300x.${nochr}.atleast${mindepth}DP_100xFA_${vartype}.tsv") as f:
    for line in f:
        vafs1.append(float(line.split()[0]))
        cnts1.append(float(line.split()[1]))
        vafs2 = [v[0] for v in zip(vafs1, cnts1) if (v[0] > 0)]
        cnts2 = [v[1] for v in zip(vafs1, cnts1) if (v[0] > 0)]
        vafs3 = [v[0] for v in zip(vafs2, cnts2) if (v[0] > 0 and v[0] < 15)]
        cnts3 = [v[1] for v in zip(vafs2, cnts2) if (v[0] > 0 and v[0] < 15)]
        
        #vafs = np.array(vafs)
        #cnts = np.array(cnts)
        #popt, pcov = curve_fit(target_func, vafs3, cnts3, p0 = np.asarray([-10,100]))
        #popt_log, pcov_log, ydatafit_log = curve_fit_log(vafs3, cnts3)
        #a=scipy.stats.powerlaw.fit(3)
expvafs = [i for i in range(1, 100, 1)]
expcnts = [(${regsize} / 1e9 * (100.0 / vaf)**(3)) for vaf in expvafs]
plt.suptitle("${hg00x} chromosome=${nochr} min_depth=${mindepth} region_size=${regsize}")
plt.title("FittedEquation: \$FalsePositiveRate=AlleleFractionPercent^{-3}\$")
ret = plt.scatter(vafs2, cnts2)
#plt.plot(vafs2, target_func(vafs2, *popt), '--')
plt.plot(expvafs, expcnts, '--')
#plt.plot(ydatafit_log)
plt.xlabel("Binned variant allele fraction in percentage for ${vartype_string}")
plt.ylabel("Number of false positive variants")
plt.yscale('log')
plt.xscale('log')
plt.grid()
plt.savefig("${ROOTDIR}/eval/tn/insilico-pl/${hg00x}_hs37d5_300x_${nochr}_atleast${mindepth}DP_100xFA_${vartype}_${version}.pdf")
print(zip(vafs1, cnts1))

EOF)
            
        echo "${c}" > "${nochrISEC}/plotcmd.py"
        python -c "$c"
    done
done

