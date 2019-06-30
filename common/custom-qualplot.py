#!/usr/bin/env python
import gzip
import sys

from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1 sample2 ...
def main():
    title = sys.argv[1].split(';')[0].replace('\\n', '\n')
    vcf = sys.argv[2]
    bed = sys.argv[3]
    pdfout = sys.argv[4]
    type2qual2cnt = ([0] * 120, [0] * 120)
    type2percfa2cnt = ([0] * (100+1), [0] * (100+1))
    
    with gzip.open(vcf) as vcffile:
        for line in vcffile:
            if line.startswith('#'): continue
            tokens = line.split()
            if tokens[4].startswith('<'): continue
            tokens = line.split()
            qual = int(tokens[5])
            type = (0 if len(tokens[3]) == len(tokens[4]) else 1)
            assert qual < 130, '{} < 130 failed for QUAL'.format(qual)
            type2qual2cnt[type][qual] += 1
            tFA = 0
            tAD = 0
            tRAD = 0
            tDP = 0
            for subtok in tokens[7].split(';'):
                name_value_pair = subtok.split('=')
                if 2 == len(name_value_pair):
                    name, value = tuple(name_value_pair)
                    if name == 'tDP': 
                        tDP += int(value)
                    if name == 'tADR':
                        dtokens = [int(x) for x in value.split(',')]
                        tRD = int(dtokens[0])
                        tAD = int(dtokens[1])
                        tRAD = tRD + tAD
                        tFA = int(100 * float(tAD) / float(tRAD))
            assert (tAD <= tDP)
            if tDP >= 200: type2percfa2cnt[type][100 * tAD / tDP] += 1
            # if tRAD >= 60: type2percfa2cnt[type][tFA] += 1
    regionsize = 0
    with open(bed) as bedfile:
        for line in bedfile:
            tokens = line.split()
            regionsize += int(tokens[2]) - int(tokens[1])
    for qual in range(120): print('SNV\tQUAL\t{}\t{}\t{}'.format(qual, type2qual2cnt[0][qual], regionsize))
    for qual in range(120): print('INDEL\tQUAL\t{}\t{}\t{}'.format(qual, type2qual2cnt[1][qual], regionsize))
    for percfa in range(100): print('SNV\tFA100\t{}\t{}\t{}'.format(percfa, type2percfa2cnt[0][percfa], regionsize))
    for percfa in range(100): print('INDEL\tFA100\t{}\t{}\t{}'.format(percfa, type2percfa2cnt[1][percfa], regionsize))
    
    obsquals = [i for i in range(len(type2qual2cnt[0])) if 0 < max((type2qual2cnt[0][i], type2qual2cnt[1][i]))]
    minqual = min(obsquals)
    maxqual = max(obsquals)
    
    obspercfas = [i for i in range(len(type2percfa2cnt[0])) if 0 < max((type2percfa2cnt[0][i], type2percfa2cnt[1][i]))]
    minpercfa = min(obspercfas)
    maxpercfa = max(obspercfas) 
    
    with PdfPages(pdfout) as pdf:
        fig, ax1 = plt.subplots()
        qcs = ([qc for qc in enumerate(type2qual2cnt[0]) if qc[1] > 0])
        ax1.scatter([x[0] for x in qcs], [x[1] for x in qcs], label = 'Number of false positive SNV calls')
        qcs = ([qc for qc in enumerate(type2qual2cnt[1]) if qc[1] > 0])
        ax1.scatter([x[0] for x in qcs], [x[1] for x in qcs], label = 'Number of false positive InDel calls')
        xs = list(range(minqual, maxqual+1))
        ys = [(regionsize * 10**(-float(x)/10)) for x in xs]
        ax1.plot(xs, ys, label = 'Theoretical maxinum number of false positive calls')
        ax1.set_xlabel('UVC QUAL')
        ax1.set_ylabel('False-positive count')
        
        ax1.set_yscale('log')
        
        ax1.legend()
        ax1.set_title(title)
        
        plt.tight_layout()
        pdf.savefig()
        plt.close()

        plt.rc('text', usetex=False)
        
        fig, ax2 = plt.subplots()
        qcs = ([qc for qc in enumerate(type2percfa2cnt[0]) if qc[1] > 0])
        ax2.scatter([x[0] for x in qcs], [x[1] for x in qcs], label = 'Number of false positive SNV calls')
        qcs = ([qc for qc in enumerate(type2percfa2cnt[1]) if qc[1] > 0])
        ax2.scatter([x[0] for x in qcs], [x[1] for x in qcs], label = 'Number of false positive InDel calls')
        xs = list(range(max(minpercfa,1), maxpercfa+1, 1))
        ys = [(regionsize / float(x**3)) for x in xs]
        ax2.plot(xs, ys, label = 'Theoretical maximum number of false positive calls')
        ax2.set_xlabel('Variant allele fraction')
        ax2.set_ylabel('False-positive count')
        
        ax2.set_xscale('log')
        ax2.set_yscale('log') 
        
        ax2.legend()
        pdf.savefig()
        plt.close()
    
if __name__ == '__main__':
    main()
    
