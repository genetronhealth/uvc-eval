#!/usr/bin/env python

import sys
from collections import defaultdict

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
#plt.rcParams["figure.figsize"] = (10,6) # width height

from matplotlib.backends.backend_pdf import PdfPages

MARKER_SEQ = 'sdo^<>vPxp*'

# header pdf tsv1 [tsv2, ...]
def main():
    mq2npos = defaultdict(int)
    snv_mq2percfa2cnt = {0 : defaultdict(int), 60: defaultdict(int)}
    nonsnv_mq2percfa2cnt = {0 : defaultdict(int), 60: defaultdict(int)}
    
    # MQ > pos > SNV-NON-SNV > 
    plotheader =  sys.argv[1]
    univ_offset = sys.argv[2]
    for tsv in sys.argv[4:]:
        mq = -1
        with open(tsv) as tsvfile:
            for lineno, line in enumerate(tsvfile):
                tokens = line.split()
                if line.startswith('##KeptReadsWithMQ'): mq = int(tokens[2])
                elif line.startswith('##num_positions'): mq2npos[mq] += int(tokens[2])
                elif tokens[0] == 'SNV':
                    snv_mq2percfa2cnt[mq][float(tokens[1])] += int(tokens[2])
                elif tokens[0] == 'NON_SNV':
                    nonsnv_mq2percfa2cnt[mq][float(tokens[1])] += int(tokens[2])
    
    print('{}\t{}\t{}\t{}'.format('KeptReadsWithMQ', 'SNV', 'percentVAF', 'count'))
    for mq in [0, 60]:
        for percfa, count in sorted(snv_mq2percfa2cnt[mq].items()):
            print('{}\t{}\t{}\t{}'.format(mq, 'SNV', percfa, count))
        for percfa, count in sorted(nonsnv_mq2percfa2cnt[mq].items()):
            print('{}\t{}\t{}\t{}'.format(mq, 'NON_SNV', percfa, count))
    
    with PdfPages(sys.argv[3]) as pdf:
        #fig = plt.figure(figsize=(6, 10))
        #ax1 = fig
        fig, ax1 = plt.subplots()
        fig.set_figheight(6)
        fig.set_figwidth(10)
        
        #thetitle = '''Cell_line=HG001 platform=Illumina/HiSeq average_sequencing_depth=300x min_depth_required=200x''' #.format(regionsize)
        #thetitle = plotheader
        plt.title(plotheader + '\n' + 'FittedEquation: $FalsePositiveRate=(10 \\times AlleleFractionPercent / c)^{{-3}} / 3$ where $c = {}$'.format(univ_offset))
        #title_eq = "FittedEquation: $FalsePositiveRate=(10 \\times AlleleFractionPercent)^{-3} / 3$"
        #plt.title(thetitle + '\n' + title_eq)
        idx = 0
        mqlist = [0] # [0, 60]
        for mq in mqlist:
            nbins = 0
            regionsize = mq2npos[mq]
            for vartype, mq2percfa2cnt in zip(['SNV', 'NON_SNV'], [snv_mq2percfa2cnt, nonsnv_mq2percfa2cnt]):
                fa100list, cntlist = zip(*sorted(mq2percfa2cnt[mq].items()))
                nbins = len(fa100list) - 1;
                pc = 1.0 * 100 / float(len(fa100list) - 1)
                ymult = 100.0 / float(len(fa100list) - 1)
                if len(mqlist) == 1:
                    label = '{} NumberOfPositionsCalled={}'.format(vartype, mq2npos[mq])
                else:
                    label = '{} KeptReadWithMQ>={} NumberOfPositionsCalled={}'.format(vartype, mq, mq2npos[mq])
                ret = ax1.scatter(
                    [(x + pc * 0.5) for x in fa100list],
                    cntlist,
                    label = label,
                    marker = MARKER_SEQ[idx],
                    alpha = 0.5)
                idx += 1
            sys.stderr.write('pc = {}, ymult = {}\n'.format(pc, ymult))
            xs = [((float(x) + pc)) for x in range(0,  100)]
            ys = [(regionsize / float((10 * x / float(univ_offset))**3) * ymult / 3.0) for x in xs]
            ax1.plot(xs, ys, label = 'Theoretical SNV with {} bins'.format(nbins))
        ax1.set_xlim((pc / 2.0, 100.0))
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlabel('Binned allele fraction in percentage')
        ax1.set_ylabel('Number of false positive variants')
        ax1.grid()
        
        plt.legend(loc = 'upper right', fontsize = plt.rcParams['font.size'] - 1)
        plt.tight_layout()
        pdf.savefig()
        plt.close()
        print('The file {} is closed'.format(sys.argv[3]))
        
if __name__ == '__main__':
    main()
    
