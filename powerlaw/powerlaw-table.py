#!/usr/bin/env python
import gzip
import sys

# stdin: call-vcf
# stdout: tsv file with number of positions satisfying 200 DP threshold \n distribution of FP allele frac perc
# param1: truth-vcf-gz

#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1 sample2 ...
def main():
    
    minMQ = (int(sys.argv[3]) if (len(sys.argv) > 3) else 0)
    nbins = (int(sys.argv[4]) if (len(sys.argv) > 4) else 100)
    mindp = (int(sys.argv[5]) if (len(sys.argv) > 5) else 200)
    
    # condMQ = ('KeptReadsWithMQ >= {}'.format(sys.argv[2]) if (len(sys.argv) > 2) else 'KeptReadsWithMQ >= 0')
    sys.stderr.write('Start reading from the file {}\n'.format(sys.argv[2]))
    chr_pos_truth_set = set([])
    with gzip.open(sys.argv[2]) as vcffile:
        for lineno, line in enumerate(vcffile):
            if line.startswith('#'): continue
            tokens = line.split()
            if tokens[4].startswith('<'): continue
            tokens = line.split()
            chrom, pos, _, ref, alts = tuple(tokens[0:5])
            if chrom != sys.argv[1] and sys.argv[1] not in ['.', '*']: continue
            #max_indel_len = 0
            #for alt in alts.split(','):
            #    if not alt.startswith('<'):
            #        max_indel_len = max([max_indel_len, abs(len(alt) - len(ref))])
            pos_beg = int(pos) - 8 # max([max_indel_len, 8])
            pos_end = int(pos) + 8 # max([max_indel_len, 8])
            for varpos in range(pos_beg, pos_end + 1, 1):
                chr_pos_truth_set.add((chrom, varpos))
            if lineno % (500*1000) == 0: sys.stderr.write('Read {} lines from the file {}\n'.format(lineno, sys.argv[2]))
    
    sys.stderr.write('Start reading from stdin\n')
    num_positions = 0
    snv_fa100_to_count = [0] * (nbins + 1)
    non_snv_fa100_to_count = [0] * (nbins + 1)
    for lineno, line in enumerate(sys.stdin):
        if line.startswith('#'): continue
        tokens = line.split()
        chr_pos = (tokens[0], int(tokens[1]))
        fmtkeys = tokens[8].split(':')
        fmtvals = tokens[9].split(':')
        
        infodp = 0
        for kv in tokens[7].split(';'):
            subtokens = kv.split('=')
            if subtokens[0] == 'DP': infodp = int(subtokens[1]) # int(v.split(',')[0])
        fmtdp = 0
        for k, v in zip(fmtkeys, fmtvals):
            if k == 'DP': fmtdp = int(v) # int(v.split(',')[0])
        
        dp = fmtdp
        
        if dp < mindp: continue
        num_positions += 1
        if tokens[4] == '<*>' or chr_pos in chr_pos_truth_set: continue
        fa100_snv_max = -1
        fa100_non_snv_max = -1
        for k, v in zip(fmtkeys, fmtvals):
            if k == 'AD':
                for alt, ad in zip([tokens[3]] + tokens[4].split(','), v.split(',')):
                    ad = int(ad)
                    assert ad <= dp, '{} <= {} failed for line {}'.format(ad, dp, line)
                    if alt.startswith('<') or alt == '*' or alt == tokens[3]: continue
                    fa100 = nbins * ad / dp;
                    #sys.stderr.write('The following line is FP: {}\n'.format(line))
                    if 1 == len(alt): fa100_snv_max = max([fa100_snv_max, fa100])
                    else: fa100_non_snv_max = max([fa100_non_snv_max, fa100])
        if fa100_snv_max > 0: snv_fa100_to_count[fa100_snv_max] += 1
        if fa100_non_snv_max > 0: non_snv_fa100_to_count[fa100_non_snv_max] += 1
    print('##KeptReadsWithMQ >= {}'.format(minMQ))
    print('##num_positions = {}'.format(num_positions))
    print('#variant-type\tvariant-allele-fraction\tfalse-positive-allele-count')
    for i in range(nbins + 1):
        print('SNV\t{}\t{}'.format(float(i) * 100 / nbins, snv_fa100_to_count[i]))
        #print('SNV\t{}\t{}\t{}'.format(i, snv_fa100_to_count[i], condMQ))
    for i in range(nbins + 1):
        print('NON_SNV\t{}\t{}'.format(float(i) * 100 / nbins, non_snv_fa100_to_count[i]))
        #print('NON_SNV\t{}\t{}\t{}'.format(i, non_snv_fa100_to_count[i], condMQ))
    
if __name__ == '__main__':
    main()
    
