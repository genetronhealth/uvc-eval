import gzip
import math
import multiprocessing
import os
import sys

def vcfpath2ss2vq(path):
    ss2vq = {}
    sample = ''
    with gzip.open(path) as vcf:
        for line in vcf:
            if line.startswith('##'): continue
            tokens = line.rstrip().split('\t')
            if line.startswith('#CHROM\t'):
                sample = tokens[-1]
                continue
            #print('processing token {}'.format(tokens))
            chrom = tokens[0]
            pos = tokens[1]
            ref = tokens[3]
            alt = tokens[4]
            if alt.startswith('<'): continue
            qual = float(tokens[5])
            ss = (sample, chrom, pos, ref, alt)
            assert ss not in ss2vq
            ss2vq[ss] = qual
    return (sample, ss2vq)
    
def main():
    header = ''
    samples = set([])
    sample2tpcount = {}
    ss2vq = {} # sample-site to variant quality
    inlines = [line for line in sys.stdin]
    for line in inlines:
        if line.startswith('Comment\t'): header = 'Original-order\tis-true-pos\tartifact-score\tUVC-QUAL\t' + line.rstrip()
        if line.startswith('Comment\t') or line.strip() == '': continue
        tokens = line.rstrip().split('\t')
        #print('processing token {}'.format(tokens))
        sample = tokens[1]
        chrom = tokens[2]
        pos = tokens[3]
        ref = tokens[4]
        alt = tokens[5]
        ss = (sample, chrom, pos, ref, alt)
        ss2vq[ss] = -1000*1000
        sample2tpcount[sample] = sample2tpcount.get(sample, 0) + 1
    
    s_ss2vq_all = multiprocessing.Pool(processes=8).map(vcfpath2ss2vq, sys.argv[1:])
    for s_ss2vq_each in s_ss2vq_all:
        sample = s_ss2vq_each[0]
        samples.add(sample)
        for (ss, vq) in s_ss2vq_each[1].items():
            ss2vq[ss] = vq
    
    is_false_pos = is_true_pos = is_liketrue_pos = 0
    print(header)
    for i, line in enumerate(inlines):
        if line.startswith('Comment\t') or line.strip() == '': continue
        tokens = line.split('\t')
        sample = tokens[1]
        chrom = tokens[2]
        pos = tokens[3]
        ref = tokens[4]
        alt = tokens[5]
        ss = (sample, chrom, pos, ref, alt)

        cmt = tokens[0].lower().strip()
        if 'inspected' != cmt:
            is_false_pos = ('not reported' in cmt)
            is_true_pos = (cmt.startswith('reported'))
        is_liketrue_pos = (('BQ too low' in cmt) and is_true_pos)
        if (is_true_pos + is_false_pos != 1):
            sys.stderr.write('Skipping the line {} because it has no manual-review result!\n'.format(line))
            is_liketrue_pos = 0
            continue 
        category = tokens[-1].lower().strip()
        if category.startswith('polymorphism'): 
            artifact_score = -0.5
        elif category.startswith('artifact') or category.startswith('probably true'):
            artifact_score = float(tokens[-1].strip().split()[-1].lstrip('(').rstrip(')'))
        else:
            sys.stderr.write('WARNING: the line ({}) is assumed to have the previous artifact score because it has a malformed category (artifact, probably true, polymorphism, etc.)!\n'.format(line))
        if ss2vq[ss] >= 0:
            print('{}\t{}\t{}\t{}\t{}'.format(i+1, (0 if is_false_pos else (2)), artifact_score, ss2vq[ss], line.rstrip()))
        elif sample in samples:
            dummyvq = 0 # -(10/math.log(10)) * math.log(1 + sample2tpcount.get(sample, 0))
            print('{}\t{}\t{}\t{}\t{}'.format(i+1, (0 if is_false_pos else (2)), artifact_score, dummyvq, line.rstrip()))

if __name__ == '__main__':
    main()

