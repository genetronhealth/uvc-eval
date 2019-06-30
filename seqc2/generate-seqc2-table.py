#!/usr/bin/env python
import sys

from collections import defaultdict

def fanme_to_sra_to_assay_dict(fname):
    ret = {}
    with open(fname) as metafile:
        for line in metafile:
            tokens = line.split(',')
            sra = tokens[0]
            assay = tokens[2]
            ret[sra] = assay
    return ret
    
def vc_fscore_to_priority(vc_fscore):
    vc, fscore = vc_fscore
    if 'UVC' in vc:
        return '0' + vc
    if 'Mutect2' in vc or 'Strelka2' in vc or 'Varscan2' in vc:
        return '1' + vc
    return '2' + vc
    
def main():
    tnsamples_variantcaller_fscores = defaultdict(dict)
    
    sra_to_assay_dict = fanme_to_sra_to_assay_dict(sys.argv[1])
    
    for fpath in sys.argv[2:]:
        fname = fpath.split('/')[-1]
        tsample = fname.split('_')[0]
        nsample = fname.split('_')[1]
        
        with open(fpath) as file:
            for line in file:
                if line.startswith('#'): continue
                tokens = line.split()
                variantcaller = tokens[2]
                fscore = tokens[0].split('=')[-1]
                if fscore == 'Not-applicable': variantcaller = tokens[1]
                fscore = (float(fscore) if fscore != 'Not-applicable' else -1)
                tnsamples_variantcaller_fscores[(tsample, nsample)][variantcaller] = fscore
    
    (tsample, nsample), variantcaller_fscores = sorted(tnsamples_variantcaller_fscores.items())[0]
    sys.stdout.write('{} \t& {} '.format('Tumor/normal accessions', 'assay-type'))
    for variantcaller, fscore in sorted(variantcaller_fscores.items(), key=vc_fscore_to_priority):
        if variantcaller == 'SomaticSniper': continue
        sys.stdout.write(' \t& {}'.format(variantcaller))
    sys.stdout.write('\\\\\n')
    
    for (tsample, nsample), variantcaller_fscores in sorted(tnsamples_variantcaller_fscores.items()):
        # sra_to_assay_dict[nsample]
        sys.stdout.write('{}/{} \t& {} '.format(tsample, nsample, sra_to_assay_dict[tsample]))
        fscore_max = max(variantcaller_fscores.values())
        for variantcaller, fscore in sorted(variantcaller_fscores.items(), key=lambda x:vc_fscore_to_priority(x)):
            if variantcaller == 'SomaticSniper': continue
            fscore_out = ('N/A' if fscore <= -1 else ('{:.4f}'.format(fscore) if (fscore < fscore_max) else '\\textbf{{{:.4f}}}'.format(fscore)))
            sys.stdout.write(' \t& {}'.format(fscore_out))
        sys.stdout.write('\\\\\n')
    
if __name__ == '__main__':
    main()

