#!/usr/bin/env python
import sys

from collections import defaultdict

def fname_to_sample_to_depth_purity_dict(fname):
    ret = {}
    with open(fname) as metafile:
        for line in metafile:
            if line.startswith('#'): continue
            tokens = line.split()
            scenario = tokens[0]
            #if int(scenario) > 32: continue
            assay = tokens[1]
            tpurity = tokens[7]
            npurity = tokens[8]
            tdepth = tokens[9]
            ndepth = tokens[10]
            tsample = 'HGM12PA{}T'.format(scenario) # HGM12PA01T
            nsample = 'HGM12PA{}N'.format(scenario)
            ret[tsample] = (float(tdepth), float(tpurity))
            ret[nsample] = (float(ndepth), float(npurity))
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
    
    sample_to_depth_purity_dict = fname_to_sample_to_depth_purity_dict(sys.argv[1])
    sys.stderr.write('{}\n'.format(sample_to_depth_purity_dict))
    for fpath in sys.argv[2:]:
        fname = fpath.split('/')[-1]
        tsample = fname.split('_')[0]
        nsample = fname.split('_')[1]
        
        with open(fpath) as file:
            for line in file:
                if line.startswith('#'): continue
                tokens = line.split()
                fscore = tokens[0].split('=')[-1]
                if fscore == 'Not-applicable':
                    variantcaller = tokens[1]
                    fscore = -1
                else:
                    variantcaller = tokens[2]
                    fscore = float(fscore)
                tnsamples_variantcaller_fscores[(tsample, nsample)][variantcaller] = fscore
    
    (tsample, nsample), variantcaller_fscores = sorted(tnsamples_variantcaller_fscores.items())[0]
    print('\\toprule')
    sys.stdout.write('{} \t& {} '.format('T/N depths', 'T/N purities'))
    vcs = []
    for variantcaller, fscore in sorted(variantcaller_fscores.items(), key=vc_fscore_to_priority):
        if 'SomaticSniper' == variantcaller: continue
        sys.stdout.write(' \t& {}'.format(variantcaller))
        vcs.append(variantcaller)
    sys.stdout.write('\\\\\n')
    print('\\midrule')
    for (tsample, nsample), variantcaller_fscores in sorted(tnsamples_variantcaller_fscores.items()):
        tdepth, tpurity = sample_to_depth_purity_dict[tsample]
        ndepth, npurity = sample_to_depth_purity_dict[nsample]
        sys.stdout.write('{}/{} \t& {:.4f}/{:.4f} '.format(tdepth, ndepth, float(tpurity), float(npurity)))
        fscore_max = max(variantcaller_fscores.values())
        for variantcaller in vcs:
            if 'SomaticSniper' == variantcaller: continue
            # for variantcaller, fscore in sorted(variantcaller_fscores.items(), key=lambda x:vc_fscore_to_priority(x)):
            fscore = variantcaller_fscores.get(variantcaller, -2)
            fscore_out = ('N/A' if fscore <= -1 else ('{:.4f}'.format(fscore) if (fscore < fscore_max) else '\\textbf{{{:.4f}}}'.format(fscore)))
            sys.stdout.write(' \t& {}'.format(fscore_out))
        sys.stdout.write('\\\\\n')
    print('\\bottomrule')

if __name__ == '__main__':
    main()

