# -*- coding: UTF-8 -*-

import codecs
import os
import sys
import gzip
#import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset, zoomed_inset_axes)

from matplotlib.lines import Line2D
from pylab import mpl

def create_dummy_line(**kwds):
    return Line2D([], [], **kwds)

plt.rcParams['font.family'] = ['sans-serif']
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 8

#print(plt.rcParams)

colorcycles = plt.rcParams['axes.prop_cycle'].by_key()['color']

outfname = sys.argv[1]
title = sys.argv[2]
tokens = sys.argv[3].split(',')
if len(tokens) == 1:
    xlower = ylower = float(sys.argv[3])
else:
    xlower = float(tokens[0])
    ylower = float(tokens[1])


getbest = (1 < len(title.split(';')) and 'isbest' in title.split(';')[1])
should_assert_baseline_var_num = (1 < len(title.split(';')) and 'disable-baseline-var-num-assertion' not in title.split(';')[1])

def get_label(fname1, getbest = 0):
    fname = fname1.lower()
    if 'invalid-name' in fname:
        return (0, fname1, 'invalid-name')
    elif   'tnscope'      in fname:
        if 'filtered' in fname:
            return ('' if getbest else 'TNscope-PASS')
        else:
            return 'TNscope'
    elif 'tnhaplotyper2' in fname:
        return 'TNhaplotyper2'
    elif 'strelka2' in fname:
        if 'disableEVS'.lower() in fname:
            if 'filtered' in fname:
                return (2, fname1, ('' if getbest else 'Strelka2-diasbleEVS-PASS'))
            else:
                return (2, fname1, ('' if getbest else 'Strelka2-diasbleEVS'))
        else:
            return (2, fname1, 'Strelka2')
    elif 'smcounter2.vcfeval' in fname:
        return (2, fname1, 'smCounter2')
    elif 'mutect2' in fname:
        return (2, fname1, ('Mutect2-all' if 'vcfeval-all' in fname else 'Mutect2-filter'))
    elif 'varscan2' in fname:
        return (2, fname1, ('Varscan2-all' if 'vcfeval-all' in fname else 'Varscan2-filter'))
    elif 'somaticsniper' in fname:
        return (3, fname1, 'SomaticSniper')
    elif 'lolopicker' in fname:
        return (99, fname1, 'LoLoPicker')
    elif 'lofreq' in fname:
        return (3, fname1, 'LoFreq*')
    
    # germline callers
    elif 'gatk4hc' in fname:
        return (1.5, fname1, 'GATK4HaplotypeCaller')
    elif 'freebayes' in fname:
        return (2.5, fname1, 'freebayes')
    elif 'bcftools' in fname:
        return (3.5, fname1, 'bcftools')
    elif 'lancet' in fname:
        return (4.5, fname1, ('Lancet-all' if 'vcfeval-all' in fname else 'Lancet-filter'))
    elif 'octopus' in fname:
        return (4.6, fname1, ('Octopus-all' if 'vcfeval-all' in fname else 'Octopus-filter'))
    elif 'uvc' in fname:
        fpvcf = fname1.replace('non_snp_roc.tsv.gz', '').replace('snp_roc.tsv.gz', '') + 'fp.vcf.gz'
        uvcversion = 'NoVersion'
        if os.path.exists(fpvcf):
            uvcheader = os.popen('bcftools view --header-only {}'.format(fpvcf)).read()
            for line in uvcheader.split('\n'):
                if line.startswith('##variantCallerVersion='): uvcversion = line.rstrip().split('.')[-1].split('-')[-1][0:7]
        if 'vcfeval-all' in fname:
            return (1, fname1, 'UVC-version-' + uvcversion)
        elif 'vcfeval-tFA5perc.outdir'.lower() in fname:
            return (1, fname1, ('UVC' if getbest else 'UVC-5%-VAF'))
        elif 'flt.vcfeval'.lower() in fname or 'vcfeval-filter'.lower() in fname:
            return (1, fname1, ('UVC' if getbest else 'UVC-version-' + uvcversion))
        elif 'vcfeval.outdir' in fname:
            return (1, fname1, ('' if getbest else 'UVC'))
        else:
            return (1, fname1, 'UVC-HF-' + fname.split('/')[-2].replace('uvc', '').replace('vcfeval', '').replace('outdir', ''))
    else:
        raise RuntimeException('Invalid fname {}'.format(fname))
    
LINESTYLES = ['solid', 'dashed', 'dashdot', 'dotted']
#LINEWIDTHS = [3, 1.5, 0.5]
CC_LEN = len(colorcycles)
MARKER_SEQ = 'sdo^<>vPxp*'
MARKER_SEQ_LEN = len(MARKER_SEQ)

s1=1
s2=2
s3=4
style1= (0, (s3, s1,  s2, s1,  s1, s1))
style2= (0, (s3, s1,  s1, s1,  s2, s1))
style3= (0, (s3, s2,  s2, s1,  s1, s1))
style4= (0, (s3, s2,  s1, s1,  s2, s1))
style5= (0, (s3, s3,  s2, s1,  s1, s1))
style6= (0, (s3, s3,  s1, s1,  s2, s1))

style1 = (0, (s3, s1, s1, s1))
style2 = (0, (s3, s2, s1, s1))
style3 = (0, (s3, s1, s2, s1))
style4 = (0, (s3, s2, s2, s1))
style5 = (0, (s3, s1, s1, s2))
style6 = (0, (s3, s2, s1, s2))
style7 = (0, (s3, s1, s2, s2))
style8 = (0, (s3, s2, s2, s2))

linestyles = ['solid', 'dashed', 'dashdot', 'dotted', style1, style2, style3, style4, style5, style6, style7, style8]
LS_LEN = len(linestyles)

linetypes = [
     {'color': colorcycles[ 0 % CC_LEN], 'linestyle': 'solid',   'marker': MARKER_SEQ[ 0 % MARKER_SEQ_LEN]},
     {'color': colorcycles[ 1 % CC_LEN], 'linestyle': 'dashed',  'marker': MARKER_SEQ[ 1 % MARKER_SEQ_LEN]},
     {'color': colorcycles[ 2 % CC_LEN], 'linestyle': 'dashdot', 'marker': MARKER_SEQ[ 2 % MARKER_SEQ_LEN]},
     {'color': colorcycles[ 3 % CC_LEN], 'linestyle': 'dotted',  'marker': MARKER_SEQ[ 3 % MARKER_SEQ_LEN]},
     {'color': colorcycles[ 4 % CC_LEN], 'linestyle': style1,    'marker': MARKER_SEQ[ 4 % MARKER_SEQ_LEN]},
     {'color': colorcycles[ 5 % CC_LEN], 'linestyle': style2,    'marker': MARKER_SEQ[ 5 % MARKER_SEQ_LEN]},
     {'color': colorcycles[ 6 % CC_LEN], 'linestyle': style3,    'marker': MARKER_SEQ[ 6 % MARKER_SEQ_LEN]},
     {'color': colorcycles[ 7 % CC_LEN], 'linestyle': style4,    'marker': MARKER_SEQ[ 7 % MARKER_SEQ_LEN]},
     {'color': colorcycles[ 8 % CC_LEN], 'linestyle': style5,    'marker': MARKER_SEQ[ 8 % MARKER_SEQ_LEN]},
     {'color': colorcycles[ 9 % CC_LEN], 'linestyle': style6,    'marker': MARKER_SEQ[ 9 % MARKER_SEQ_LEN]},
     {'color': colorcycles[10 % CC_LEN], 'linestyle': style7,    'marker': MARKER_SEQ[10 % MARKER_SEQ_LEN]},
     {'color': colorcycles[11 % CC_LEN], 'linestyle': style8,    'marker': MARKER_SEQ[11 % MARKER_SEQ_LEN]},
     {'color': colorcycles[12 % CC_LEN], 'linestyle': 'solid',   'marker': MARKER_SEQ[12 % MARKER_SEQ_LEN]},
]

linetypes = [
    {'color': colorcycles[i % CC_LEN], 'linestyle': linestyles[i % LS_LEN],   'marker': MARKER_SEQ[i % MARKER_SEQ_LEN]}
    for i in range(20)
]

MARKERS = ([m for m, func in Line2D.markers.items() if func != 'nothing' and m not in Line2D.filled_markers][::-1]) # ['o', 'x', '+', '']

full_labels = []

fig, ax1 = plt.subplots()
ax1.set_xlim(0, 1)
ax1.set_ylim(0, 1)

#ax1=plt

AX2LIM=0.945
AX3LIM=0.995

ax2 = plt.axes([0, 0, 1, 1], label = 'inset1')
ax2.set_xlim(AX2LIM, 1)
ax2.set_ylim(AX2LIM, 1)

ax3 = plt.axes([0, 0, 1, 1], label = 'inset2')
ax3.set_xlim(AX3LIM, 1)
ax3.set_ylim(AX3LIM, 1)

ax2.set_xticklabels(ax2.get_xticks(), backgroundcolor='w', fontsize=8-1)
ax2.set_yticklabels(ax2.get_yticks(), backgroundcolor='w', fontsize=8-1)
ax3.set_xticklabels(ax3.get_xticks(), backgroundcolor='w', fontsize=8-2)
ax3.set_yticklabels(ax3.get_yticks(), backgroundcolor='w', fontsize=8-2)

ip1 = InsetPosition(ax1, [0.2, 0.4, 0.5, 0.5])
ax2.set_axes_locator(ip1)

ip2 = InsetPosition(ax2, [0.2, 0.4, 0.5, 0.5])
ax3.set_axes_locator(ip2)

#mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')
is_inset1_marked = False #True
is_inset2_marked = False #True

#fig, ax = plt.subplots()
fidx = 0
str_base_vars = ''
str_vartype = ''
uvc_version = ''

def fname_cmpval(x):
    if ('uvc1.vcfeval-flt' in x or 'uvc1.vcfeval-filter' in x or 'uvc1-flt' in x or 'uvc1-flter' in x):
        return (2, x)
    elif 'uvc1.vcfeval-tFA5perc' in x:
        return (1, x)
    elif ('uvc1.vcfeval.outdir' in x):
        return (0, x)
    elif ('uvc1.vcfeval-' in x):
        return (4, x) # for test 
    else:
        return (3, x)

auc_best = fmeasure_max_best = 0
auc_best_prog = fmeasure_max_best_prog = 'NoVariantCaller'
auc_best_vcname = fmeasure_max_best_vcname = 'NoVariantCaller'

with codecs.open(outfname + '.tsv', 'w', 'utf-8') as tsvoutfile:
    fignum = 1
    fnames = sorted(sys.argv[4:], key = lambda x:fname_cmpval(x))
    labels = [get_label(fname, getbest) for fname in fnames]
    for fidx1, (ordinal, fname, legend_label) in enumerate(sorted(labels)):
        legend_label = legend_label.decode(encoding='UTF-8') 
        if fname.endswith('uvc1.version.info'):
            with open(fname, 'r') as file:
                uvc_version_2 = ('UVC_version=' if uvc_version == '' else ',') + file.next().strip()[0:7]
                # assert uvc_version == uvc_version_2 or uvc_version == '', '{} != {}'.format(uvc_version_2, uvc_version)
                uvc_version += uvc_version_2
            continue
        if not os.path.exists(fname):
            full_label = 'Runtime-error'
            full_labels.append((full_label + ' ' + legend_label + ' '))
            fidx += 1
            continue
        with gzip.open(fname, 'r') as file:
            print('Start processing {}'.format(fname))
            #legend_label = get_label(fname, getbest).decode(encoding='UTF-8')
            if '' == legend_label: 
                print('skipping {}'.format(fname))
                continue
            best_precision = best_sensitivity = fmeasure_max = -1
            prec_sens_list = [(1, 0)]
            auc = 0
            for line in file:
                if line.startswith('#'): 
                    if line.startswith('#total baseline variants'):
                        str_base_vars_2 = line.split()[-1]
                        if should_assert_baseline_var_num:
                            assert str_base_vars_2 == str_base_vars or str_base_vars == '', '{} != {}'.format(str_base_vars_2, str_base_vars)
                        elif not (str_base_vars_2 == str_base_vars or str_base_vars == ''):
                            print('WARNING: {} = {} failed for str_base_vars equality for the new file {}'.format(str_base_vars_2, str_base_vars, fname))
                        str_base_vars = str_base_vars_2
                    if line.startswith('#selection:'):
                        str_vartype_2 = line.split()[1].replace('SNP', 'SNV')
                        assert str_vartype_2 == str_vartype or str_vartype == ''
                        str_vartype = str_vartype_2
                    if line.startswith('#score field:'):
                        score_field = line.split()[2]
                        if score_field in ['QSS_QSI_NT', 'QSI_QSS_NT']:
                            score_field = ('QSS_NT' if str_vartype.startswith('SNV') else 'QSI_NT')
                    continue
                tokens = line.strip().split()
                sensitivity = float(tokens[-2])
                precision = float(tokens[-3])
                fmeasure = float(tokens[-1])
                if fmeasure > fmeasure_max:
                    best_precision, best_sensitivity = (precision, sensitivity)
                    fmeasure_max = fmeasure
                prev_prec, prev_sens = prec_sens_list[-1]
                auc += precision * (sensitivity - prev_sens)
                prec_sens_list.append((precision, sensitivity))
            if fmeasure_max > -0.5:
                partlabel = 'Fscore={:.4f} AUC={:.4f} '.format(fmeasure_max, auc)
            else:
                partlabel = 'Not-applicable '
            full_label = partlabel + legend_label + ' ' + score_field
            full_labels.append(full_label)
            tsvoutfile.write(full_label)
            tsvoutfile.write(' file=' + outfname + '\n');
            if fmeasure_max > fmeasure_max_best:
                fmeasure_max_best = fmeasure_max
                fmeasure_max_best_prog = full_label
                fmeasure_max_best_vcname = legend_label
            if auc > auc_best:
                auc_best = auc
                auc_best_prog = full_label
                auc_best_vcname = legend_label
            for i, (ax, xylim) in enumerate([(ax1, -1), (ax2, AX2LIM), (ax3, AX3LIM)]):
                prec_sens_list2 = [x for x in prec_sens_list if (x[0] >= xylim and x[1] >= xylim)];
                if len(prec_sens_list2) > 0:
                    prec_sens_list2.append((0, prec_sens_list[-1][1]))
                    
                    if i == 1 and not is_inset1_marked:
                        mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')
                        is_inset1_marked = True
                    if i == 2 and not is_inset2_marked:
                        mark_inset(ax2, ax3, loc1=2, loc2=4, fc="none", ec='0.5')
                        is_inset2_marked = True
                    if best_sensitivity >= xylim and best_precision >= xylim:
                        ax.scatter(best_sensitivity, best_precision, 
                                s=(1.25 * plt.rcParams['lines.markersize'])**2, 
                                color=linetypes[fidx]['color'], 
                                marker=linetypes[fidx]['marker'], 
                                label=full_label)
                    ax.plot([x[1] for x in prec_sens_list2], [x[0] for x in prec_sens_list2], #linewidth = 1.25, #label=legend_label,
                            color=linetypes[fidx]['color'], linestyle=linetypes[fidx]['linestyle'])
            fidx += 1
            #ax.annotate('X', (best_sensitivity, best_precision))
    #if fignum == 1:
    #    pass
    #    #ax = fig.add_subplot(122)
    #    #fig, ax = plt.subplots(figsize=[5, 4])
    #    #axins2 = zoomed_inset_axes(ax, 0.4, loc=1)
    #    #break

    try:
        outstring = (u"#Best-Fscore: {} for {} {}\n").format(fmeasure_max_best_prog, title, outfname)
        tsvoutfile.write(outstring)
        outstring = (u"#Best-PR-AUC: {} for {} {}\n").format(auc_best_prog, title, outfname)
        tsvoutfile.write(outstring)
    except Exception as err:
        sys.stderr.write(str(err))
        sys.stderr.write('\n\n')
        sys.stderr.write(fmeasure_max_best_prog + '\n')
        sys.stderr.write(auc_best_prog + '\n')
        sys.stderr.write(title + '\n')
        sys.stderr.write(outfname + '\n')

ax1.text(0.70, 0.025, ("Best Fscore: {}\n    Best AUC: {}").format(fmeasure_max_best_vcname, auc_best_vcname), fontsize = 8 - 3)

if not is_inset1_marked:
    ip0 = InsetPosition(ax1, [10, 10, 10, 10])
    ax2.set_axes_locator(ip0)
if not is_inset2_marked:
    ip0 = InsetPosition(ax1, [10, 10, 10, 10])
    ax3.set_axes_locator(ip0)

#plt.title(title.split(';')[0].decode(encoding='UTF-8'))
ax1.set_title(str_vartype + ' baseline_total=' + str_base_vars + '\n' + title.split(';')[0].replace('\\n', '\n'), fontdict = {'fontsize' : 8-1.5})
ax1.text(xlower, ylower, 'Evaluated {} callers'.format(fidx), fontsize=3)
ax1.legend(
    [create_dummy_line(**l) for l in linetypes[0:len(full_labels)]],
    full_labels,
    loc='lower left', markerscale=1.25, handlelength=5.25, framealpha=0.50, fontsize=plt.rcParams['font.size']-3)
#ax1.grid()
ax1.set_xlabel('Recall')
ax1.set_ylabel('Precision')

plt.tight_layout()
plt.savefig(outfname)
exit(0)

