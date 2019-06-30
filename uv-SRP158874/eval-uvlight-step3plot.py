## # -*- coding: UTF-8 -*-

import sys
import gzip
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from matplotlib.lines import Line2D
from pylab import mpl

#plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 11
plt.tight_layout()
plt.subplots_adjust(bottom=0.17)
plt.subplots_adjust(left=0.15)
plt.subplots_adjust(top=0.92)
plt.subplots_adjust(right=0.96)

letter2inc = {'A': -0.24, 'C': -0.08, 'G': 0.08, 'T': 0.24}
letter2idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
letter2dat = {'A': [], 'C': [], 'G': [], 'T': []}

#>3:16306496-16306546
PROM_BEG = 16306496
PROM_SEQ = 'GGACTAGCCCTTCCGGCGCATAGGCAATGACGCAACTCCGCCCTGCG' #'GGACTAGCCCTTCCGGCGCATAGGCAATGACGCAACTCCGCCCTGCGCGGC'
PROM_END = PROM_BEG + len(PROM_SEQ) - 1

def samples2title(s1, s2):
    ss = s1 + '\t' + s2
    ret = ''
    if 'SRR7757437' in ss and 'SRR7757438' in ss:
        ret = 'XPA'
    if 'SRR7757439' in ss and 'SRR7757440' in ss:
        ret = 'CSA'
    if 'SRR7757441' in ss and 'SRR7757442' in ss:
        ret = 'CSB'
    if 'SRR7757443' in ss and 'SRR7757444' in ss:
        ret = 'XPC'
    return ret

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(1,1,1)
fmin = fmax = 0
for line in sys.stdin:
    if line.startswith('#CHROM'): 
        ptitle = samples2title(line.split()[-2], line.split()[-1])
    if line.startswith('##variantCallerVersion='): 
        uvc_version = line.split()[0].split('=')[1].split('.')[-1]
        uvc_version = (uvc_version[0:7] if ('dirty' not in uvc_version and '-' not in uvc_version and '_' not in uvc_version) else uvc_version)
    if line.startswith('#'): continue
    tokens = line.split()
    pos = int(tokens[1])
    if pos < PROM_BEG or pos > PROM_END: continue
    ref = tokens[3]
    alt = tokens[4]
    if len(ref) != len(alt) or ref == alt or alt not in letter2inc: continue
    idx = letter2idx[alt]
    info   = tokens[7]
    fmtkey = tokens[8]
    nfmtval = tokens[9]
    tfmtval = tokens[10]
    tDP = tAD = nDP = nAD = 0
    for k, tv, nv in zip(fmtkey.split(':'), tfmtval.split(':'), nfmtval.split(':')):
        if k == 'cDP2f' or k == 'cDP2r':
            tAD += float(tv.split(',')[1])
            tDP += float(tv.split(',')[0]) + float(tv.split(',')[1])
            nAD += float(nv.split(',')[1])
            nDP += float(nv.split(',')[0]) + float(nv.split(',')[1])
        if k == 'cADTT': 
            tAD += float(tv.split(',')[0] + tv.split(',')[1]) 
            nAD += float(nv.split(',')[0] + nv.split(',')[1])
        if k == 'cADT1' or k == 'cADTN': 
            tAD -= float(tv.split(',')[0] + tv.split(',')[1]) 
            nAD -= float(nv.split(',')[0] + nv.split(',')[1])
        if k == 'cDPTT': 
            tDP += float(tv.split(',')[0] + tv.split(',')[1]) 
            nDP += float(nv.split(',')[0] + nv.split(',')[1])
        if k == 'cDPT1' or k == 'cDPTN': 
            tDP -= float(tv.split(',')[0] + tv.split(',')[1]) 
            nDP -= float(nv.split(',')[0] + nv.split(',')[1])
    tf = tAD/float(tDP+sys.float_info.min) * 100.0
    nf = nAD/float(nDP+sys.float_info.min) * 100.0
    po = letter2inc[alt] + float(pos-PROM_BEG)
    letter2dat[alt].append((po,  tf))
    letter2dat[alt].append((po, -nf))
    fmin = min((fmin, -nf))
    fmax = max((fmax,  tf))

yticks = [0] #[-6e-2, -4e-2, -2e-2, 0, 2e-2, 4e-2, 6e-2, 8e-2]
fcur = 0
while fcur < fmax:
    fcur += 0.02
    yticks.append(fcur)
fcur = 0
while fcur > fmin:
    fcur -= 0.02
    yticks.append(fcur)
yticks = sorted(yticks)
yticklabels = ['{:.2f}'.format(abs(y)) for y in yticks]

plt.title(ptitle)

ax1.add_patch(matplotlib.patches.Rectangle((8-0.5,-100),1,200,linewidth=0, edgecolor='none',facecolor='gainsboro'))
ax1.set(ylim=(fmin*1.1, fmax*1.1))
ax1.text(0, 0, 'UVC_version='+uvc_version, fontsize=4, transform=ax1.transAxes)
ax1.set_xticks([i for i in range(len(PROM_SEQ))])
ax1.set_xticklabels(['{}\n{}'.format(base,(i if (i%4==0) else '')) for (i, base) in enumerate(PROM_SEQ)])
ax1.set_yticks(yticks)
ax1.set_yticklabels(yticklabels)
ax1.set_xticks([(i-0.5) for i in range(len(PROM_SEQ))], minor=True)
for letter, data in sorted(letter2dat.items()):
    x  = [ a[0] for a in data]
    y1 = [ a[1] for a in data]
    #print('{} {}'.format(len(x), len(y1)))
    ax1.bar(x, y1, 0.16, label=letter)
    print('Base {}, data {}'.format(letter, zip(x, y1)))
    print('Base {}, range {} to {}'.format(letter, min(y1), max(y1)))
ax1.set_xlabel("Base and its distance from the TSS of DPH3 gene in its promoter")
ax1.set_ylabel("Mutation frequency (%) for\n-UV control {} (at bottom) and {} +UV{} (at top)".format(''*5, ''*5, ''*25))
ax1.grid(True, which='minor', axis='x', linewidth=0.01)
ax1.grid(True, which='both',  axis='y', linewidth=0.01)
ax1.add_line(matplotlib.lines.Line2D([-100, 100], [0, 0], linewidth=0.1, linestyle=':', color='black'))
ax1.legend(loc = 'upper right', numpoints = 1, framealpha=0.25)
plt.savefig(sys.argv[1])

