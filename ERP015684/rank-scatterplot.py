# -*- coding: UTF-8 -*-

import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pylab import mpl

PATIENT_IDS='''21
236
245
282
300
469
575
610
664
708
726
735
772
781
815
861'''.strip().split()

plt.rcParams['font.family'] = ['sans-serif']
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 8

#print(plt.rcParams)
fig, ax1 = plt.subplots()

artq_x_uvcq_list_tfpair = ([], [], [])
for line in sys.stdin:
    tokens = line.split()
    if len(tokens) < 4: continue
    istp = -1
    if tokens[1] == '2': istp = 2
    if tokens[1] == '1': istp = 1
    if tokens[1] == '0': istp = 0
    if istp < 0: continue
    print('The line {} is processed'.format(line))
    artq_x_uvcq_list_tfpair[istp].append((float(tokens[2]), float(tokens[3])))

# MARKER_SEQ = 'sdo^<>vPxp*'

ax1.axhline(y = 60, linewidth=1, color='grey')
ax1.axvline(x = 0,  linewidth=1, color='grey')

x = [pt[0] for pt in artq_x_uvcq_list_tfpair[2]]
y = [pt[1] for pt in artq_x_uvcq_list_tfpair[2]]
ax1.scatter(x, y,
    label = 'True positive calls (n = {})'.format(len(x)),
    color = 'green', marker = 'o', alpha=0.25, s = 80)

'''
ax1.scatter(
        [pt[0] for pt in artq_x_uvcq_list_tfpair[1]], 
        [pt[1] for pt in artq_x_uvcq_list_tfpair[1]],
        label = 'Likely true positive calls',
        color = 'orange', marker = 'x', alpha=0.25, s = 80)
'''

x = [pt[0] for pt in artq_x_uvcq_list_tfpair[0]]
y = [pt[1] for pt in artq_x_uvcq_list_tfpair[0]]
ax1.scatter(x, y, 
        label = 'False positive calls (n = {})'.format(len(x)),
        color = 'red', marker = 'd', alpha=0.25, s = 80)

ax1.set_title('Project ERP015684\npatient IDs: ' + ', '.join(PATIENT_IDS), fontsize = 8)

ax1.set_xlabel('appreci8 artifact-score\n(higher means more likely to be false positive)')
ax1.set_ylabel('UVC-version-{} QUAL\n(higher means more likely to be true positive)'.format(sys.argv[2].split('-')[-1]))

plt.legend()
plt.tight_layout()
plt.savefig(sys.argv[1])
exit(0)

