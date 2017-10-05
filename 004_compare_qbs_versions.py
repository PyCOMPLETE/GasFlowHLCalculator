from __future__ import division
import os
import re
import argparse

import matplotlib.pyplot as plt

import h5_storage
import qbs_fill as qf
from config_qbs import config_qbs
from compute_QBS_LHC import compute_qbs

import LHCMeasurementTools.mystyle as ms

parser = argparse.ArgumentParser()
parser.add_argument('filln', type=int)
parser.add_argument('--include-current', action='store_true', help='Recompute using the current code')
args = parser.parse_args()

plt.close('all')
ms.mystyle()

filln = args.filln
show_current = args.include_current

re_dir = re.compile('recalculated_qbs_(nodP_)?v(\d+)')
dirs = filter(re_dir.match, os.listdir(h5_storage.recalc_dir))

dir_dp_version = []
for dir_ in dirs:
    dir_dp_version.append((dir_,)+re_dir.match(dir_).groups())

sps = []
sp = None
for ctr, arc in enumerate(config_qbs.arc_list):
    sp_ctr = ctr%4 + 1
    if sp_ctr == 1:
        fig = ms.figure('Comparison of versions')
    sp = plt.subplot(2,2,sp_ctr, sharex=sp)
    sp.set_title(arc)
    sp.set_ylabel('Time [h]')
    sp.set_xlabel('Heat load [W/hc]')
    sp.grid(True)
    sps.append(sp)

for ctr, (dir_, nodp, version) in enumerate(dir_dp_version):
    use_dP = (nodp == None)
    if use_dP:
        label = 'V%s with dP' % version
    else:
        label = 'V%s without dP' % version
    qbs_ob = qf.compute_qbs_fill(filln, version=int(version), use_dP=use_dP, recompute_if_missing=True)
    arc_averages = qf.compute_qbs_arc_avg(qbs_ob)
    tt = (qbs_ob.timestamps - qbs_ob.timestamps[0]) / 3600.

    color = ms.colorprog(ctr, len(dir_dp_version)+2)
    for ctr, (arc, arr) in enumerate(sorted(arc_averages.dictionary.items())):
        sps[ctr].plot(tt, arr, label=label, lw=2, color=color)

# Current version
if show_current:
    atd = h5_storage.load_data_file(filln)
    qbs = compute_qbs(atd, use_dP=True)
    qbs_nodp = compute_qbs(atd, use_dP=False)

    for ctr, (qbs_ob, label) in enumerate(zip((qbs, qbs_nodp),('Current with dP', 'Current without dP'))):
        arc_averages = qf.compute_qbs_arc_avg(qbs_ob)
        tt = (qbs_ob.timestamps - qbs_ob.timestamps[0]) / 3600.

        color = ms.colorprog(ctr+len(dir_dp_version), len(dir_dp_version)+2)
        for ctr, (arc, arr) in enumerate(sorted(arc_averages.dictionary.items())):
            sps[ctr].plot(tt, arr, label=label, lw=2, color=color, ls='--')

for sp in sps[1::4]:
    sp.legend(bbox_to_anchor=(1.2,1))

plt.show()

