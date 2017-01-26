from __future__ import division
import sys
import os
import cPickle as pickle
import argparse

import matplotlib.pyplot as plt

from qbs_fill import compute_qbs_arc_avg
from data_qbs import arc_list
from compute_QBS_LHC import compute_qbs
import h5_storage

if '..' not in sys.path: sys.path.append('..')
import LHCMeasurementTools.TimberManager as tm
import LHCMeasurementTools.mystyle as ms

trusted_version = 3

parser = argparse.ArgumentParser()
parser.add_argument('--noshow', help='Do not call plt.show().', action='store_true')
parser.add_argument('--ref', help='Create pickle for reference run.', action='store_true')
args = parser.parse_args()
reference_run = args.ref

plt.close('all')

# Only set to true if you want to create the pickle for a new version

for use_dP in (True, False):

    if use_dP:
        dp = 'with_dP'
    else:
        dp = 'without_dP'

    ref_run_file = 'reference_5416_v%i_%s.pkl' % (trusted_version, dp)

    filename = os.path.abspath(os.path.dirname(__file__)) + '/TIMBER_DATA_Fill5416_LHCBEAMSCREEN_TT84x_injec.csv'
    atd_ob = tm.parse_aligned_csv_file(filename)
    tt = atd_ob.timestamps - atd_ob.timestamps[0]

    qbs_ob = compute_qbs(atd_ob, use_dP, strict=False)
    qbs_arc_avg = compute_qbs_arc_avg(qbs_ob).dictionary

    if reference_run:
        if os.path.isfile(ref_run_file):
            raise ValueError('Reference file already exists!')
        with open(ref_run_file, 'w') as f:
            pickle.dump(qbs_ob, f, protocol=-1)
    else:
        with open(ref_run_file, 'r') as f:
            qbs_ref = pickle.load(f)
        ref_arc_avg = compute_qbs_arc_avg(qbs_ref).dictionary

    fig = plt.figure()
    fig.canvas.set_window_title('QBS LHC reference %s' % dp)
    fig.patch.set_facecolor('w')
    sp = plt.subplot(2,1,1)
    for ctr, arc in enumerate(arc_list):
        color = ms.colorprog(ctr, arc_list)
        sp.plot(tt/3600., qbs_arc_avg[arc], lw=2, label=arc, color=color)
        if not reference_run:
            sp.plot(tt/3600., ref_arc_avg[arc], lw=2, color=color, ls='--')
    sp.set_xlabel('time [hr]')
    sp.set_ylabel('Qdbs [W]')
    sp.set_title('Average beam screen heat load per ARC %s' % dp)
    sp.legend()
    sp.grid(True)
    sp = plt.subplot(2,1,2)
    sp.plot(tt/3600,qbs_ob.data)
    sp.set_xlabel('time (hr)')
    sp.set_ylabel('Qdbs (W)')
    sp.set_title('Beam screen heat loads over LHC')
    sp.grid(True)

    fig.savefig('./cell_qbs_reference_%s.png' % dp, dpi=200)

if not args.noshow:
    plt.show()