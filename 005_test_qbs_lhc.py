from __future__ import division
import os
import cPickle as pickle
import argparse

import matplotlib.pyplot as plt

from qbs_fill import compute_qbs_arc_avg
from config_qbs import arc_list
from compute_QBS_LHC import compute_qbs
import h5_storage

import LHCMeasurementTools.mystyle as ms

trusted_version = 5
filln = 5219

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

    ref_run_file = 'reference_%i_v%i_%s.pkl' % (filln, trusted_version, dp)

    atd_ob = h5_storage.load_data_file(filln)
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
        tt_ref = (qbs_ref.timestamps - qbs_ref.timestamps[0])

    fig = plt.figure()
    fig.canvas.set_window_title('QBS LHC reference %s' % dp)
    fig.patch.set_facecolor('w')
    sp = plt.subplot(2,1,1)
    for ctr, arc in enumerate(arc_list):
        color = ms.colorprog(ctr, arc_list)
        sp.plot(tt/3600., qbs_arc_avg[arc], lw=2, label=arc, color=color)
        if not reference_run:
            sp.plot(tt_ref/3600., ref_arc_avg[arc], lw=2, color=color, ls='--')
    sp.set_xlabel('time [hr]')
    sp.set_ylabel('Qdbs [W]')
    sp.set_title('Average beam screen heat load per ARC %s' % dp)
    sp.legend()
    sp.grid(True)
    sp = plt.subplot(2,1,2)
    sp.plot(tt/3600,qbs_ob.data)
    sp.set_xlabel('time [hr]')
    sp.set_ylabel('Qdbs [W]')
    sp.set_title('Beam screen heat loads over LHC')
    sp.grid(True)

    if reference_run:
        fig.savefig('./cell_qbs_reference_%i_%s.png' % (trusted_version, dp), dpi=200)

if not args.noshow:
    plt.show()
