from __future__ import division
import os
import re
import time
import argparse
import random

import LHCMeasurementTools.myfilemanager as mfm
import LHCMeasurementTools.TimberManager as tm
import compute_QBS_LHC as cql
import h5_storage

use_dP = True

parser = argparse.ArgumentParser()
parser.add_argument('-r', help='reversed', action='store_true')

args = parser.parse_args()

h5_dir = '/eos/user/l/lhcscrub/timber_data_h5/cryo_heat_load_data/'
curr_dir = os.getcwd()
os.chdir(h5_dir)
try:
    re_file = re.compile('cryo_data_fill_(\d{4}).h5')

    atd_files = os.listdir(h5_dir)
    if args.r:
        random.shuffle(atd_files)

    for atd_file in atd_files:
        info = re_file.search(atd_file)
        if info is not None:
            filln = int(info.group(1))
            this_qbs_file = h5_storage.get_qbs_file(filln)
            if not os.path.isfile(this_qbs_file):
                time_0 = time.time()
                atd = mfm.h5_to_obj(atd_file)
                atd = tm.atd_from_h5(atd)
                qbs_ob = cql.compute_qbs(atd, use_dP)
                h5_storage.store_qbs(filln, qbs_ob, use_dP)
                dt = time.time() - time_0
                n_timesteps = len(qbs_ob.timestamps)
                print('Calculation for fill %i with %i timesteps finished in %i s.' % (filln, n_timesteps, dt))
finally:
    os.chdir(curr_dir)
