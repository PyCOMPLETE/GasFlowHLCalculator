from __future__ import division
import sys
import os
import re

sys.path.append('..')
import LHCMeasurementTools.myfilemanager as mfm
import compute_QBS_LHC as cql
import h5_storage

use_dP = True

h5_dir = '/eos/user/l/lhcscrub/timber_data_h5/cryo_heat_load_data/'
os.chdir(h5_dir)
qbs_file = 'recalculated_qbs_v%i' % h5_storage.version + '_%i.h5'
re_file = re.compile('cryo_data_fill_(\d{4}).h5')

atd_files = os.listdir(h5_dir)
for atd_file in atd_files:
    info = re_file.search(atd_file)
    if info is not None:
        filln = int(info.group(1))
        this_qbs_file = qbs_file % filln
        if not os.path.isfile(this_qbs_file):
            print('Starting calculation for fill %i.' % filln)
            atd = mfm.h5_to_obj(atd_file)
            qbs_ob = cql.compute_qbs(atd, use_dP)
            h5_storage.store_qbs(filln, qbs_ob, use_dP)
            print('Calculation for fill %i saved.' % filln)
