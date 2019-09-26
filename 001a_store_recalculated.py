from __future__ import division
import os
import re
import time
import argparse
import random

import compute_QBS_LHC as cql

from GasFlowHLCalculator.h5_storage import H5_storage

h5_storage = H5_storage(h5_dir = '/eos/user/l/lhcecld/heatload_data_storage/')
calibration_csv_file = '/afs/cern.ch/work/l/lhcecld/tools/LHCCryoHeatLoadCalibration/CryoBeamScreenData.csv'

data_dir = h5_storage.data_dir

use_dPs = (True,False)

parser = argparse.ArgumentParser()
parser.add_argument('-r', help='random', action='store_true')
parser.add_argument('--filln', help='specify fill number')
args = parser.parse_args()

re_file = re.compile('cryo_data_fill_(\d{4,}).h5')

atd_files = os.listdir(data_dir)
if args.r:
    random.shuffle(atd_files)

for atd_file in atd_files:
    info = re_file.search(atd_file)
    if info is not None:
        filln = int(info.group(1))

    if args.filln:
        if not int(filln)==int(args.filln):
            #print 'Skipped fill', filln
            continue

    for use_dP in use_dPs:
        kwargs = {'use_dP': use_dP}
        this_qbs_file = h5_storage.get_qbs_file(filln, **kwargs)
        if not os.path.isfile(this_qbs_file):

            print('\nCalculation for fill %i (usedP: %s) started...' % (filln, use_dP))
            time_0 = time.time()
            atd_ob = h5_storage.load_data_file(filln)
            qbs_ob = cql.compute_qbs(atd_ob,
                        calibration_csv_file=calibration_csv_file, **kwargs)

            n_tries = 5
            while n_tries > 0:
                try:
                    h5_storage.store_qbs(filln, qbs_ob, **kwargs)
                    break
                except IOError:
                    n_tries -= 1
                    time.sleep(60)
            else:
                raise IOError('Saving failed for fill %i!' % filln)
            dt = time.time() - time_0
            n_timesteps = len(qbs_ob.timestamps)
            print('Calculation for fill %i (usedP: %s) with %i timesteps finished in %i s.' % (filln, use_dP, n_timesteps, dt))

