from __future__ import division
import os
import re
import time
import argparse
import random

import compute_QBS_LHC as cql
import h5_storage
from h5_storage import data_dir

new_version = 7

parser = argparse.ArgumentParser()
parser.add_argument('-r', help='random', action='store_true')
args = parser.parse_args()

re_file = re.compile('cryo_data_fill_(\d{4,}).h5')

atd_files = os.listdir(data_dir)
if args.r:
    random.shuffle(atd_files)

for atd_file in atd_files:
    info = re_file.search(atd_file)
    if info is not None:
        filln = int(info.group(1))
        for use_dP in (True, False):
            kwargs = {'use_dP': use_dP, 'version': new_version}
            this_qbs_file = h5_storage.get_qbs_file(filln, **kwargs)
            if not os.path.isfile(this_qbs_file):
                time_0 = time.time()
                atd_ob = h5_storage.load_data_file(filln)
                qbs_ob = cql.compute_qbs(atd_ob, **kwargs)

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

