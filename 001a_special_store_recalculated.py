from __future__ import division
import os
import re
import time
import argparse
import random

import compute_QBS_special as cqs
import h5_storage
from h5_storage import special_data_dir, special_version

parser = argparse.ArgumentParser()
parser.add_argument('-r', help='random', action='store_true')
parser.add_argument('--reverse', action='store_true')
args = parser.parse_args()

re_file = re.compile('special_data_fill_(\d{4,}).h5')

atd_files = os.listdir(special_data_dir)
if args.r:
    random.shuffle(atd_files)
elif args.reverse:
    atd_files = atd_files[::-1]

for atd_file in atd_files:
    info = re_file.search(atd_file)
    if info != None:
        filln = int(info.group(1))
        this_qbs_file = h5_storage.get_special_qbs_file(filln)
        if not os.path.isfile(this_qbs_file):
            time_0 = time.time()
            atd_ob = h5_storage.load_special_data_file(filln)
            new_cell = filln > 5600
            qbs_ob = cqs.compute_qbs_special(atd_ob, new_cell=new_cell, separate=True, aligned=True)
            n_tries = 5
            while n_tries > 0:
                try:
                    h5_storage.store_special_qbs(filln, qbs_ob, special_version)
                    break
                except IOError:
                    n_tries -= 1
                    time.sleep(5)
            else:
                raise IOError('Saving failed for fill %i!' % filln)
            dt = time.time() - time_0
            n_timesteps = len(qbs_ob.timestamps)
            print('Calculation for fill %i with %i timesteps finished in %i s.' % (filln, n_timesteps, dt))

