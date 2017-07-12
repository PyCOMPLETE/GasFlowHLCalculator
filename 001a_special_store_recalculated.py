from __future__ import division
import os
import re
import time
import argparse
import random

from compute_QBS_special import compute_qbs_special
import qbs_fill as qf
import h5_storage
from h5_storage import special_data_dir

parser = argparse.ArgumentParser()
parser.add_argument('-r', help='random', action='store_true')
args = parser.parse_args()

re_file = re.compile('special_data_fill_(\d{4,}).h5')

atd_files = os.listdir(special_data_dir)
if args.r:
    random.shuffle(atd_files)

for atd_file in atd_files:
    info = re_file.search(atd_file)
    if info != None:
        filln = int(info.group(1))
        this_qbs_file = h5_storage.get_special_qbs_file(filln)
        if not os.path.isfile(this_qbs_file):
            time_0 = time.time()
            atd_ob = h5_storage.load_special_data_file(filln)
            new_cell = filln > 5600
            qbs_dict = compute_qbs_special(atd_ob, new_cell=new_cell)
            qbs_ob = qf.dict_to_aligned(qbs_dict)
            n_tries = 5
            while n_tries > 0:
                try:
                    h5_storage.store_special_qbs(filln, qbs_ob)
                    break
                except IOError:
                    n_tries -= 1
                    time.sleep(5)
            else:
                raise IOError('Saving failed for fill %i!' % filln)
            dt = time.time() - time_0
            n_timesteps = len(qbs_ob.timestamps)
            print('Calculation for fill %i with %i timesteps finished in %i s.' % (filln, n_timesteps, dt))

