
import os
import re
import time
import argparse
import random

from GasFlowHLCalculator.h5_storage import H5_storage
import GasFlowHLCalculator.recalc_multiple_circuits as rmc

from GasFlowHLCalculator.calibration_config import calibration_config
from GasFlowHLCalculator.calibration import Calibration, CalibrationManager

h5_storage = H5_storage(h5_dir = '/home/kparasch/workspace/Instrumented_HL_calc/heatload_data_storage/')
cal_manager = CalibrationManager(calibration_config=calibration_config)

special_data_dir = h5_storage.special_data_dir

parser = argparse.ArgumentParser()
parser.add_argument('-r', help='random', action='store_true')
parser.add_argument('--reverse', action='store_true')
parser.add_argument('--filln', help='specify fill number')
args = parser.parse_args()

re_file = re.compile('special_data_fill_(\d{4,}).h5')

atd_files = os.listdir(special_data_dir)
if args.r:
    random.shuffle(atd_files)
elif args.reverse:
    atd_files = atd_files[::-1]
print(atd_files)
for atd_file in atd_files:
    info = re_file.search(atd_file)
    if info != None:
        filln = int(info.group(1))

        if args.filln:
            if not int(filln)==int(args.filln):
                #print 'Skipped fill', filln
                continue

        this_qbs_file = h5_storage.get_special_qbs_file(filln)
        if not os.path.isfile(this_qbs_file):
            time_0 = time.mktime(time.localtime())
            atd_ob = h5_storage.load_special_data_file(filln)

            calibration = cal_manager.get_calibration(atd_ob.timestamps[0])

            qbs_ob, other = rmc.recalc_multiple_circuits(atd_ob,
                calibration, circuit_selection='all_instrumented',
                with_P_drop=True)

            n_tries = 5
            while n_tries > 0:
                try:
                    h5_storage.store_special_qbs(filln, qbs_ob)
                    break
                except IOError as err:
                    n_tries -= 1
                    print(err)
                    print("Retrying...")
                    time.sleep(5)
            else:
                raise IOError('Saving failed for fill %i!' % filln)
            dt = time.mktime(time.localtime()) - time_0
            n_timesteps = len(qbs_ob.timestamps)
            print(('Calculation for fill %i with %i timesteps finished in %i s.' % (filln, n_timesteps, dt)))

