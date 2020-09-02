import os
import json
import argparse
import time

import LHCMeasurementTools.lhc_log_db_query as lldb
from LHCMeasurementTools.SetOfHomogeneousVariables import SetOfHomogeneousNumericVariables
import LHCMeasurementTools.myfilemanager as mfm
from LHCMeasurementTools.LHC_Fill_LDB_Query import load_fill_dict_from_json

from GasFlowHLCalculator.h5_storage import H5_storage

import GasFlowHLCalculator

default_json_filename = None
h5_storage = H5_storage(h5_dir = '/eos/user/l/lhcecld/heatload_data_storage/')

# Config
dt_seconds = 60
max_fill_hrs = 35
blacklist = []
blacklist.append(4948) # 116 hour long fill, exceeds memory
blacklist.append(5488) # 40 hour long fill, also exceeds memory

parser = argparse.ArgumentParser()
parser.add_argument('-r', help='reversed', action='store_true')
parser.add_argument('--year', choices=[0, 2012, 2015, 2016, 2017, 2018, 2019], type=int, default=0)

args = parser.parse_args()
year = args.year

# File names
h5_dir_0 = h5_storage.data_dir

# [For all cells, for 3 special cells]
gfolder = '/'.join(GasFlowHLCalculator.__file__.split('/')[:-1])
variable_files = [gfolder+'/variable_list_complete.txt',
	gfolder+'/variable_list_special.txt']
h5_dirs = [h5_dir_0, h5_dir_0 + 'special_cells/']
file_names = ['cryo_data_fill', 'special_data_fill']
temp_filepaths = ['./tmp/' + f for f in file_names]
temp_files = [t + '_%i.csv' for t in temp_filepaths]
data_file_funcs = [h5_storage.get_data_file, h5_storage.get_special_data_file]

if year == 0:
    fills_json_name = default_json_filename
elif year == 2012:
    fills_json_name = '/afs/cern.ch/project/spsecloud/LHC_2012_selected_periods/fills_and_bmodes.json'
elif year == 2015:
    fills_json_name = '/afs/cern.ch/project/spsecloud/LHC_2015_PhysicsAfterTS2/fills_and_bmodes.json'
elif year == 2016:
    fills_json_name = '/afs/cern.ch/project/spsecloud/LHC_2016_25ns/LHC_2016_25ns_beforeTS1/fills_and_bmodes.json'
elif year == 2017:
    fills_json_name = '/afs/cern.ch/project/spsecloud/LHC_2017_operation/LHC_2017_operation/fills_and_bmodes.json'
elif year == 2018:
    fills_json_name = '/afs/cern.ch/work/l/lhcscrub/LHC_2018_followup/fills_and_bmodes.json'
elif year == 2019:
    fills_json_name = '/afs/cern.ch/work/l/lhcecld/run3_setup/LHC_followup_download_scripts/fills_and_bmodes.json'
else:
    raise ValueError('Invalid year')

dict_fill_bmodes = load_fill_dict_from_json(fills_json_name)

for variable_file, h5_dir, file_name, temp_filepath, temp_file, data_file_func in \
        zip(variable_files, h5_dirs, file_names, temp_filepaths, temp_files, data_file_funcs):

    fill_sublist = sorted(list(dict_fill_bmodes.keys()), reverse=args.r)
    fill_sublist_2 = []
    data_files = os.listdir(os.path.dirname(data_file_func(0)))
    for filln in fill_sublist:
        if os.path.basename(data_file_func(filln)) in data_files:
            pass
        elif filln in blacklist:
            print(('Fill %i is blacklisted' % filln))
        elif not dict_fill_bmodes[filln]['flag_complete']:
            print(('Fill %i is not completed' % filln))
        else:
            t_start_fill = dict_fill_bmodes[filln]['t_startfill']
            t_end_fill   = dict_fill_bmodes[filln]['t_endfill']
            fill_hrs = (t_end_fill - t_start_fill)/3600.
            if fill_hrs < max_fill_hrs:
                fill_sublist_2.append(filln)
            else:
                print(('Fill %i exceeds %i hours and is skipped' % (filln, max_fill_hrs)))
    print(('Processing %i fills!' % len(fill_sublist_2)))
    time.sleep(5)

    with open(variable_file, 'r') as f:
        varlist = f.read().splitlines()[0].split(',')
    print('%i Timber variables' % len(varlist))
    for ii, var in enumerate(varlist):
        if '.POSST' not in var:
            raise ValueError('%s does not have a .POSST' % var)

    if not os.path.isdir(h5_dir):
        os.mkdir(h5_dir)

    for filln in fill_sublist_2:
        h5_file = data_file_func(filln)
        if h5_file in os.listdir(os.path.dirname(h5_file)):
            continue
        this_temp_file = temp_file % filln
        print(('Downloading csv for fill %i' % filln))
        t_start_fill = dict_fill_bmodes[filln]['t_startfill']
        t_end_fill   = dict_fill_bmodes[filln]['t_endfill']
        lldb.dbquery(varlist, t_start_fill, t_end_fill, this_temp_file)
        print(('Aligning data for fill %i' % filln))
        htd_ob = SetOfHomogeneousNumericVariables(varlist, this_temp_file).aligned_object(dt_seconds)
        print(('Creating h5 file for fill %i' % filln))
        n_tries_max = 5
        for n_try in range(n_tries_max):
            try:
                mfm.aligned_obj_to_h5(htd_ob, h5_file)
                break
            except Exception as e:
                print('Saving of h5 failed')
                time.sleep(10)
        else:
            print(('Raise error after trying to save the h5 file %i times' % n_tries_max))
            raise

        if os.path.isfile(h5_file) and os.path.getsize(h5_file) > 500:
            os.remove(this_temp_file)
            print(('Deleted temporary file %s!' % (this_temp_file)))
        else:
            print(('Warning! Something went wrong for file %s!\nKeeping temporary file %s.' % (h5_file % filln, temp_file % filln)))

