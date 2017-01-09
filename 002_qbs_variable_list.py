import sys
if '..' not in sys.path: sys.path.append('..')

import LHCMeasurementTools.myfilemanager as  mfm

filln = 5219

h5_filename = '/eos/user/l/lhcscrub/timber_data_h5/cryo_heat_load_data/cryo_data_fill_%i.h5' % filln

atd = mfm.h5_to_obj(h5_filename)

variable_list = atd.variables

with open('./variable_list_complete.txt', 'w') as f:
    f.write(','.join(variable_list))

