import h5py
import time
import sys
import numpy as np
sys.path.append('..')
import LHCMeasurementTools.TimberManager as tm
import LHCMeasurementTools.myfilemanager as mfm

version = 2
h5_dir = '/eos/user/l/lhcscrub/timber_data_h5/cryo_heat_load_data/'
qbs_file = 'recalculated_qbs_v%i' % version + '_%i.h5'

def store_qbs(filln, qbs_ob, use_dP):
    with h5py.File(h5_dir+qbs_file % filln, 'w') as h5_handle:
        h5_handle.create_dataset('timestamps', data=qbs_ob.timestamps)
        h5_handle.create_dataset('variables', data=qbs_ob.variables)
        qbs_dataset = h5_handle.create_dataset('data', data=qbs_ob.data)
        qbs_dataset.attrs['version'] = version
        qbs_dataset.attrs['with_dP'] = use_dP
        qbs_dataset.attrs['time_created'] = tm.UnixTimeStamp2UTCTimberTimeString(time.time())

def load_qbs(filln):
    return mfm.h5_to_obj(h5_dir+qbs_file % filln)
