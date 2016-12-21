import h5py
import time
import sys
import os
sys.path.append('..')
import LHCMeasurementTools.TimberManager as tm
import LHCMeasurementTools.myfilemanager as mfm

# latest version, default
version = 3
h5_dir = '/eos/user/l/lhcscrub/timber_data_h5/cryo_heat_load_data/'

def get_qbs_file(filln, version=version):
    return h5_dir + 'recalculated_qbs_v%i/recalculated_qbs_v%i_%i.h5' % (version, version, filln)

def store_qbs(filln, qbs_ob, use_dP, version=version):
    qbs_file = get_qbs_file(filln, version)
    if not os.path.isdir(os.path.dirname(qbs_file)):
        os.mkdir(os.path.dirname(qbs_file))

    with h5py.File(qbs_file, 'w') as h5_handle:
        h5_handle.create_dataset('timestamps', data=qbs_ob.timestamps)
        h5_handle.create_dataset('variables', data=qbs_ob.variables)
        qbs_dataset = h5_handle.create_dataset('data', data=qbs_ob.data)
        qbs_dataset.attrs['version'] = version
        qbs_dataset.attrs['with_dP'] = use_dP
        qbs_dataset.attrs['time_created'] = tm.UnixTimeStamp2UTCTimberTimeString(time.time())

def load_qbs(filln, version=version):
    qbs_file = get_qbs_file(filln, version)
    return mfm.h5_to_obj(qbs_file)
