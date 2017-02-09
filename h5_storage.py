import h5py
import time
import os
import LHCMeasurementTools.TimberManager as tm
import LHCMeasurementTools.myfilemanager as mfm

# latest, default version
version = 5
h5_dir = '/eos/user/l/lhcscrub/timber_data_h5/'
data_dir = h5_dir + '/cryo_heat_load_data/'

def get_qbs_file(filln, version=version):
    return h5_dir + '/recalculated_qbs/recalculated_qbs_v%i/recalculated_qbs_v%i_%i.h5' % (version, version, filln)

def get_data_file(filln):
    return data_dir + '/cryo_data_fill_%i.h5' % filln

def load_data_file(filln):
    ob =  mfm.h5_to_obj(get_data_file(filln))
    return tm.AlignedTimberData(ob.timestamps, ob.data, ob.variables)

def get_special_data_file(filln):
    return data_dir + '/special_cells/special_data_fill_%i.h5' % filln

def load_special_data_file(filln):
    ob = mfm.h5_to_obj(get_special_data_file(filln))
    return tm.AlignedTimberData(ob.timestamps, ob.data, ob.variables)

def store_qbs(filln, qbs_ob, use_dP, version=version):
    """
    Arguments:
        -filln
        - qbs_ob
        - use_dP
        - version=version
    """
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
    qbs_ob = mfm.h5_to_obj(qbs_file)
    return tm.AlignedTimberData(qbs_ob.timestamps, qbs_ob.data, qbs_ob.variables)
