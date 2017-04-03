import h5py
import time
import os
import LHCMeasurementTools.TimberManager as tm
import LHCMeasurementTools.myfilemanager as mfm

# latest, default version (cell data only)
version = 7

# Directories
h5_dir = '/eos/user/l/lhcscrub/timber_data_h5/'
data_dir = h5_dir + '/cryo_heat_load_data/'
special_data_dir = h5_dir + 'cryo_special_cell_data/'
recalc_dir = h5_dir + 'recalculated_qbs/'

# Filenames for recomputed dada
def get_qbs_file(filln, use_dP, version=version):
    if use_dP:
        return recalc_dir + 'recalculated_qbs_v%i/recalculated_qbs_v%i_%i.h5' % (version, version, filln)
    else:
        if version == 7:
            version = 6
        return recalc_dir + '/recalculated_qbs_nodP_v%i/recalculated_qbs_nodP_v%i_%i.h5' % (version, version, filln)

def get_special_qbs_file(filln):
    return h5_dir + '/recalculated_special_qbs/recalculated_special_qbs_%i.h5' % filln


# Filename for raw data
def get_data_file(filln):
    return data_dir + '/cryo_data_fill_%i.h5' % filln

def get_special_data_file(filln):
    return special_data_dir + '/special_data_fill_%i.h5' % filln


# Load raw data
def load_data_file(filln):
    ob =  mfm.h5_to_obj(get_data_file(filln))
    return tm.AlignedTimberData(ob.timestamps, ob.data, ob.variables)

def load_special_data_file(filln):
    ob = mfm.h5_to_obj(get_special_data_file(filln))
    return tm.AlignedTimberData(ob.timestamps, ob.data, ob.variables)

# Load recomputed data
def load_qbs(filln, use_dP, version=version):
    qbs_file = get_qbs_file(filln, version=version, use_dP=use_dP)
    qbs_ob = mfm.h5_to_obj(qbs_file)
    return tm.AlignedTimberData(qbs_ob.timestamps, qbs_ob.data, qbs_ob.variables)

def load_special_qbs(filln):
    qbs_file = get_special_qbs_file(filln)
    qbs_ob = mfm.h5_to_obj(qbs_file)
    return tm.AlignedTimberData(qbs_ob.timestamps, qbs_ob.data, qbs_ob.variables)

# Store recomputed data
def store_qbs(filln, qbs_ob, use_dP, version=version):
    """
    Arguments:
        - filln
        - qbs_ob
        - use_dP
        - version=version
    """
    qbs_file = get_qbs_file(filln, version=version, use_dP=use_dP)
    if not os.path.isdir(os.path.dirname(qbs_file)):
        os.mkdir(os.path.dirname(qbs_file))

    with h5py.File(qbs_file, 'w') as h5_handle:
        h5_handle.create_dataset('timestamps', data=qbs_ob.timestamps)
        h5_handle.create_dataset('variables', data=qbs_ob.variables)
        qbs_dataset = h5_handle.create_dataset('data', data=qbs_ob.data)
        qbs_dataset.attrs['version'] = version
        qbs_dataset.attrs['with_dP'] = use_dP
        qbs_dataset.attrs['time_created'] = tm.UnixTimeStamp2UTCTimberTimeString(time.time())

def store_special_qbs(filln, qbs_ob):
    """
    Arguments:
        - filln
        - qbs_ob
    """
    qbs_file = get_special_qbs_file(filln)
    if not os.path.isdir(os.path.dirname(qbs_file)):
        os.mkdir(os.path.dirname(qbs_file))

    with h5py.File(qbs_file, 'w') as h5_handle:
        h5_handle.create_dataset('timestamps', data=qbs_ob.timestamps)
        h5_handle.create_dataset('variables', data=qbs_ob.variables)
        qbs_dataset = h5_handle.create_dataset('data', data=qbs_ob.data)
        qbs_dataset.attrs['time_created'] = tm.UnixTimeStamp2UTCTimberTimeString(time.time())

