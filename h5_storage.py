import h5py
import time
import os

import numpy as np

import LHCMeasurementTools.TimberManager as tm
import LHCMeasurementTools.myfilemanager as mfm

def decode_if_needed(ss):
    if hasattr(ss, 'decode'):
        return ss.decode('utf-8')
    else:
        return ss

class H5_storage(object):

    def __init__(self, h5_dir):

        self.h5_dir = h5_dir + '/'

        self.data_dir = self.h5_dir + '/cryo_heat_load_data/'
        self.special_data_dir = self.h5_dir + '/cryo_special_cell_data/'
        self.recalc_dir = h5_dir + '/recalculated_qbs/'

    # Filenames for recomputed dada
    def get_qbs_file(self, filln, use_dP):

        if use_dP:
            return self.recalc_dir + 'recalculated_qbs/recalculated_qbs_%i.h5' % (filln)
        else:
            return self.recalc_dir + '/recalculated_qbs_nodP/recalculated_qbs_nodP_%i.h5' % (filln)

    def get_special_qbs_file(self, filln):

        return self.h5_dir + '/recalculated_special_qbs/recalculated_special_qbs/recalculated_special_qbs_%i.h5' % (filln)

    # Filename for raw data
    def get_data_file(self, filln):
        return self.data_dir + '/cryo_data_fill_%i.h5' % filln

    def get_special_data_file(self, filln):
        return self.special_data_dir + '/special_data_fill_%i.h5' % filln


    # Load raw data
    def load_data_file(self, filln):
        ob =  mfm.h5_to_obj(self.get_data_file(filln))
        variables = [decode_if_needed(vv) for vv in ob.variables]
        return tm.AlignedTimberData(ob.timestamps, ob.data, variables)

    def load_special_data_file(self, filln):
        ob = mfm.h5_to_obj(self.get_special_data_file(filln))
        variables = [decode_if_needed(vv) for vv in ob.variables]
        return tm.AlignedTimberData(ob.timestamps, ob.data, variables)

    # Load recomputed data
    def load_qbs(self, filln, use_dP):
        qbs_file = self.get_qbs_file(filln, use_dP=use_dP)
        qbs_ob = mfm.h5_to_obj(qbs_file)
        #print('Loaded file %s' % qbs_file)
        variables = [decode_if_needed(vv) for vv in qbs_ob.variables]
        return tm.AlignedTimberData(qbs_ob.timestamps, qbs_ob.data, variables)

    def load_special_qbs(self, filln):
        qbs_file = self.get_special_qbs_file(filln)
        qbs_ob = mfm.h5_to_obj(qbs_file)
        variables = [decode_if_needed(vv) for vv in qbs_ob.variables]
        return tm.AlignedTimberData(qbs_ob.timestamps, qbs_ob.data, variables)


    # Store recomputed data
    def store_qbs(self, filln, qbs_ob, use_dP):
        """
        Arguments:
            - filln
            - qbs_ob
            - use_dP
        """
        qbs_file = self.get_qbs_file(filln, use_dP=use_dP)
        if not os.path.isdir(os.path.dirname(qbs_file)):
            os.mkdir(os.path.dirname(qbs_file))

        with h5py.File(qbs_file, 'w') as h5_handle:
            h5_handle.create_dataset('timestamps', data=qbs_ob.timestamps)
            h5_handle.create_dataset('variables', data=np.string_(
                        qbs_ob.variables))
            qbs_dataset = h5_handle.create_dataset('data', data=qbs_ob.data)
            qbs_dataset.attrs['with_dP'] = use_dP
            qbs_dataset.attrs['time_created'] = tm.UnixTimeStamp2UTCTimberTimeString(time.time())

    def store_special_qbs(self, filln, qbs_ob):
        """
        Arguments:
            - filln
            - qbs_ob
        """
        qbs_file = self.get_special_qbs_file(filln)
        if not os.path.isdir(os.path.dirname(qbs_file)):
            os.mkdir(os.path.dirname(qbs_file))

        with h5py.File(qbs_file, 'w') as h5_handle:
            h5_handle.create_dataset('timestamps', data=qbs_ob.timestamps)
            h5_handle.create_dataset('variables',
                    data=np.string_(qbs_ob.variables))
            qbs_dataset = h5_handle.create_dataset('data', data=qbs_ob.data)
            qbs_dataset.attrs['time_created'] = tm.UnixTimeStamp2UTCTimberTimeString(time.time())

