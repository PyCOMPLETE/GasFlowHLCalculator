import h5py
import time
import os
import LHCMeasurementTools.TimberManager as tm
import LHCMeasurementTools.myfilemanager as mfm

class H5_storage(object):

    def __init__(self, version = 7, special_version = 2, 
            h5_dir = '/eos/user/l/lhcscrub/timber_data_h5/'):
        
        self.version = version
        self.special_version = special_version
        self.h5_dir = h5_dir + '/'
        
        self.data_dir = self.h5_dir + '/cryo_heat_load_data/'
        self.special_data_dir = self.h5_dir + '/cryo_special_cell_data/'
        self.recalc_dir = h5_dir + '/recalculated_qbs/'
    
    # Filenames for recomputed dada
    def get_qbs_file(self, filln, use_dP, version=None):

        if version is None:
            version = self.version
        
        if use_dP:
            return self.recalc_dir + 'recalculated_qbs_v%i/recalculated_qbs_v%i_%i.h5' % (version, version, filln)
        else:
            # version 6->7 only changed things for with dP
            if version == 7:
                version = 6
            return self.recalc_dir + '/recalculated_qbs_nodP_v%i/recalculated_qbs_nodP_v%i_%i.h5' % (version, version, filln)
    
    def get_special_qbs_file(self, filln, special_version=None):
        if special_version is None:
            special_version=self.special_version	

        if special_version == 0:
            return self.h5_dir + '/recalculated_special_qbs/recalculated_special_qbs_v%i/recalculated_special_qbs_%i.h5' % (special_version, filln)
        else:
            return self.h5_dir + '/recalculated_special_qbs/recalculated_special_qbs_v%i/recalculated_special_qbs_v%i_%i.h5' % (special_version, special_version, filln)
    
    
    # Filename for raw data
    def get_data_file(self, filln):
        return self.data_dir + '/cryo_data_fill_%i.h5' % filln
    
    def get_special_data_file(self, filln):
        return self.special_data_dir + '/special_data_fill_%i.h5' % filln
    
    
    # Load raw data
    def load_data_file(self, filln):
        ob =  mfm.h5_to_obj(self.get_data_file(filln))
        return tm.AlignedTimberData(ob.timestamps, ob.data, ob.variables)
    
    def load_special_data_file(self, filln):
        ob = mfm.h5_to_obj(self.get_special_data_file(filln))
        return tm.AlignedTimberData(ob.timestamps, ob.data, ob.variables)
    
    # Load recomputed data
    def load_qbs(self, filln, use_dP, version=None):
        if version is None:
            vesion = self.version
        qbs_file = self.get_qbs_file(filln, version=version, use_dP=use_dP)
        qbs_ob = mfm.h5_to_obj(qbs_file)
        #print('Loaded file %s' % qbs_file)
        return tm.AlignedTimberData(qbs_ob.timestamps, qbs_ob.data, qbs_ob.variables)
    
    def load_special_qbs(self, filln, special_version=None):
        if special_version is None:
            special_version = self.special_version
        qbs_file = self.get_special_qbs_file(filln, special_version=special_version)
        qbs_ob = mfm.h5_to_obj(qbs_file)
        return tm.AlignedTimberData(qbs_ob.timestamps, qbs_ob.data, qbs_ob.variables)
    
    
    # Store recomputed data
    def store_qbs(self, filln, qbs_ob, use_dP, version):
        """
        Arguments:
            - filln
            - qbs_ob
            - use_dP
            - version=version
        """
        qbs_file = self.get_qbs_file(filln, version=version, use_dP=use_dP)
        if not os.path.isdir(os.path.dirname(qbs_file)):
            os.mkdir(os.path.dirname(qbs_file))
    
        with h5py.File(qbs_file, 'w') as h5_handle:
            h5_handle.create_dataset('timestamps', data=qbs_ob.timestamps)
            h5_handle.create_dataset('variables', data=qbs_ob.variables)
            qbs_dataset = h5_handle.create_dataset('data', data=qbs_ob.data)
            qbs_dataset.attrs['version'] = version
            qbs_dataset.attrs['with_dP'] = use_dP
            qbs_dataset.attrs['time_created'] = tm.UnixTimeStamp2UTCTimberTimeString(time.time())
    
    def store_special_qbs(self, filln, qbs_ob, special_version):
        """
        Arguments:
            - filln
            - qbs_ob
            - special_version=special_version
        """
        qbs_file = self.get_special_qbs_file(filln, special_version=special_version)
        if not os.path.isdir(os.path.dirname(qbs_file)):
            os.mkdir(os.path.dirname(qbs_file))
    
        with h5py.File(qbs_file, 'w') as h5_handle:
            h5_handle.create_dataset('timestamps', data=qbs_ob.timestamps)
            h5_handle.create_dataset('variables', data=qbs_ob.variables)
            qbs_dataset = h5_handle.create_dataset('data', data=qbs_ob.data)
            qbs_dataset.attrs['time_created'] = tm.UnixTimeStamp2UTCTimberTimeString(time.time())

