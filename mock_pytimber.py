import numpy as np
import pickle

from LHCMeasurementTools.LHC_Fills import Fills_Info
import LHCMeasurementTools.TimestampHelpers as TH

import GasFlowHLCalculator.qbs_fill as qf
from GasFlowHLCalculator.h5_storage import H5_storage

class MockPytimber(object):

    def __init__(self, data_folder_list, recalc_h5_folder):
        self.data_folder_list = data_folder_list
        self.recalc_h5_folder = recalc_h5_folder

    def get(self, variables, t1, t2):

        # merge pickles and add info on location
        dict_fill_bmodes={}
        for df in self.data_folder_list:
            with open(df+'/fills_and_bmodes.pkl', 'rb') as fid:
                this_dict_fill_bmodes = pickle.load(fid)
                for kk in this_dict_fill_bmodes:
                    this_dict_fill_bmodes[kk]['data_folder'] = df
                dict_fill_bmodes.update(this_dict_fill_bmodes)

        fill_info = Fills_Info(dict_fill_bmodes)
        fill_list = fill_info.fills_in_time_window(t1, t2)

        output = {}
        for vv in variables:
            output[vv] = [[],[]]

        for ff in sorted(fill_list):
            fdict = qf.get_fill_dict(ff,
                      h5_storage=H5_storage(self.recalc_h5_folder),
                      use_dP=True)
            for vv in variables:
                mask = (fdict[vv].t_stamps>=t1) & (fdict[vv].t_stamps<=t2)
                output[vv][0] += list(fdict[vv].t_stamps[mask]
                                    +1e-3*fdict[vv].ms[mask])
                output[vv][1] += list(fdict[vv].values[mask])

        for vv in variables:
            for ii in [0, 1]:
                output[vv][ii] = np.array(output[vv][ii])

        return(output)
