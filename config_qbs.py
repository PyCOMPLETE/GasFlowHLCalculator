import numpy as np
import csv
import os
from h5_storage import version

arc_list = ['S12','S23','S34','S45','S56','S67','S78','S81']
Radius = 3.7e-3/2.  #internal radius of beam screen cooling pipe
rug = 1e-5         #beam screen cooling circuit roughness
arc_index = np.array(
      [[  5,  56],
       [ 68, 119],
       [125, 176],
       [186, 237],
       [249, 300],
       [307, 358],
       [364, 415],
       [427, 478]])

latest_config_file_version=3
def get_config_file(version):
    if version == -1:
        return '/LHCCryoHeatLoadCalibration/CryoBeamScreenData_beforeLS1_arcs.csv'
    if version == 2:
        return '/config_qbs_lhc_2.csv'
    elif 2 < version < 8:
        return '/config_qbs_lhc_3.csv'
    elif version == 8:
        return '/LHCCryoHeatLoadCalibration/CryoBeamScreenData.csv'
    else:
        raise ValueError('Config file not defined!')

def get_delimiter(version):
    if version == -1:
        return ','
    if version < 8:
        return '\t'
    else:
        return ','


class Config_qbs(object):
    def __init__(self, version=version):
        csv_file_name = os.path.dirname(os.path.abspath(__file__)) + get_config_file(version)
        with open(csv_file_name, 'r') as f:
            config_qbs = []
            tsv = csv.reader(f, delimiter=get_delimiter(version))
            for ctr, row in enumerate(tsv):
                if ctr == 0:
                    first_row = row
                    for ii in xrange(len(first_row)):
                        config_qbs.append([])
                else:
                    for ii, item in enumerate(row):
                        config_qbs[ii].append(item)

        for key, value in zip(first_row, config_qbs):
            if key != 'Sector_list': # Sector list remains a list of strings
                try:
                    value = np.array(value, float)
                except:
                    pass
            setattr(self, key, value)

        if get_delimiter(version) == ',':

            self.Sector_list    = self.Sector
            self.Type_list      = self.Type
            self.Cell_list      = self.Loop
            self.EH84x_list     = self.QEH
            self.TT84x_list     = self.T2
            self.CV94x_list     = self.CV1
            self.PT961_list     = self.P1
            self.PT991_list     = self.P4
            self.TT94x_list     = self.T3
            self.TT961_list     = self.T1
            self.R_list         = self.R
            self.Qs_list        = self.QS
            self.Kv_list        = self.Kvmax
            self.nc_list        = self.nc
            self.L_list         = self.L

        self.arc_index = arc_index
        self.arc_list = arc_list
        self.Radius = Radius
        self.rug = rug
        

    def get_varnames_for_cell(self,cell):
        index = self.Cell_list.index(cell)
        

        cq = self
        var_data_dict = {
            'T1': {
                'vars': cq.TT961_list,
                'correct': True,
                'negative': False,
            },
            'T3': {
                'vars': cq.TT94x_list,
                'correct': False,
                'negative': False,
            },
            'CV': {
                'vars': cq.CV94x_list,
                'correct': False,
                'negative': False,
            },
            'EH': {
                'vars': cq.EH84x_list,
                'correct': False,
                'negative': True,
            },
            'P1': {
                'vars': cq.PT961_list,
                'correct': True,
                'negative': False,
            },
            'P4': {
                'vars': cq.PT991_list,
                'correct': True,
                'negative': False,
            },
            'T2': {
                'vars': cq.TT84x_list,
                'correct': False,
                'negative': False,
            },
        }
        
        out = {}
        for attr in var_data_dict:
            out[attr] = var_data_dict[attr]['vars'][index]
        return out

# Default object
config_qbs = Config_qbs()

def assert_arc_index(config_qbs=config_qbs):
    arc_index_2 = np.zeros_like(arc_index)
    j = 0 #sector number
    Type_list = config_qbs.Type_list
    for i in xrange(len(Type_list)):
        if Type_list[i-1] == 'LSS' and  Type_list[i] == 'ARC':   #begining of ARC
            arc_index_2[j,0] = i
        elif Type_list[i-1] == 'ARC' and Type_list[i] == 'LSS':  #end of ARC
            arc_index_2[j,1] = i-1
            j += 1
    assert np.all(arc_index_2 == arc_index)

def assert_correct_05L4_05R4(config_qbs=config_qbs):
    index_R = config_qbs.Cell_list.index('05R4_947')
    index_L = config_qbs.Cell_list.index('05L4_947')

    # Make sure that the correct (QRLEB) cells come first and second in case of 05R4 and 05L4
    # also in future versions of the config qbs objects (as from configuration file)
    assert 'QRLEB' in config_qbs.CV94x_list[index_R]
    assert 'QRLEB' in config_qbs.CV94x_list[index_L+1]


assert_arc_index()
assert_correct_05L4_05R4()


## This is how a csv file can be created from python

#import data_QBS_LHC as dql
#list_to_save = [dql.__dict__[ss] for ss in variable_list]

#def save():
#    with open('config_qbs_lhc.csv', 'w') as f:
#        ww = csv.writer(f, delimiter='\t')
#        ww.writerow(variable_list)
#        for items in zip(*list_to_save):
#            ww.writerow(items)
#    variable_list = [
#            'Cell_list',
#            'Type_list',
#            'Sector_list',
#            'EH84x_list',
#            'TT84x_list',
#            'CV94x_list',
#            'PT961_list',
#            'PT991_list',
#            'TT94x_list',
#            'TT961_list',
#            'R_list',
#            'Qs_list',
#            'Kv_list',
#            'nc_list',
#            'L_list'
#            ]

