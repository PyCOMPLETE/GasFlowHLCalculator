import numpy as np
import csv
import os


arc_list = ['S12','S23','S34','S45','S56','S67','S78','S81']
Radius = 3.7e-3/2.  #internal radius of beam screen cooling pipe
rug = 1e-5         #beam screen cooling circuit roughness
arc_index = np.array(
      [[5,  56],
       [68, 119],
       [125, 176],
       [186, 237],
       [249, 300],
       [307, 358],
       [364, 415],
       [427, 478]])

def get_config_file():
    return '/config_qbs_lhc_3.csv'

def get_delimiter():
    return '\t'

class Config_qbs(object):
    def __init__(self):
        csv_file_name = os.path.dirname(os.path.abspath(__file__)) + get_config_file()
        with open(csv_file_name, 'r') as f:
            config_qbs = []
            tsv = csv.reader(f, delimiter=get_delimiter())
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
                except ValueError:
                    pass
            setattr(self, key, value)

        if get_delimiter() == ',':
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

        self.assert_correct_05L4_05R4()
        self.assert_arc_index()

        # Handle duplicates in Cell_list
        for ii, cell_name in enumerate(self.Cell_list[:]):
            if cell_name in ('05R4_947', '05L4_947'):
                if self.CV94x_list[ii].startswith('QRLEB'):
                    self.Cell_list[ii] += '_comb'
                elif self.CV94x_list[ii][:5] in ('QRLFE', 'QRLFF'):
                    self.Cell_list[ii] += '_quad'
                else:
                    raise ValueError('Illegal entry for %s: %s' % (cell_name, self.CV94x_list[ii]))


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


    def assert_arc_index(self):
        """
        Assert that the arc indices, which determine at which cells the arcs begin and end, remain unchanged in future versions of the csv file.
        If it changes, this will at least not go unnoticed.
        """
        arc_index_2 = np.zeros_like(arc_index)
        j = 0 #sector number
        Type_list = self.Type_list
        for i in xrange(len(Type_list)):
            if Type_list[i-1] == 'LSS' and Type_list[i] == 'ARC':   #begining of ARC
                arc_index_2[j,0] = i
            elif Type_list[i-1] == 'ARC' and Type_list[i] == 'LSS':  #end of ARC
                arc_index_2[j,1] = i-1
                j += 1
        assert np.all(arc_index_2 == arc_index)

    def assert_correct_05L4_05R4(self):
        """
        Make sure that the correct (QRLEB) cells come first and second in
        case of 05R4 and 05L4. Also in future versions of the config qbs
        objects (as from configuration file)
        """
        index_R = self.Cell_list.index('05R4_947')
        index_L = self.Cell_list.index('05L4_947')

        assert 'QRLEB' in self.CV94x_list[index_R]
        assert 'QRLEB' in self.CV94x_list[index_L+1]


