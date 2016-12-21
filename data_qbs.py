import numpy as np
import csv
import sys
import os
from h5_storage import version

arc_list = ['ARC12','ARC23','ARC34','ARC45','ARC56','ARC67','ARC78','ARC81']
Radius = 3.7e-3/2.  #internal radius of beam screen cooling pipe (D=3.7 mm)
rug = 1e-5         #beam screen cooling circuit roughness (10 um)
arc_index = np.array(
      [[  5,  56],
       [ 68, 119],
       [125, 176],
       [186, 237],
       [249, 300],
       [307, 358],
       [364, 415],
       [427, 478]])

class Data_qbs(object):
    def __init__(self, version=version):
        csv_file_name = os.path.dirname(os.path.abspath(__file__)) + '/data_qbs_lhc_%i.csv' % version
        with open(csv_file_name, 'r') as f:
            data_qbs = {}
            ww = csv.reader(f, delimiter='\t')
            for ctr, row in enumerate(ww):
                if ctr == 0:
                    first_row = row
                    for ii, item in enumerate(first_row):
                        data_qbs[item] = []
                else:
                    for ii, item in enumerate(row):
                        data_qbs[first_row[ii]].append(item)

        for key, value in data_qbs.iteritems():
            if key != 'Sector_list': # Sector list remains a list of strings
                try:
                    value = np.array(value, float)
                except:
                    pass
            setattr(self, key, value)

        self.arc_index = arc_index
        self.arc_list = arc_list
        self.Radius = Radius
        self.rug = rug

## This is how a csv file can be created from python

#import data_QBS_LHC as dql
#list_to_save = [dql.__dict__[ss] for ss in variable_list]

#def save():
#    with open('data_qbs_lhc.csv', 'w') as f:
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

## This is how the arc_index is created:

#j = 0 #sector number
#        Type_list = self.Type_list
#        for i in xrange(len(Type_list)):
#            if Type_list[i-1] == 'LSS' and  Type_list[i] == 'ARC':   #begining of ARC
#                arc_index[j,0] = i
#            elif Type_list[i-1] == 'ARC' and Type_list[i] == 'LSS':  #end of ARC
#                arc_index[j,1] = i-1
#                j += 1
