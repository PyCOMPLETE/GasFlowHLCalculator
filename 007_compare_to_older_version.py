from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

import qbs_fill as qf
import LHCMeasurementTools.mystyle as ms
import LHCMeasurementTools.LHC_Heatloads as HL
import LHCMeasurementTools.TimberManager as tm
from LHCMeasurementTools.SetOfHomogeneousVariables import SetOfHomogeneousNumericVariables

ms.mystyle_arial()

plt.close('all')
filln = 5219

version = 5
old_version = 4

arc_key_list = HL.variable_lists_heatloads['AVG_ARC']
quad_key_list = []
for key, list_ in HL.variable_lists_heatloads.iteritems():
    if key[:2] in ('Q6', 'Q5'):
        quad_key_list.extend(list_)

qbs_ob = qf.compute_qbs_fill(filln, version=version)
qbs_ob_old = qf.compute_qbs_fill(filln, version=old_version)

# Logged data
fill_dict = {}
fill_dict.update(tm.parse_timber_file('../fill_heatload_data_csvs/heatloads_fill_%d.csv' % filln, verbose=False))
heatloads = SetOfHomogeneousNumericVariables(variable_list=arc_key_list, timber_variables=fill_dict)

arc_data = qf.compute_qbs_arc_avg(qbs_ob)
arc_data_old = qf.compute_qbs_arc_avg(qbs_ob_old)

def tt_hrs(arr):
    return (arr - arr[0])/3600.

tt = tt_hrs(arc_data.timestamps)

# Arcs
fig = ms.figure('Comparison of arc and quad heat load for fill %i' % filln)
fig.subplots_adjust(left=0.05, wspace=0.37)
sp = plt.subplot(2,2,1)
sp.grid('on')
sp.set_ylabel('Heat load [W]')
sp.set_xlabel('Time [h]')

for ctr, (arc, data) in enumerate(arc_data.dictionary.iteritems()):
    color = ms.colorprog(ctr, arc_data.dictionary)
    sp.plot(tt, data, label=arc+' v%i'%version, color=color)

    if ctr == 0:
        label_old = 'v%i' % old_version
        label_logged = 'logged'
    else:
        label_old, label_logged = None, None

    sp.plot(tt, arc_data_old.dictionary[arc], ls ='--', color=color, label=label_old)

    for key in arc_key_list:
        if arc in key:
            sp.plot(tt_hrs(heatloads.timber_variables[key].t_stamps), heatloads.timber_variables[key].values, color=color, ls=':', label=label_logged)

sp.legend(bbox_to_anchor=(1.2,1))

# Quads
sp_5 = plt.subplot(2,2,2)
sp_5.set_title('Q5')

sp_6 = plt.subplot(2,2,4)
sp_6.set_title('Q6')
for sp in sp_5, sp_6:
    sp.grid('on')
    sp.set_ylabel('Heat load [W]')
    sp.set_xlabel('Time [h]')

fill_dict_v4 = qf.get_fill_dict(qbs_ob_old)
fill_dict_v5 = qf.get_fill_dict(qbs_ob)

for fd, title, ls in zip((fill_dict, fill_dict_v4, fill_dict_v5), ('logged', 'v4', 'v5'), (':', '--', '-')):
    heatloads = SetOfHomogeneousNumericVariables(variable_list=quad_key_list, timber_variables=fd)
    for ctr, (key, tvl) in enumerate(heatloads.timber_variables.iteritems()):
        if key[7] == '5':
            sp = sp_5
        elif key[7] == '6':
            sp = sp_6
        else:
            raise ValueError
        if ls == '-' or ctr == 0:
            label = key[6:10]
        else:
            label = None

        color = ms.colorprog(ctr, heatloads.timber_variables)

        sp.plot(tt_hrs(tvl.t_stamps), tvl.values, ls=ls, label=label, color=color)

for sp in sp_5, sp_6:
    sp.legend(bbox_to_anchor=(1.2,1))



plt.show()
