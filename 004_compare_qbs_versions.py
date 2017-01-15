# OUTDATED

from __future__ import division
import sys
import matplotlib.pyplot as plt
import numpy as np
import h5_storage
import qbs_fill as qf

import data_QBS_LHC as dql
import data_qbs as dqbs

dq = dqbs.data_qbs
this_module = sys.modules[__name__]
for key,value in dq.__dict__.iteritems():
    setattr(this_module, key, value)

if '..' not in sys.path: sys.path.append('..')
import LHCMeasurementTools.mystyle as ms
variable_list = [
#        'Cell_list',
#        'Type_list',
#        'Sector_list',
#        'EH84x_list',
#        'TT84x_list',
#        'CV94x_list',
#        'PT961_list',
#        'PT991_list',
#        'TT94x_list',
#        'TT961_list',
        'R_list',
        'Qs_list',
        'Kv_list',
        'nc_list',
        'L_list'
        ]

plt.close('all')
arc_list = ['ARC12','ARC23','ARC34','ARC45','ARC56','ARC67','ARC78','ARC81']
arc_nr_list = [arc[-2:] for arc in arc_list]

versions = [2, 3]
filln = 5219

fig = plt.figure()
title = 'Absolute'
fig.canvas.set_window_title(title)
fig.patch.set_facecolor('w')
fig.set_size_inches(15., 8.)

qbs_dict = {}
arc_avg_dict = {}
for ctr,version in enumerate(versions):
    sp = plt.subplot(len(versions), 1, ctr+1)
    sp.set_title('Version %i' % version)
    qbs_ob = qf.compute_qbs_fill(filln, version=version)
    qbs_dict[version] = qbs_ob
    arc_avg_dict[version] = qf.compute_qbs_arc_avg(qbs_ob)

    tt = (qbs_ob.timestamps - qbs_ob.timestamps[0])/3600.
    for cell in xrange(qbs_ob.data.shape[1]):
        sector_ctr = arc_nr_list.index(Sector_list[cell])
        sp.plot(tt, qbs_ob.data[:,cell], color=ms.colorprog(sector_ctr,8))

fig = plt.figure()
title = 'Delta'
fig.canvas.set_window_title(title)
fig.patch.set_facecolor('w')
fig.set_size_inches(15., 8.)

delta_data = qbs_dict[3].data - qbs_dict[2].data
delta_arc = arc_avg_dict[3] - arc_avg_dict[2]

mean_delta = np.mean(delta_data, axis=0)


sorted_delta = sorted(zip(mean_delta, Cell_list, Sector_list, Type_list), key=lambda a:a[0])

for ctr, i in enumerate(sorted_delta):
    if abs(i[0]) > 1:
        print(i)

print('')
for ctr, i in enumerate(sorted_delta):
    if i[2] == '78' and i[3] == 'ARC':
        print(i)

mask_arc = np.array(Type_list) == 'ARC'

sp = plt.subplot(2,1,1)
sp.set_title('Cells')
for cell in xrange(qbs_ob.data.shape[1]):
    sector_ctr = arc_nr_list.index(Sector_list[cell])
    if np.mean(delta_data[:,cell]) > 10:
        print('Max at cell', cell, Sector_list[cell], Cell_list[cell])
    if mask_arc[cell]:
        sp.plot(tt, delta_data[:,cell], color=ms.colorprog(sector_ctr,8))

sp = plt.subplot(2,1,2)
sp.set_title('Arcs')
for arc_ctr, arc in enumerate(arc_list):
    sp.plot(tt, delta_arc[:,arc_ctr], color=ms.colorprog(arc_ctr,8), label=arc[-2:])
sp.legend()

plt.show()
