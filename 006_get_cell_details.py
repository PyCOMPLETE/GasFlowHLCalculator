from __future__ import division, print_function
import sys
import numpy as np
import matplotlib.pyplot as plt

import h5_storage
from compute_QBS_LHC import HeatLoadComputer

import LHCMeasurementTools.mystyle as ms

plt.close('all')

cells = ['13R4', '23L2']
filln = 5219

dict_cell_dict = {}
hl_dict = {}
data_ob = h5_storage.load_data_file(filln)

hlc = HeatLoadComputer(data_ob)

hl_cells_type = zip(np.mean(hlc.computed_values['qbs'],axis=0), hlc.dq.Cell_list, hlc.dq.Type_list)
hl_cells = []
for hl, cell, type_ in hl_cells_type:
    if type_ == 'ARC' and not np.isnan(hl):
        hl_cells.append((hl, cell))


hl_cells.sort()

cell_nrs = [0,1,-2,-1]
cells = [hl_cells[i][1] for i in cell_nrs]

for cell in cells:
    dict_cell_dict[cell] = hlc.get_single_cell_data(cell)

tt = (data_ob.timestamps - data_ob.timestamps[0])/3600

fig = ms.figure('Best and worst cell for fill %i' % filln)

sp = None
sp_t = plt.subplot(2,2,1, sharex=sp)
sp = sp_t
sp.set_title('Temperatures')
sp.set_xlabel('Time [h]')
sp.set_ylabel('T [K]')

sp_p = plt.subplot(2,2,2, sharex=sp)
sp = sp_p
sp.set_title('Pressures')
sp.set_xlabel('Time [h]')
sp.set_ylabel('Pressure [bar]')

sp_cv = plt.subplot(2,2,3, sharex=sp)
sp = sp_cv
sp.set_title('CV')
sp.set_xlabel('Time [h]')
sp.set_ylabel('CV [AU]')

sp_eh = sp_cv.twinx()
sp = sp_eh
sp.set_ylabel('EH [AU]')

sp_hl = plt.subplot(2,2,4, sharex=sp)
sp = sp_hl
sp.set_title('Heat loads')
sp.set_xlabel('Time [h]')

sp_mf = sp_hl.twinx()
sp = sp_mf
sp.set_ylabel('Mass flow')

for cell_ctr, cell in enumerate(cells):
    cell_dict = dict_cell_dict[cell]
    ls = ['-', '--', '-.', ':'][cell_ctr]
    for ctr, (key, data) in enumerate(cell_dict.iteritems()):
        if key[0] == 'T':
            sp = sp_t
        elif key[0] == 'P':
            sp = sp_p
        elif key == 'CV':
            sp = sp_cv
        elif key == 'qbs':
            sp = sp_hl
        elif key == 'EH':
            sp = sp_eh
        elif key == 'm_L':
            sp = sp_mf
        else:
            continue

        color = ms.colorprog(ctr, cell_dict)
        label = ' '.join([cell, key])

        sp.plot(tt, data, lw=2, label=label, ls=ls, color=color)

for sp in sp_t, sp_p:
    sp.legend(bbox_to_anchor=(1.1,1))

ms.comb_legend(sp_cv, sp_eh, bbox_to_anchor=(1.1,1))
ms.comb_legend(sp_hl, sp_mf, bbox_to_anchor=(1.1,1))
plt.show()

