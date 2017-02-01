from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt

import h5_storage
from compute_QBS_LHC import HeatLoadComputer

import LHCMeasurementTools.mystyle as ms

plt.close('all')

def plot_cell_details(cells, hlc, title):
    for cell in cells:
        dict_cell_dict[cell] = hlc.get_single_cell_data(cell)

    tt = (data_ob.timestamps - data_ob.timestamps[0])/3600

    fig = ms.figure(title)
    fig.subplots_adjust(right=0.86, wspace=0.28, left=0.04, top=0.93, bottom=0.08, hspace=0.23)

    sp = None
    sp_t = plt.subplot(2,2,1, sharex=sp)
    sp = sp_t
    sp.grid(True)
    sp.set_title('Temperatures')
    sp.set_xlabel('Time [h]')
    sp.set_ylabel('T [K]')

    sp_p = plt.subplot(2,2,2, sharex=sp)
    sp = sp_p
    sp.grid(True)
    sp.set_title('Pressures')
    sp.set_xlabel('Time [h]')
    sp.set_ylabel('Pressure [bar]')

    sp_cv = plt.subplot(2,2,3, sharex=sp)
    sp = sp_cv
    sp.grid(True)
    sp.set_title('CV')
    sp.set_xlabel('Time [h]')
    sp.set_ylabel('CV [AU]')

    sp_hl = plt.subplot(2,2,4, sharex=sp)
    sp = sp_hl
    sp.grid(True)
    sp.set_title('Heat loads')
    sp.set_xlabel('Time [h]')

    sp_mf = sp_hl.twinx()
    sp = sp_mf
    sp.set_ylabel('Mass flow')

    fig = ms.figure(title)
    fig.subplots_adjust(right=0.86, wspace=0.28, left=0.04, top=0.93, bottom=0.08, hspace=0.23)

    sp_h = plt.subplot(2,2,1, sharex=sp)
    sp = sp_h
    sp.grid(True)
    sp.set_ylabel('h3, hC')
    sp.set_xlabel('Time [h]')

    sp_r = sp.twinx()
    sp_r.set_ylabel('ro')

    for cell_ctr, cell in enumerate(cells):
        cell_dict = dict_cell_dict[cell]

        if len(cells) <= 4:
            ls = ['-', '--', '-.', ':'][cell_ctr]
        else:
            ls = '-'

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
                sp = sp_hl
            elif key == 'm_L':
                sp = sp_mf
            elif key in ('h3', 'hC'):
                sp = sp_h
            else:
                sp = sp_r

            color = ms.colorprog(ctr, cell_dict)
            if len(cells) <= 4 or cell_ctr == 0:
                label = ' '.join([cell, key])
            else:
                label=None

            sp.plot(tt, data, lw=2, label=label, color=color, ls=ls)

    bta = (1.2,1)

    sp_t.legend(bbox_to_anchor=bta)
    sp_p.legend(bbox_to_anchor=bta)
    sp_cv.legend(bbox_to_anchor=bta)
    ms.comb_legend(sp_hl, sp_mf, bbox_to_anchor=bta)
    ms.comb_legend(sp_h, sp_r, bbox_to_anchor=bta)


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

cell_nrs = [0, 200, -1]
cells = [hl_cells[i][1] for i in cell_nrs]
plot_cell_details(cells, hlc, 'Best, average and worst cell for fill %i' % filln)

cell_nrs = xrange(0,10)
cells = [hl_cells[i][1] for i in cell_nrs]
plot_cell_details(cells, hlc, 'Best cells for fill %i' % filln)

cell_nrs = xrange(len(hl_cells)-10, len(hl_cells))
cells = [hl_cells[i][1] for i in cell_nrs]
plot_cell_details(cells, hlc, 'Worst cells for fill %i' % filln)

plt.show()

