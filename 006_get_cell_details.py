from __future__ import division, print_function
import os
import argparse
import cPickle as pickle

import numpy as np
import matplotlib.pyplot as plt

import h5_storage
from compute_QBS_LHC import HeatLoadComputer
import data_S45_details as dsd
import config_qbs as cq

import LHCMeasurementTools.mystyle as ms

plt.close('all')

parser = argparse.ArgumentParser()
parser.add_argument('--load', help='Load from pickle instead of calculating.', action='store_true')
parser.add_argument('--save', help='Save obj to pickle so you can load it later.', action='store_true')
parser.add_argument('--cells', help='Cells to plot', nargs='+')
parser.add_argument('--special', help='Special cells', action='store_true')
parser.add_argument('--best-worst', help='Best and worst hl cells', action='store_true')
parser.add_argument('--use_dP', action='store_true')
parser.add_argument('filln', type=int)
args = parser.parse_args()

filln = args.filln
use_dP = args.use_dP
cells = args.cells

if use_dP:
    storage = os.path.abspath(os.path.dirname(__file__))+ '/hlc_%i.pkl' % filln
else:
    storage = os.path.abspath(os.path.dirname(__file__))+ '/hlc_%i_nodP.pkl' % filln

load_pkl = args.load
save_pkl = args.save

new_cell = filln > 5700

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
    sp.set_ylabel('Heat load [W]')

    sp_mf = sp_hl.twinx()
    sp = sp_mf
    sp.set_ylabel('Mass flow [g/s]')

    fig = ms.figure(title)
    fig.subplots_adjust(right=0.86, wspace=0.28, left=0.08, top=0.93, bottom=0.08, hspace=0.23)

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
            yy_factor = 1
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
                yy_factor = 1e3
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

            sp.plot(tt, data*yy_factor, lw=2, label=label, color=color, ls=ls)

    bta = (1.2,1)

    sp_t.legend(bbox_to_anchor=bta)
    sp_p.legend(bbox_to_anchor=bta)
    sp_cv.legend(bbox_to_anchor=bta)
    ms.comb_legend(sp_hl, sp_mf, bbox_to_anchor=bta)
    ms.comb_legend(sp_h, sp_r, bbox_to_anchor=bta)


dict_cell_dict = {}
hl_dict = {}
data_ob = h5_storage.load_data_file(filln)

if load_pkl:
    with open(storage, 'r') as f:
        hlc = pickle.load(f)
else:
    hlc = HeatLoadComputer(data_ob, compute_Re=True, details=True, use_dP=use_dP)
    if save_pkl:
        with open(storage, 'w') as f:
            pickle.dump(hlc, f, -1)
            print('Pickle saved to %s' % storage)

hl_cells_type = zip(np.mean(hlc.computed_values['qbs'],axis=0), hlc.cq.Cell_list, hlc.cq.Type_list)
hl_cells = []
for hl, cell, type_ in hl_cells_type:
    if type_ == 'ARC' and not np.isnan(hl):
        hl_cells.append((hl, cell))
hl_cells.sort()

plot_cell_details(cells, hlc, 'Selected cells for fill %i' % filln)

if args.best_worst:
    cell_nrs = [0, 200, -1]
    cells = [hl_cells[i][1] for i in cell_nrs]
    plot_cell_details(cells, hlc, 'Best, average and worst cell for fill %i' % filln)

    cell_nrs = xrange(0,10)
    cells = [hl_cells[i][1] for i in cell_nrs]
    plot_cell_details(cells, hlc, 'Best cells for fill %i' % filln)

    cell_nrs = xrange(len(hl_cells)-10, len(hl_cells))
    cells = [hl_cells[i][1] for i in cell_nrs]
    plot_cell_details(cells, hlc, 'Worst cells for fill %i' % filln)

if args.special:
    cell_special_0 = dsd.cell_list
    if not new_cell:
        cell_special_0 = cell_special_0[:-1]

    cells_special = []
    for cell_ctr, cell in enumerate(cell_special_0):
        eh = dsd.EH84x_list[cell_ctr]
        index = cq.config_qbs.EH84x_list.index(eh)
        cells_special.append(cq.config_qbs.Cell_list[index])

    plot_cell_details(cells_special, hlc, 'Special cells for fill %i' % filln)

plt.show()

