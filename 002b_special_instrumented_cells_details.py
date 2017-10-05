from __future__ import division
import sys
import time
import re
import argparse
import cPickle as pickle

import numpy as np
import matplotlib.pyplot as plt

if '..' not in sys.path:
    sys.path.append('..')

import LHCMeasurementTools.TimberManager as tm
import LHCMeasurementTools.mystyle as ms
import LHCMeasurementTools.LHC_Heatloads as HL
from LHCMeasurementTools.SetOfHomogeneousVariables import SetOfHomogeneousNumericVariables
import LHCMeasurementTools.LHC_Energy as Energy
import LHCMeasurementTools.savefig as sf

import HeatLoadCalculators.impedance_heatload as ihl
import HeatLoadCalculators.synchrotron_radiation_heatload as srhl
import HeatLoadCalculators.FillCalculator as fc

from data_folders import data_folder_list

import GasFlowHLCalculator.qbs_fill as qf
from GasFlowHLCalculator.config_qbs import config_qbs
import GasFlowHLCalculator.data_S45_details as dsd
Cell_list = config_qbs.Cell_list

# Config
avg_pm_hrs = 0.1
fills_bmodes_name = './fills_and_bmodes.pkl'
separate_fig_for_hist = True


parser = argparse.ArgumentParser()
parser.add_argument('filln', type=int)
parser.add_argument('-a', help='Point in time where to calculate the heat load', type=float, default=-1.)
parser.add_argument('-d', help='Plot for all fills in 2015/2016', action='store_true')
parser.add_argument('-w', help='Histogram bin width', default=20., type=float)
parser.add_argument('--nolog', help='Do not show logged data', action='store_true')
parser.add_argument('--hist', help='Show histograms', action='store_true')
parser.add_argument('--nodetails', help='Show details of input data', action='store_true')
parser.add_argument('--pdsave', help='Save fig in pdijksta folder', action='store_true')
#parser.add_argument('--simple-titles', help='No confusing plot titles', action='store_true')
parser.add_argument('--nobroken', help='Do not plot the cell with broken sensor', action='store_true')
parser.add_argument('--separate', help='Show heat load for BS separately', action='store_true')
parser.add_argument('--noshow', help='Do not call plt.show', action='store_true')
parser.add_argument('--contributions', help='Show contributions', action='store_true')
parser.add_argument('--no-use-dP', help='Load cell data without dP', action='store_true')
parser.add_argument('--force-recalc', help='Recalc special qbs instead of loading from h5.', action='store_true')
parser.add_argument('--crude-cell-swap', help='Swap cell /timber var relation..', action='store_true')
parser.add_argument('--xlim', help='Set xlim for plots in hours', nargs=2, type=float, default=(None, None))
args = parser.parse_args()

filln = args.filln
avg_time_hrs = args.a
show_dict = args.d
binwidth = args.w
logged = not args.nolog
hist = args.hist
details = not args.nodetails
use_dP = not args.no_use_dP
use_new_cell = filln > 5700
swap_cell_for_old_fills = args.crude_cell_swap
force_xlim = args.xlim

myfontsz = 12
ms.mystyle_arial(fontsz=myfontsz, dist_tick_lab=8)
re_dev = re.compile('^QRLA[AB]_(\d\d[RL]\d)_QBS\d{3}_([QD]\d).POSST$')
plt.close('all')

# Definitions
if args.nobroken:
    cells = ['13L5', '13R4']
else:
    cells = ['13L5', '33L5', '13R4']
new_cell = '31L2'

if use_new_cell:
    cells.append(new_cell)

cell_notation_recalc_logged = {
    '33L5': '33L5',
    '13R4': '13L5',
    '13L5': '13R5',
    new_cell: new_cell,
}

#if args.simple_titles:
cell_title_dict = {
    '13L5': '13L5',
    '33L5': '33L5',
    '13R4': '13R4',
    '31L2': '31L2',
    }
#else:
#    cell_title_dict = {
#            '13L5': '13L5 / 12R4',
#            '33L5': '33L5 / 32R4 (broken sensor)',
#            '13R4': '13R4 / 13L5 (reversed gas flow)',
#            '31L2': '31L2 / 32L2 (new cell)',
#            }
cells_and_new = cells + [new_cell]
affix_list = ['Q1', 'D2', 'D3', 'D4']
beam_colors = {1: 'b', 2: 'r'}

# Which cells are the special ones in the qbs data?
cell_index_dict = {}
for cell_ctr, cell in enumerate(dsd.cell_list):
    eh = dsd.EH84x_list[cell_ctr]
    index = config_qbs.EH84x_list.index(eh)
    cell_index_dict[cell] = index

keys = ['special_HC_Q1', 'special_HC_D2', 'special_HC_D3', 'special_HC_D4', 'special_total']
variable_list = HL.variable_lists_heatloads['special_total'][:]
for list_ in dsd.QBS_12R4, dsd.QBS_13L5, dsd.QBS_32R4:
    variable_list.extend(list_)

if new_cell not in cells:
    variable_list = filter(lambda x: not(new_cell in x), variable_list)
else:
    variable_list.extend(dsd.QBS_32L2)


# Dictionary for variables
hl_dict_logged = {
    '13R4': {
        'Q1': 'QRLAA_13R4_QBS947_Q1.POSST',
        'D2': 'QRLAA_13R4_QBS947_D2.POSST',
        'D3': 'QRLAA_13R4_QBS947_D3.POSST',
        'D4': 'QRLAA_13R4_QBS947_D4.POSST',
        'Cell': 'QRLAA_13R4_QBS947.POSST',
    },
    '33L5': {
        'Q1': 'QRLAA_33L5_QBS947_Q1.POSST',
        'D2': 'QRLAA_33L5_QBS947_D2.POSST',
        'D3': 'QRLAA_33L5_QBS947_D3.POSST',
        'D4': 'QRLAA_33L5_QBS947_D4.POSST',
        'Cell': 'QRLAA_33L5_QBS947.POSST',
    },
    '13L5': {
        'Q1': 'QRLAA_13L5_QBS943_Q1.POSST',
        'D2': 'QRLAA_13L5_QBS943_D2.POSST',
        'D3': 'QRLAA_13L5_QBS943_D3.POSST',
        'D4': 'QRLAA_13L5_QBS943_D4.POSST',
        'Cell': 'QRLAA_13L5_QBS943.POSST',
    },
    '31L2': {
        'Q1': 'QRLAB_31L2_QBS943_Q1.POSST',
        'D2': 'QRLAB_31L2_QBS943_D2.POSST',
        'D3': 'QRLAB_31L2_QBS943_D3.POSST',
        'D4': 'QRLAB_31L2_QBS943_D4.POSST',
        'Cell': 'QRLAB_31L2_QBS943.POSST',
    },
}

if swap_cell_for_old_fills and filln < 5700:
    hl_dict_logged = {
        '13R4': { # swapped magnets with cell
            'Q1': 'QRLAA_13L5_QBS943_Q1.POSST',
            'D2': 'QRLAA_13L5_QBS943_D2.POSST',
            'D3': 'QRLAA_13L5_QBS943_D3.POSST',
            'D4': 'QRLAA_13L5_QBS943_D4.POSST',
            'Cell': 'QRLAA_13R4_QBS947.POSST',
        },
        '33L5': {
            'Q1': 'QRLAA_33L5_QBS947_Q1.POSST',
            'D2': 'QRLAA_33L5_QBS947_D2.POSST',
            'D3': 'QRLAA_33L5_QBS947_D3.POSST',
            'D4': 'QRLAA_33L5_QBS947_D4.POSST',
            'Cell': 'QRLAA_33L5_QBS947.POSST',
        },
        '13L5': { # swapped magnets with cell
            'Q1': 'QRLAA_13R4_QBS947_Q1.POSST',
            'D4': 'QRLAA_13R4_QBS947_D2.POSST',  # swapped D2-D4
            'D3': 'QRLAA_13R4_QBS947_D3.POSST',
            'D2': 'QRLAA_13R4_QBS947_D4.POSST',  # swapped D2-D4
            'Cell': 'QRLAA_13L5_QBS943.POSST',
        },
        '31L2': {
            'Q1': 'QRLAB_31L2_QBS943_Q1.POSST',
            'D2': 'QRLAB_31L2_QBS943_D2.POSST',
            'D3': 'QRLAB_31L2_QBS943_D3.POSST',
            'D4': 'QRLAB_31L2_QBS943_D4.POSST',
            'Cell': 'QRLAB_31L2_QBS943.POSST',
        },
    }

# merge pickles and add info on location
dict_fill_bmodes={}
for df in data_folder_list:
    with open(df+'/fills_and_bmodes.pkl', 'rb') as fid:
        this_dict_fill_bmodes = pickle.load(fid)
        for kk in this_dict_fill_bmodes:
            this_dict_fill_bmodes[kk]['data_folder'] = df
        dict_fill_bmodes.update(this_dict_fill_bmodes)

if avg_time_hrs == -1.:
    if dict_fill_bmodes[filln]['t_start_STABLE'] != -1:
        avg_time_hrs = (dict_fill_bmodes[filln]['t_start_STABLE'] - dict_fill_bmodes[filln]['t_startfill'])/3600.
    else:
        print('Warning: Avg time hrs = 0.5')
        avg_time_hrs = 0.5


# get location of current data
data_folder_fill = dict_fill_bmodes[filln]['data_folder']

fill_dict = {}
if filln < 4857:
    fill_dict.update(tm.parse_timber_file('/afs/cern.ch/project/spsecloud/LHC_2015_PhysicsAfterTS2/fill_csvs/fill_%d.csv' % filln, verbose=False))
    fill_dict.update(tm.parse_timber_file('/afs/cern.ch/project/spsecloud/LHC_2015_PhysicsAfterTS2/heatloads_fill_h5s/heatloads_all_fill_%i.h5' % filln, verbose=False))
else:
    fill_dict.update(tm.parse_timber_file(data_folder_fill + '/fill_basic_data_csvs/basic_data_fill_%d.csv' % filln, verbose=False))
    fill_dict.update(tm.parse_timber_file(data_folder_fill + '/fill_heatload_data_csvs/heatloads_fill_%d.csv' % filln, verbose=False))

energy = Energy.energy(fill_dict, beam=1)
energy.t_stamps = (energy.t_stamps - energy.t_stamps[0])/3600.
t_ref = dict_fill_bmodes[filln]['t_startfill']
tref_string = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime(t_ref))

heatloads = SetOfHomogeneousNumericVariables(variable_list=variable_list, timber_variables=fill_dict).aligned_object(dt_seconds=60)
## Swap 13L5 D2 and D4, as the correct gas flow is in reverse direction
#for ctr,var in enumerate(heatloads.variables):
#    if '13R4' in var:
#        if 'D2' in var:
#            ctr_2 = ctr
#            var_2 = var
#        elif 'D4' in var:
#            ctr_4 = ctr
#            var_4 = var
#heatloads.variables[ctr_2] = var_4
#heatloads.variables[ctr_4] = var_2

y_min, y_max = -10, np.max(heatloads.data)+5
timestamps = (heatloads.timestamps - heatloads.timestamps[0])/3600.


# Recalculated objects
special_hl = qf.special_qbs_fill(filln, force_recalc=args.force_recalc)
special_tt = (special_hl['timestamps'] - special_hl['timestamps'][0]) / 3600.
qbs_ob = qf.compute_qbs_fill(filln, use_dP=use_dP)
qbs_tt = (qbs_ob.timestamps - qbs_ob.timestamps[0])/3600.
atd_mask_mean = np.abs(qbs_tt - avg_time_hrs) < avg_pm_hrs

lhc_hist_dict = qf.lhc_histograms(qbs_ob, avg_time_hrs, avg_pm_hrs)
arc_hist_total = lhc_hist_dict['total']
arc_hist_dict = lhc_hist_dict['arcs']

# Plots

title = 'Special instrumented cells for fill %i' % filln
fig = ms.figure(title)
fig.subplots_adjust(left=.1, right=.75, hspace=.38, wspace=.45)
sp = None

# Cells
var_mean_dict = {} # unused
for cell_ctr, cell in enumerate(cells):
    cell_vars = hl_dict_logged[cell]
    sp_ctr = cell_ctr +1
    sp = plt.subplot(2,2,sp_ctr, sharex=sp)
    sp.set_title(cell_title_dict[cell])
    sp.grid(True)
    sp.set_xlabel('Time [h]')
    sp.set_ylabel('Heat load [W]')
#    if filln == 5277 and cell == '13L5':
#        sp.set_ylim(-150, y_max)
    sp2 = sp.twinx()
    sp2.set_ylim(0,7)
    sp2.set_ylabel('Energy [TeV]')
    sp2.plot(energy.t_stamps, energy.energy/1e3, c='black', lw=2., label='Energy')

    summed_log, summed_re = 0., 0.
    for ctr, affix in enumerate(affix_list):
        var = cell_vars[affix]
        row_ctr = heatloads.variables.index(var)
        values = heatloads.data[:,row_ctr]
        summed_log += values
        summed_re += special_hl[cell][affix]['Sum']
        color = ms.colorprog(ctr, cell_vars)
        if logged:
            sp.plot(timestamps, values, label=affix+' logged', ls='--', lw=2., color=color)
        sp.plot(special_tt, special_hl[cell][affix]['Sum'], ls='-', lw=2., color=color, label=affix)
    #sp.axvline(avg_time_hrs, color='black')
    # summed
    if logged:
        sp.plot(timestamps, summed_log, ls='--', lw=2., color='blue', label='Sum logged')
    sp.plot(special_tt, summed_re, ls='-', lw=2., color='blue', label='Sum')

    # celss
    sp.plot(qbs_tt, qbs_ob.data[:,cell_index_dict[cell]], label='Cell recalc.', ls='-', lw=2., c='orange')
    cell_index = heatloads.variables.index(cell_vars['Cell'])
    if logged:
        sp.plot(timestamps, heatloads.data[:,cell_index], label='Cell logged', ls='--', lw=2., c='orange')
    sp.set_ylim(-20, None)
    if sp_ctr == 2:
        ms.comb_legend(sp,sp2,bbox_to_anchor=(1.8,1), fontsize=myfontsz)
    sp.set_xlim(force_xlim)

# Also show LHC hist
def round_to(arr, precision):
    return np.round(arr/precision)*precision

bins = np.arange(round_to(arc_hist_total.min(),binwidth)-binwidth, round_to(arc_hist_total.max(),binwidth)+binwidth*3/2, binwidth)

if separate_fig_for_hist:
    fig = ms.figure('Special cells in histogram of heatloads')
    sp = plt.subplot(2,2,1)
else:
    sp = plt.subplot(2,2,4)
sp.set_title('LHC cell heat load - all arcs')
sp.grid(True)
sp.set_xlabel('Time [h]')
sp.set_ylabel('Heat load [W]')
sp.hist(arc_hist_total, bins=bins, alpha=0.5, color='blue')
colors=['red', 'green', 'orange', 'black', 'blue']
for cell_ctr, cell in enumerate(cells_and_new):
    mean = np.mean(qbs_ob.data[atd_mask_mean,cell_index_dict[cell]])
    sp.axvline(mean, label=cell, color=colors[cell_ctr], lw=2)
sp.legend(bbox_to_anchor=(1.6,1), title='Recalc. cells')

# Histograms
if hist:
    # 1 for each arc
    #bins = np.arange(-50, 251, 300./11.)
    for ctr, (arc, data) in enumerate(arc_hist_dict.iteritems()):
        sp_ctr = ctr % 4 + 1
        if sp_ctr == 1:
            title = 'Fill %i: Heat loads at %.1f hours' % (filln, avg_time_hrs)
            fig = ms.figure(title)
        sp = plt.subplot(2,2,sp_ctr)
        sp.hist(arc_hist_total, bins=bins, alpha=0.5, color='blue', weights=1./len(arc_hist_total)*np.ones_like(arc_hist_total))
        sp.hist(data, bins=bins, color='green', alpha=0.5, weights=1./len(data)*np.ones_like(data))
        sp.axvline(np.mean(data), lw=2., color='green')
        sp.axvline(np.mean(arc_hist_total), lw=2., color='blue')
        sp.grid(True)
        sp.set_xlabel('Heat load [W]')
        sp.set_ylabel('# Half cells (normalized)')
        sp.set_title('Arc %s' % arc)

        if arc == '45':
            colors = ['red', 'orange', 'brown']
            for cell_ctr, cell in enumerate(cells):
                mean = np.mean(qbs_ob.data[atd_mask_mean,cell_index_dict[cell]])
                sp.axvline(mean, label=cell, color=colors[cell_ctr])
            sp.legend(bbox_to_anchor=(1.2,1))
        elif arc == '12':
            mean = np.mean(qbs_ob.data[atd_mask_mean,cell_index_dict[new_cell]])
            sp.axvline(mean, label=new_cell, color='red')
            sp.legend(bbox_to_anchor=(1.2,1))

#    # 1 plot for all sectors
#    title = 'Fill %i at %.1f h: LHC Arcs histograms' % (filln, avg_time_hrs)
#    fig = plt.figure()
#    fig.canvas.set_window_title(title)
#    fig.patch.set_facecolor('w')
#    sp_hist = plt.subplot(2,2,1)
#    sp_hist.set_xlabel('Heat load [W]')
#    sp_hist.set_ylabel('# Half cells')
#    sp_hist.set_title(title)
#    sp_hist.grid(True)
#    for ctr, (arc, data) in zip(xrange(len(arc_hist_dict)), arc_hist_dict.iteritems()):
#        hist, null = np.histogram(data, bins=bins)
#        sp_hist.step(bins[:-1]+10, hist, label='Arc %s' % arc, color=ms.colorprog(ctr, arc_hist_dict), lw=2)
#        #sp_hist.hist(data, bins=bins, color=ms.colorprog(ctr, arc_hist_dict), alpha=0.2, lw=2., label='Arc %s' % arc)
#
#    #hist, null = np.histogram(arc_hist_total, bins=bins)
#    #sp_hist.plot(bins[:-1]+10, hist,'.', label='All', markersize=3.)
#    sp_hist.legend(bbox_to_anchor=(1.2,1))


# Compare dipoles to quads
fig_dev = ms.figure('Compare devices')
fig_dev.subplots_adjust(left=.06, right=.84, top=.93, hspace=.38, wspace=.46)

# Logged data
sp_dip = plt.subplot(2,2,1)
sp_quad = plt.subplot(2,2,3)
dip_list = []
quad_list = []

ls_list = ['-', '--', '-.', ':']

for var in variable_list:
    for affix in affix_list:
        if affix in var:
            if 'D' in affix:
                dip_list.append(var)
            elif 'Q' in affix:
                quad_list.append(var)

for sp, dev_list, title in zip((sp_dip, sp_quad), (dip_list, quad_list), ('Dipoles', 'Quadrupoles')):
    dev_list.sort()
    for ctr, dev in enumerate(dev_list):
        index = heatloads.variables.index(dev)
        values = heatloads.data[:,index]
        info = re_dev.search(dev).groups()
        if args.nobroken and info[0] == '33L5':
            continue
        if info == ('33L5', 'D4') or info == ('33L5', 'D3'):
            continue
        color = ms.colorprog(affix_list.index(info[1]), affix_list)
        ls = ls_list[cells.index(info[0])]
        label = ' '.join(info)
        sp.plot(timestamps, values, label=label, lw=2., color=color, ls=ls)
    sp.set_title(title+' logged')
    sp.legend(bbox_to_anchor=(1.2,1))
    sp.set_ylim(-20,None)
    sp.set_ylabel('Heat load [W]')
    sp.set_xlabel('Time [h]')
    sp.grid(True)

# Recalculated data

sp_dip = plt.subplot(2,2,2)
sp_quad = plt.subplot(2,2,4)
for cell_ctr, cell in enumerate(cells):
    ls = ls_list[cell_ctr]
    for a_ctr, affix in enumerate(affix_list):
        if cell == '33L5' and affix in ['D4', 'D3']:
            continue
        color = ms.colorprog(a_ctr, affix_list)
        if affix == 'Q1':
            sp_quad.plot(special_tt, special_hl[cell][affix]['Sum'], color=color, ls=ls, label=cell, lw=2)
        else:
            sp_dip.plot(special_tt, special_hl[cell][affix]['Sum'], color=color, ls=ls, label=cell+ ' ' + affix, lw=2)

for sp, title in zip((sp_dip, sp_quad), ('Dipoles', 'Quadrupoles')):
    sp.grid(True)
    sp.set_title(title+' recalculated')
    sp.set_ylabel('Heat load [W]')
    sp.set_xlabel('Time [h]')
    sp.set_ylim(-20,None)
    sp.legend(bbox_to_anchor=(1.2,1))



# From large HL dict
if show_dict:
    from hl_dicts.LHC_Heat_load_dict import hl_dict, mask_dict, main_dict_2016

    mask = hl_dict['stable_beams']['n_bunches']['b1'] > 1000
    ff_2016 = main_dict_2016['filln'][0]
    hl_dict = mask_dict(hl_dict, mask)

    int_norm_factor = hl_dict['stable_beams']['intensity']['total']
    ylims = [(-.1e-13,4.5e-13), (-10, 120), (-0.5,5), (-0.5,9)]
    norm_factors = [int_norm_factor, 1., 0., 0.]
    titles = ['Normalized by intensity', 'Heat loads', 'Normalized to 1', 'Normalized by intensity and to 1']
    units = ['[W/p]', '[W]', 'arb. units', 'arb. units']

    for ctr, norm_factor, ylim, title, unit in zip(xrange(len(norm_factors)), norm_factors, ylims, titles, units):
        sp = None
        for cell_ctr, cell in enumerate(cells):
            sp_ctr = cell_ctr % 3 + 1
            if sp_ctr == 1:
                fig = ms.figure('HL dict' + title)
                fig.subplots_adjust(left=.06, right=.84, top=.93, hspace=.38, wspace=.42)
            sp = plt.subplot(3,1,sp_ctr, sharex=sp)
            sp.set_xlabel('Fill number')
            sp.set_ylabel('Heat load %s' % unit)
            sp.set_title(cell_title_dict[cell])
            sp.grid(True)
            names = [cell]
            labels = [cell]
            for var in variable_list:
                if cell in var:
                    info = re_dev.search(var)
                    if info != None:
                        names.append(''.join(re_dev.search(var).groups()))
                        labels.append(info.group(2).replace('_',' '))
                    else:
                        print('None: ', cell, var)

            for name, label in zip(names, labels):
                if ctr == 2:
                    norm_factor_2 = np.mean(hl_dict['stable_beams']['heat_load'][name][-20:])
                elif ctr == 3:
                    norm_factor_2 = np.mean(hl_dict['stable_beams']['heat_load'][name][-20:]/int_norm_factor[-20:])*int_norm_factor
                else:
                    norm_factor_2 = norm_factor
                sp.plot(hl_dict['filln'], hl_dict['stable_beams']['heat_load'][name]/norm_factor_2,'.', lw=2., label=label, markersize=12)

            sp.axvline(ff_2016, color='black', lw=2.)
            sp.legend(bbox_to_anchor=(1.05,1))
            sp.set_ylim(*ylim)
            sp.grid(True)

# Show details on input data variables
if details:
    import GasFlowHLCalculator.data_S45_details as dsd
    import GasFlowHLCalculator.h5_storage as h5_storage
    import GasFlowHLCalculator.compute_QBS_special as cqs

    #special_atd = mfm.h5_to_obj(h5_dir + 'special_cells/special_data_fill_%i.h5' % filln)
    special_atd = h5_storage.load_special_data_file(filln)
    special_tt = (special_atd.timestamps - special_atd.timestamps[0])/3600.
    special_vars = list(special_atd.variables)
    # Mask first and last hour
    mask = np.logical_and(special_tt > 0, special_tt < special_tt[-1] -0)
    special_tt = special_tt[mask]

    #alternate_notation = {
    #        '13L5': '12R4',
    #        '33L5': '32R4',
    #        '13R4': '13L5',
    #        '31L2': '32L2',
    #        }
    alternate_notation = {
            '13L5': '13L5',
            '33L5': '32R4',
            '13R4': '12R4',
            '31L2': '32L2',
            }

    re_var = re.compile('^\w{4}_\w?(\d\d[RL]\d_TT\d{3})\.POSST$')
    title = 'Fill %i temperature sensors' % filln
    fig = ms.figure(title)
    fig.subplots_adjust(left=.06, right=.84, top=.93, hspace=.38, wspace=.42)
    for cell_ctr, cell in enumerate(cells):
        acell = alternate_notation[cell]
        lists = {}
        for key, item in dsd.__dict__.iteritems():
            if type(item) is list and acell in key:
                if ('Tin' in key or 'Tout' in key):
                    lists[key] = item

        sp = plt.subplot(2,2,cell_ctr+1)
        sp.set_ylabel('Temperature [K]')
        sp.set_xlabel('Time [h]')
        sp.grid(True)
        sp.set_title(cell_title_dict[cell])
        used_vars = set()
        # Temperatures
        for cc, name in enumerate(sorted(lists.keys())):
            ll = lists[name]
            color = ms.colorprog(cc, lists)
            for ctr, var in enumerate(ll):
                if var not in used_vars:
                    used_vars.add(var)
                    if ctr == 0:
                        label = name
                    else:
                        label = None
                    ls = ls_list[ctr]
                    row_ctr = special_vars.index(var)
                    sp.plot(special_tt, special_atd.data[mask,row_ctr], lw=2, label=label, ls=ls, color=color)
        sp.legend(bbox_to_anchor=(1.3,1))

    # Separate hl for b1, b2
    if args.separate:
        title = 'Fill %i separate beam screens' % filln
        fig = ms.figure(title)
        fig.subplots_adjust(left=.06, right=.84, top=.93, hspace=.38, wspace=.42)

        qbs_special = cqs.compute_qbs_special(special_atd, separate=True)
        qbs_tt = (qbs_special['timestamps']-qbs_special['timestamps'][0])/3600.
        for cell_ctr, cell in enumerate(cells):
            sp_ctr = cell_ctr+1
            sp = plt.subplot(2,2,sp_ctr, sharex=sp)
            sp.set_title(cell_title_dict[cell])
            sp.grid(True)
            sp.set_xlabel('Time [h]')
            sp.set_ylabel('Heat load [W]')
            sp2 = sp.twinx()
            sp2.set_ylim(0,7)
            sp2.set_ylabel('Energy [TeV]')
            sp2.plot(energy.t_stamps, energy.energy/1e3, c='black', lw=2., label='Energy')

            for actr,affix in enumerate(affix_list[1:]):
                for bctr, beam in enumerate(('b1','b2')):
                    key = affix+'_'+beam
                    if key in qbs_special[cell]:
                        color = ms.colorprog(actr, affix_list)
                        ls = ['-','--'][bctr]
                        sp.plot(qbs_tt, qbs_special[cell][key], ls=ls, lw=2., label=key, color=color)
            sp.set_ylim(-10, None)
            ms.comb_legend(sp,sp2,bbox_to_anchor=(1.3,1), fontsize=myfontsz)

if args.contributions:
    fig = ms.figure('Heat load contributions')
    fig.subplots_adjust(right=0.75, wspace=0.45, hspace=0.4)
    sp = None

    fill_dict.update(tm.parse_timber_file('./fill_bunchbybunch_data_csvs/bunchbybunch_data_fill_%d.csv' % filln, verbose=False))

    hli_calculator  = ihl.HeatLoadCalculatorImpedanceLHCArc()
    hlsr_calculator  = srhl.HeatLoadCalculatorSynchrotronRadiationLHCArc()

    hl_imped_fill = fc.HeatLoad_calculated_fill(fill_dict, hli_calculator)
    hl_sr_fill = fc.HeatLoad_calculated_fill(fill_dict, hlsr_calculator)

    model_time = (hl_sr_fill.t_stamps-t_ref)/3600

    t_start_injphys = (dict_fill_bmodes[filln]['t_start_INJPHYS'] -t_ref)/3600.
    def compute_offset(tt, yy):
        mask_offset = np.logical_and(tt < t_start_injphys, tt > t_start_injphys - 600)
        offset = np.mean(yy[mask_offset])
        return offset



    for cell_ctr, cell in enumerate(cells):
        sp_ctr = cell_ctr % 4 + 1
        index = cell_index_dict[cell]

        sp = plt.subplot(2,2,sp_ctr, sharex=sp)
        sp.set_xlabel('Time [h]')
        sp.set_ylabel('Heat load [W/hc]')
        sp.set_title(cell)
        sp.grid(True)

        sp2 = sp.twinx()
        sp2.set_ylabel('Energy [TeV]')
        sp2.plot(energy.t_stamps, energy.energy/1e3, c='black', lw=2., ls='--', label='Energy')

        # Total
        xx_time = special_tt
        #yy = qbs_ob.data[:,index]
        #yy_all = yy - compute_offset(xx_time, yy)
        #sp.plot(xx_time, yy_all, label='Total', color='black', lw=3)

        # Imp
        bottom = np.zeros_like(xx_time)
        data = np.interp(xx_time, model_time, hl_imped_fill.heat_load_calculated_total)*53.45
        offset = compute_offset(xx_time, data)
        top = data-offset
        sp.fill_between(xx_time, bottom, top, label='Impedance', alpha=0.5, color='r')
        bottom = np.copy(top)

        # SR
        data = np.interp(xx_time, model_time, hl_sr_fill.heat_load_calculated_total)*53.45
        offset = compute_offset(xx_time, data)
        top += data-offset
        sp.fill_between(xx_time, bottom, top, label='SR', alpha=0.5, color='b')
        bottom = np.copy(top)

        # Quad
        data = special_hl[cell]['Q1']['Sum']
        offset = compute_offset(xx_time, data)
        top += data-offset
        sp.fill_between(xx_time, bottom, top, label='Average quad', alpha=0.5, color='orange')
        bottom = np.copy(top)

        # Rest
        # Make sure the rest is not smaller than imp+sr+quad
        data = special_hl[cell]['Sum'] - special_hl[cell]['Q1']['Sum']
        offset = compute_offset(xx_time, data)
        top = data-offset
        sp.fill_between(xx_time, bottom, top, label='Dipoles and drifts', alpha=0.5, color='green')
        bottom = np.copy(top)

        #Unexplained
        #sp.fill_between(xx_time, bottom, yy_all, label='Unexplained', alpha=0.5, color='brown')

        sp.set_ylim(0, 120)

    sp.legend(loc='upper left', bbox_to_anchor=(1.12,1))





if args.pdsave:
    sf.saveall_pdijksta()

if not args.noshow:
    plt.show()
