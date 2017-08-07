from __future__ import division
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
import cPickle as pickle
import argparse

import Helium_properties as hp
from Helium_properties import interp_P_T_hPT, interp_P_T_DPT, interp_P_T_mu
import compute_QBS_LHC as qbl
import h5_storage

import LHCMeasurementTools.mystyle as ms
#import LHCMeasurementTools.savefig as sf

parser = argparse.ArgumentParser()
parser.add_argument('--calc', help='Calculate instead of loading from pickle.', action='store_true')
parser.add_argument('--savefig', help='Save in file.')
parser.add_argument('--noshow', help='Do not call plt.show', action='store_true')
parser.add_argument('--onlyarcs', help='Only show arc_cells', action='store_true')
parser.add_argument('--onlyinterp', help='Only show interp values', action='store_true')
args = parser.parse_args()

recompute = args.calc
only_interp = args.onlyinterp
figs = []


ms.mystyle()

plt.close('all')

filln = 5219
pkl_file = os.path.abspath(os.path.dirname(__file__)) + '/hlc_%i.pkl' % filln
interp_P_T_gamma = interp2d(hp.P, hp.T, hp.gamma_PT)
#t_arr = np.arange(4,30,0.5)
mask_t_arr = np.logical_and(hp.T <= 30, hp.T >= 4)
t_arr = hp.T[mask_t_arr]
mask_p_arr = np.logical_and(hp.P <= 5, hp.P >= 2)
p_arr = hp.P[mask_p_arr]

with open(pkl_file) as f:
    hlc = pickle.load(f)

gamma = np.zeros_like(hlc.computed_values['m_L'])
P3 = hlc.computed_values['P3']
T3 = hlc.data_dict['T3']
for j in xrange(hlc.Nvalue):
    for i in xrange(hlc.Ncell):
        gamma[j,i] = interp_P_T_gamma(P3[j,i], T3[j,i])

fig = ms.figure('Interpolated Data', figs)

interps = (interp_P_T_hPT, interp_P_T_DPT, interp_P_T_mu, interp_P_T_gamma)
titles  = ('Enthalpy', 'Density', 'Viscosity', 'Heat capacity ratio')
pressures = np.arange(2.5,5.1,0.5)

for ctr, (interp, title) in enumerate(zip(interps, titles)):
    sp_ctr = ctr+1

    sp = plt.subplot(2,2,sp_ctr)
    sp.grid(True)
    sp.set_title(title)
    sp.set_xlabel('Temperature [K]')
    sp.set_ylabel(title)
    if title in ('Enthalpy', 'Viscosity'):
        ms.sciy()
    for pressure in p_arr:
        values = np.empty_like(t_arr)
        for i, temp in enumerate(t_arr):
            values[i] = interp(pressure, temp)
        sp.plot(t_arr, values, lw=2., label='P=%.1f bar' % pressure, marker='o')

    if interp is interp_P_T_gamma:
        sp.axhline(5./3., lw=2, ls='--', label='5/3')

    #if sp_ctr == 4:
    sp.legend(loc=0)

if only_interp:
    plt.show()
    sys.exit()

if recompute:
    atd = h5_storage.load_data_file(filln)
    hlc = qbl.HeatLoadComputer(atd, compute_Re=True)
else:
    with open(pkl_file) as f:
        hlc = pickle.load(f)

combined_dict = hlc.data_dict.copy()
combined_dict.update(hlc.computed_values)
combined_dict['Ratio $P_4$/$P_3$'] = combined_dict['P4'] / combined_dict['P3']
combined_dict['gamma'] = gamma
combined_dict['$\Delta P$ / P1'] = (combined_dict['P3'] - combined_dict['P1'])/combined_dict['P1']

del combined_dict['EH']

tt = hlc.atd_ob.timestamps
index_tt = np.argmin(np.abs(tt - tt[0] - 3600*2))

fig_ctr = 1
key_sp_dict = {}
for cc in (1,2):
    for ctr, (key, data) in enumerate(sorted(combined_dict.items())):

        if key.startswith('P'):
            unit = '$P_%s$ [bar]' % key[1]
        elif key.startswith('T'):
            unit = '$T_%s$ [K]' % key[1]
        elif key == 'CV':
            unit = 'u [%]'
        elif key == 'gamma':
            unit = '$\gamma$'
        elif key == 'qbs':
            unit = '$Q_{BS}$ [W]'
        elif key == 'm_L':
            unit = '$m_L [kg/s]$'
        elif key == 'h3':
            unit = '$h_3$ [J/kg]'
        elif key == 'hC':
            key = 'h1'
            unit = '$h_1$ [J/kg]'
        else:
            unit = key

        if cc == 1:
            data2_pure = np.nan_to_num(data)
            plot_title = ('Occurrence of raw data values for fill %i' % filln)
            label = "All data"
        else:
            data2_pure = np.nan_to_num(data[index_tt,:])
            plot_title = ('Occurrence of raw data values for fill %i after 2 hrs' % filln)
            label = "After 2 hours"

        if args.onlyarcs:
            plot_title += ' - Only Arcs'
            for type_ctr, type_ in enumerate(hlc.cq.Type_list):
                if type_ != 'ARC':
                    if cc == 1:
                        data2_pure[:,type_ctr] = 0
                    else:
                        data2_pure[type_ctr] = 0

            data2_pure = data2_pure.flatten()


        data2 = data2_pure[data2_pure != 0]
        hist, bin_edges = np.histogram(data2, bins=10, normed=True)
        new_hist = []
        new_bin_edges = []
        for h_ctr, val in enumerate(hist):
            edge1, edge2 = bin_edges[h_ctr], bin_edges[h_ctr+1]
            binwidth = edge2 - edge1
            if val*binwidth > 0.01:
                if not new_bin_edges:
                    new_bin_edges.append(edge1)
                new_bin_edges.append(edge2)
                new_hist.append(val)

        if key == 'EH':
            new_bin_edges = [0,20]
        data2 = data2[np.logical_and(data2 > new_bin_edges[0], data2 < new_bin_edges[-1])]

        #if key == 'T2':
        #    data2 = data2[data2 < 30]
        #elif key == 'gamma':
        #    data2 = data2[data2 < 3]
        #elif key == 'm_L':
        #    data2 = data2[data2 < 0.0025]
        #elif key == 'ro':
        #    data2 = data2[data2 < 60]

        if key in hlc.data_dict:
            affix = '(data)'
        else:
            affix = '(computed)'
        sp_ctr = ctr % 4 +1
        if cc == 1 and sp_ctr == 1:
            fig = ms.figure(plot_title, figs)
            fig.subplots_adjust(right=0.75, )
            fig_ctr += 1
        if cc == 1:
            sp = plt.subplot(2,2,sp_ctr)
            sp.set_title(key+' '+affix)
            sp.set_xlabel(unit)
            sp.set_ylabel('Occurence')
            #sp.set_yticklabels([])
            if key in ('h3', 'm_L', 'h1', 'Re', 'h3'):
                ms.scix()
            ms.sciy()
            sp.grid(True)
            key_sp_dict[key] = sp
        else:
            sp = key_sp_dict[key]

        sp.hist(data2, normed=True, label=label, alpha=0.5)
        if key[0] == '$':
            sp.set_xticks(np.arange(-0.12, 0.01, 0.03))

#        if cc == 1 and key == 'Re':
#            for xx, label in zip([3e3, 1e5], ['Discontinuities', None]):
#                sp.axvline(xx, lw=2, color='red', ls='--', label=label)

        if cc == 2 and sp_ctr == 2:
            sp.legend(bbox_to_anchor=(1.05,1))
        if cc == 2 and key in ('Re', 'P1'):
            sp.legend(loc=0)

if args.savefig:
    for fig in figs:
        fig.savefig(os.path.expanduser(args.savefig) + '_%i.png' % fig.number)

if not args.noshow:
    plt.show()
