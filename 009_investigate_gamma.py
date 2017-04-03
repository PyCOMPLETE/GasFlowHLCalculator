from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
import cPickle as pickle
import time
import argparse

import Helium_properties as hp
from Helium_properties import interp_P_T_hPT, interp_P_T_DPT, interp_P_T_mu
import compute_QBS_LHC as qbl
import h5_storage

import LHCMeasurementTools.mystyle as ms
import LHCMeasurementTools.savefig as sf

parser = argparse.ArgumentParser()
parser.add_argument('--calc', help='Calculate instead of loading from pickle.', action='store_true')
parser.add_argument('--pdsave', help='Save in pdijksta plot dir.', action='store_true')
args = parser.parse_args()

recompute = args.calc
savefig   = args.pdsave
figs = []

try:
    from RcParams import init_pyplot
    init_pyplot(labelsize=22)
except ImportError:
    pass

plt.close('all')
ms.mystyle()

time_0 = time.time()
def timer(n):
    print(n, time.time() - time_0)

filln = 5219
interp_P_T_gamma = interp2d(hp.P, hp.T, hp.gamma_PT)
t_arr = np.arange(4,50,0.5)

with open('hlc_%i.pkl' % filln) as f:
    hlc = pickle.load(f)

gamma = np.zeros_like(hlc.computed_values['m_L'])
P3 = hlc.computed_values['P3']
T3 = hlc.data_dict['T3']
for j in xrange(hlc.Nvalue):
    for i in xrange(hlc.Ncell):
        gamma[j,i] = interp_P_T_gamma(P3[j,i], T3[j,i])

fig = ms.figure('Interpolated Data')


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
    for pressure in pressures:
        values = np.empty_like(t_arr)
        for i, temp in enumerate(t_arr):
            values[i] = interp(pressure, temp)
        sp.plot(t_arr, values, lw=2., label='P=%.1f bar' % pressure)

    if interp is interp_P_T_gamma:
        sp.axhline(5./3., lw=2, ls='--', label='5/3')

    if sp_ctr == 4:
        sp.legend(loc=1)

if recompute:
    atd = h5_storage.load_data_file(filln)
    hlc = qbl.HeatLoadComputer(atd, compute_Re=True)
else:
    with open('hlc_%i.pkl' % filln) as f:
        hlc = pickle.load(f)

combined_dict = hlc.data_dict.copy()
combined_dict.update(hlc.computed_values)
combined_dict['Pressure_ratio'] = combined_dict['P4'] / combined_dict['P3']
combined_dict['gamma'] = gamma

tt = hlc.atd_ob.timestamps
index_tt = np.argmin(np.abs(tt - tt[0] - 3600*2))

fig_ctr = 1
key_sp_dict = {}
for cc in (1,2):
    for ctr, (key, data) in enumerate(sorted(combined_dict.items())):

        if cc == 1:
            data2 = np.nan_to_num(data.flatten())
            plot_title = ('Occurrence of raw data values for fill %i' % filln)
            label = "All data"
        else:
            data2 = np.nan_to_num(data[index_tt,:])
            plot_title = ('Occurrence of raw data values for fill %i after 2 hrs' % filln)
            label = "After 2 hours"

        data2 = data2[data2 != 0]
        hist, bin_edges = np.histogram(data2, bins=10, normed=True)
        new_bin_edges = []
        for h_ctr, val in enumerate(hist):
            edge1, edge2 = bin_edges[h_ctr], bin_edges[h_ctr+1]
            binwidth = edge2 - edge1
            if val*binwidth > 0.01:
                new_bin_edges.extend([edge1, edge2])
        data2 = data2[np.logical_and(data2 > new_bin_edges[0], data2 < new_bin_edges[-1])]

        if key in hlc.data_dict:
            affix = '(data)'
        else:
            affix = '(computed)'
        sp_ctr = ctr % 4 +1
        if cc == 1 and sp_ctr == 1:
            fig = ms.figure(plot_title+'_%i' % fig_ctr, figs)
            fig_ctr += 1
        if cc == 1:
            sp = plt.subplot(2,2,sp_ctr)
            key_sp_dict[key] = sp
        else:
            sp = key_sp_dict[key]

        sp.set_title(key+' '+affix)
        sp.set_xlabel(key)
        sp.set_ylabel('Rel. Occurence')
        sp.set_yticklabels([])
        sp.grid(True)

        sp.hist(data2[data2 != 0], normed=True, label=label, alpha=0.5)

        if cc == 2 and sp_ctr == 2:
            sp.legend(bbox_to_anchor=(1.1,1))

        if cc == 1 and key == 'Re':
            for xx in 3e3, 1e5:
                sp.axvline(xx, lw=2, color='red', ls='--')

if savefig:
    sf.pdijksta(figs, dpi=150)

plt.show()
