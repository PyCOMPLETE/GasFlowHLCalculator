from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
import cPickle as pickle
import time

import Helium_properties as hp
from Helium_properties import interp_P_T_hPT, interp_P_T_DPT, interp_P_T_mu

import LHCMeasurementTools.mystyle as ms

plt.close('all')
ms.mystyle()

time_0 = time.time()
def timer(n):
    print(n, time.time() - time_0)

filln = 5219
interp_P_T_gamma = interp2d(hp.P, hp.T, hp.gamma_PT)
t_arr = np.arange(4,50,0.5)

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
        sp.legend(bbox_to_anchor=(1.2,1))

with open('hlc_%i.pkl' % filln) as f:
    hlc = pickle.load(f)

combined_dict = hlc.data_dict.copy()
combined_dict.update(hlc.computed_values)
combined_dict['Pressure_ratio'] = combined_dict['P4'] / combined_dict['P3']


tt = hlc.atd_ob.timestamps
index_tt = np.argmin(np.abs(tt - tt[0] - 3600*2))

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

        if key == 'T2':
            data2 = data2[data2 < 30]

        if key in hlc.data_dict:
            affix = '(data)'
        else:
            affix = '(computed)'
        sp_ctr = ctr % 4 +1
        if cc == 1 and sp_ctr == 1:
            fig = ms.figure(plot_title)
        if cc == 1:
            sp = plt.subplot(2,2,sp_ctr)
            key_sp_dict[key] = sp
        else:
            sp = key_sp_dict[key]

        sp.set_title(key + ' ' + affix)
        sp.set_xlabel(key)
        sp.set_ylabel('Frequency')
        sp.grid(True)

        sp.hist(data2[data2 != 0], normed=True, label=label)

        if cc == 2 and sp_ctr == 2:
            sp.legend(bbox_to_anchor=(1.1,1))

#
#def A(x, gamma):
#    limit = (2/(gamma+1))**(gamma/(gamma-1))
#    a = 2*gamma/(gamma-1)
#    b = 2/gamma
#    c = (gamma-1) / gamma
#
#    A = 1.379 * np.sqrt(a*x**b * (1-x**c))
#    return np.where(x < limit, 1, A)
#
#fig = ms.figure('Samson A function')
#sp = plt.subplot(2,2,1)
#sp.set_title('A')
#sp.set_xlabel('Pressure ratio')
#sp.grid(True)
#
#gamma_arr = np.arange(1.5, 2.2, 0.1)
#xx_arr = np.arange(0.2, 0.61, 0.01)
#
#for ctr, gamma in enumerate(gamma_arr):
#    sp.plot(xx_arr, A(xx_arr, gamma), lw=2, label=gamma, color=ms.colorprog(ctr, gamma_arr))
#sp.plot(xx_arr, A(xx_arr, 5./3.), lw=2, label='5/3', color='black')
#sp.legend(bbox_to_anchor=(1.2,1))
#

plt.show()
