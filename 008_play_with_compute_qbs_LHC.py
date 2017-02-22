import numpy as np
import matplotlib.pyplot as plt
#import sys
import pickle

import h5_storage
import qbs_fill as qf
import compute_QBS_LHC as cql
import LHCMeasurementTools.mystyle as ms

ms.mystyle_arial()

plt.close('all')
filln = 5219
atd = h5_storage.load_data_file(filln)
#hlc = cql.HeatLoadComputer(atd, use_dP=True, details=True)
with open('hlc_%i.pkl' % filln) as f:
    hlc = pickle.load(f)

EH = hlc.data_dict['EH']
P1 = hlc.data_dict['P1']
P3 = hlc.computed_values['P3']
P4 = hlc.data_dict['P4']

x = np.zeros_like(P1)
n_tt, n_cell = P1.shape

n_sub = []
n_both = []
n_else = []
for ctr, isnan in enumerate(hlc.nan_arr):
    if not isnan:
        x[:,ctr] = P4[:,ctr]/P3[:,ctr]
        if np.any(x[:,ctr] > .42):
            n_sub.append(ctr)
            if np.any(x[:,ctr] <= .42):
                n_both.append(ctr)
        else:
            n_else.append(ctr)

print(np.sum(x >= .42))
print(np.sum(x <= .42))
qbs_ob = hlc.qbs_atd

arc_data = qf.compute_qbs_arc_avg(qbs_ob)

fig = ms.figure('Compare for fill %i' % filln)

titles = ['Subcritical', 'Switching', 'Critical']

hl_tt = qbs_ob.nearest_older_sample(qbs_ob.timestamps[0]+2*3600)

for ctr, (list_, title) in enumerate(zip([n_sub, n_both, n_else],titles)):
    sp_ctr = ctr+1
    sp = plt.subplot(2,2,sp_ctr)
    sp.grid('on')
    sp.set_title(title + ' %i cells' % len(list_))
    sp.set_xlabel('Heat load [W/hc]')
    sp.set_ylabel('N cells')

    hl_arr = np.zeros(len(list_))
    for ctr_arr, ii in enumerate(list_):
        hl_arr[ctr_arr] = hl_tt[ii]

    sp.hist(hl_arr)

sp = plt.subplot(2,2,4)
sp.grid('on')
sp.set_xlabel('Heat load [W/hc]')
sp.set_ylabel('N cells')
hl_tt = filter(np.isfinite, hl_tt)
sp.set_title('LHC %i cells' % len(hl_tt))
sp.hist(hl_tt)

plt.show()
