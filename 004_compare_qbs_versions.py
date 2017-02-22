from __future__ import division
import sys
import matplotlib.pyplot as plt
import numpy as np
import h5_storage
import qbs_fill as qf

import LHCMeasurementTools.mystyle as ms

import compute_QBS_LHC as cql

plt.close('all')
ms.mystyle()

filln = 5219
use_dPs = (True, False)

qbs_obs = []
qbs_obs_loaded = []
cql_obs = []

atd_ob = h5_storage.load_data_file(filln)

for dp in use_dPs:

    qbs_obs.append(qf.test_compute_qbs(filln, use_dP=dp))
    cql_obs.append(cql.compute_qbs(atd_ob, use_dP=dp))
    qbs_obs_loaded.append(qf.compute_qbs_fill(filln, use_dP=dp))


arc_averages = [qf.compute_qbs_arc_avg(qo) for qo in qbs_obs]


fig = ms.figure('Compare dp')

sp = plt.subplot(2,2,1)
sp.grid(True)

for ctr in xrange(len(use_dPs)):
    arc_average = arc_averages[ctr]
    qbs_ob = qbs_obs[ctr]
    ls = ['-', '--', '-.'][ctr]

    tt = (qbs_ob.timestamps - qbs_ob.timestamps[0])/3600.

    for ctr, (key, hl) in enumerate(arc_average.dictionary.iteritems()):
        color = ms.colorprog(ctr, arc_average.dictionary)
        sp.plot(tt, hl, label=key, ls=ls)


plt.show()
