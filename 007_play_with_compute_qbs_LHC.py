import numpy as np
import matplotlib.pyplot as plt

import h5_storage
import qbs_fill as qf
import compute_QBS_LHC as cql
import LHCMeasurementTools.mystyle as ms

from data_qbs import data_qbs

ms.mystyle_arial()

plt.close('all')
filln = 5219
atd = h5_storage.load_data_file(filln)
hlc = cql.HeatLoadComputer(atd, data_qbs, use_dP=True, report=True)

EH = hlc.data_dict['EH']
print sum([np.any(EH[:,i] <= 0) for i in xrange(EH.shape[1])])


qbs_ob = hlc.qbs_atd
qbs_ob_old = qf.compute_qbs_fill(filln)

arc_data = qf.compute_qbs_arc_avg(qbs_ob)
arc_data_old = qf.compute_qbs_arc_avg(qbs_ob_old)

fig = ms.figure('Compare for fill %i' % filln)
sp = plt.subplot(1,1,1)
sp.grid('on')

for ctr, (arc, data) in enumerate(arc_data.dictionary.iteritems()):
    color = ms.colorprog(ctr, arc_data.dictionary)
    sp.plot(arc_data.timestamps, data, label=arc, color=color)
    sp.plot(arc_data_old.timestamps, arc_data_old.dictionary[arc], ls ='--', color=color, label='Stored')

sp.legend(bbox_to_anchor=(1.1,1))

plt.show()
