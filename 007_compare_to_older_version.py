from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

import qbs_fill as qf
import LHCMeasurementTools.mystyle as ms

ms.mystyle_arial()

plt.close('all')
filln = 5219

version = 5
old_version = 4

qbs_ob = qf.compute_qbs_fill(filln, version=version)
qbs_ob_old = qf.compute_qbs_fill(filln, version=old_version)

arc_data = qf.compute_qbs_arc_avg(qbs_ob)
arc_data_old = qf.compute_qbs_arc_avg(qbs_ob_old)

fig = ms.figure('Compare for fill %i' % filln)
sp = plt.subplot(2,2,1)
sp.grid('on')
sp.set_ylabel('Heat load [W/m]')
sp.set_xlabel('Time [h]')

tt = (arc_data.timestamps - arc_data.timestamps[0])/3600.

for ctr, (arc, data) in enumerate(arc_data.dictionary.iteritems()):
    color = ms.colorprog(ctr, arc_data.dictionary)
    sp.plot(tt, data, label=arc+' v%i'%version, color=color)
    sp.plot(tt, arc_data_old.dictionary[arc], ls ='--', color=color, label='v%i' % old_version)

sp.legend(bbox_to_anchor=(1.2,1))

plt.show()
