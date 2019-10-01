import numpy as np
import LHCMeasurementTools.TimberManager as tm

from calibration_config import calibration_config
from calibration import Calibration, CalibrationManager
from h5_storage import H5_storage
import heatload_recalc as hlr

import GasFlowHLCalculator.qbs_fill as qf

cal_manager = CalibrationManager(calibration_config=calibration_config)

with_P_drop = True
compute_instrumented = True

filln = 6737
filln = 6966
filln = 6967

circuit = 'QRLAB_23L2_QBS947.POSST' # Missing P4 (same result as logginh)
#circuit = 'QRLAB_15L2_QBS943.POSST' # Missing T1 (same result as logging)
#circuit = 'QRLAB_27L4_QBS943.POSST' # Missing P1 (result different from logging)
circuit = 'QRLAB_31L2_QBS943.POSST' # Instrumented cell
circuit = 'QRLAA_13L5_QBS943.POSST' # Instrumented cell
circuit = 'QRLAA_13R4_QBS947.POSST' # Instrumented cell
circuit = 'QRLAA_33L5_QBS947.POSST' # Instrumented cell

h5_storage = H5_storage(h5_dir='/eos/user/l/lhcecld/heatload_data_storage')
# Some plots
import matplotlib.pyplot as plt
plt.close('all')

fill_dict = qf.get_fill_dict(filln, h5_storage=h5_storage, use_dP=with_P_drop)

t_h = fill_dict[circuit].t_stamps

figtot = plt.figure(100)
axtot = figtot.add_subplot(111)

axtot.plot(t_h, fill_dict[circuit].values, color = 'black', label="Total")

axlist = [axtot]

if compute_instrumented:
    for i_mag, name_mag in enumerate('Q1 D2 D3 D4'.split()):
        fig = plt.figure(i_mag+1)
        ax = fig.add_subplot(111, sharex=axtot, sharey=axtot)

        nn = circuit.replace('.POSST', '_%s.POSST'%name_mag)
        nnb1 = circuit.replace('.POSST', '_%sB1.POSST'%name_mag)
        nnb2 = circuit.replace('.POSST', '_%sB2.POSST'%name_mag)

        ax.plot(t_h, fill_dict[nn].values, color='k', label='Total')
        ax.plot(t_h, fill_dict[nnb1].values, color='b', label='B1')
        ax.plot(t_h, fill_dict[nnb2].values, color='r', label='B2')

        fig.suptitle(circuit + ' - ' + name_mag)

        axtot.plot(t_h, fill_dict[nn].values, label=name_mag)
        axlist.append(ax)

for aa in axlist:
    aa.legend(loc='best')
    aa.set_xlabel('t [h]')
    aa.set_ylabel('Heat load [W]')
    aa.grid(True, linestyle='--', alpha=.5)

figtot.suptitle(circuit)
plt.show()

