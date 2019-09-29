import numpy as np
import LHCMeasurementTools.TimberManager as tm

from calibration_config import calibration_config
from calibration import Calibration, CalibrationManager
from h5_storage import H5_storage
import heatload_recalc as hlr

cal_manager = CalibrationManager(calibration_config=calibration_config)

with_P_drop = True

filln = 6737
circuit = 'QRLAB_23L2_QBS947.POSST' # Missing P4 (same result as logginh)
#circuit = 'QRLAB_15L2_QBS943.POSST' # Missing T1 (same result as logging)
#circuit = 'QRLAB_27L4_QBS943.POSST' # Missing P1 (result different from logging)
#circuit = 'QRLAB_31L2_QBS943.POSST'

h5_storage = H5_storage(h5_dir='/eos/user/l/lhcecld/heatload_data_storage')
obraw = h5_storage.load_data_file(filln=filln)

calibration = cal_manager.get_calibration(obraw.timestamps[0])


cell_calib = calibration.get_circuit(circuit)

T1 = obraw.dictionary[cell_calib['T1']]
T3 = obraw.dictionary[cell_calib['T3']]
P1 = obraw.dictionary[cell_calib['P1']]
P4 = obraw.dictionary[cell_calib['P4']]
CV= obraw.dictionary[cell_calib['CV']]
EH = obraw.dictionary[cell_calib['EH']]

# T2 = obraw.dictionary[cell_calib['T2']]

Q_bs, other = hlr.compute_heat_load(P1, T1, T3, P4, CV, EH,
        Qs_calib=cell_calib['Qs_calib'],
        Kv_calib=cell_calib['Kv_calib'],
        R_calib=cell_calib['R_calib'],
        cell_length=cell_calib['length'],
        n_channels=cell_calib['n_channels_tot'],
        channel_radius=cell_calib['channel_radius'],
        channel_roughness=cell_calib['roughness'],
        with_P_drop=with_P_drop, N_iter_max=100, scale_correction=0.3,
        iter_toll=1e-3)





