import numpy as np
import LHCMeasurementTools.TimberManager as tm

from calibration_config import calibration_config
from calibration import Calibration, CalibrationManager
from h5_storage import H5_storage
import heatload_recalc as hlr

cal_manager = CalibrationManager(calibration_config=calibration_config)

filln = 6737

h5_storage = H5_storage(h5_dir='/eos/user/l/lhcecld/heatload_data_storage')
obraw = h5_storage.load_data_file(filln=filln)

calibration = cal_manager.get_calibration(obraw.timestamps[0])

qbs_recalc = []
issues = []
for ii, circuit in enumerate(calibration.circuits):
    print(ii, circuit)

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
            with_P_drop=True, N_iter_max=100, scale_correction=0.3,
            iter_toll=1e-3)

    qbs_recalc.append(Q_bs)

    if len(other['issues'])>0:
        print('\n'.join(other['issues']))
        issues.append([circuit, other['issues']])

qbs_recalc = np.array(qbs_recalc)

obhl = tm.AlignedTimberData(timestamps=obraw.timestamps,
        data=qbs_recalc.T, variables=calibration.circuits)

