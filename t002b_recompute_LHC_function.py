import numpy as np

import LHCMeasurementTools.TimberManager as tm
import LHCMeasurementTools.LHC_Heatloads as HL

from calibration_config import calibration_config
from calibration import Calibration, CalibrationManager
from h5_storage import H5_storage
import heatload_recalc as hlr

cal_manager = CalibrationManager(calibration_config=calibration_config)
h5_storage = H5_storage(h5_dir='/eos/user/l/lhcecld/heatload_data_storage')
with_P_drop = True

filln = 6737

obraw = h5_storage.load_data_file(filln=filln)
calibration = cal_manager.get_calibration(obraw.timestamps[0])


def recalc_multiple_circuits(raw_data_object, calibration,
        circuit_selection, with_P_drop):

    if circuit_selection == 'full_lhc':
        circuits = calibration.circuits
    else:
        raise ValueError('Not implemented!')

    qbs_recalc = []
    issues = []
    for ii, circuit in enumerate(circuits):

        if len(circuits) > 100:
            if np.mod(ii, 20) == 0:
                print('Circuit %d/%d'%(ii, len(circuits)))
        else:
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
                with_P_drop=with_P_drop, N_iter_max=100, scale_correction=0.3,
                iter_toll=1e-3)

        qbs_recalc.append(Q_bs)

        if len(other['issues'])>0:
            print('Issues found for circuit %s:'%circuit)
            print('\n'.join(other['issues']))
            issues.append([circuit, other['issues']])

    avg_loads = []
    avg_varnames = []

    if circuit_selection == 'full_lhc':
        # Build temporary object to compute arc averages
        obhl = tm.AlignedTimberData(timestamps=obraw.timestamps,
                data=np.array(qbs_recalc).T, variables=calibration.circuits)

        # Compute arc averages
        for arc in '12 23 34 45 56 67 78 81'.split():
            arc_circuits =  HL.arc_cells_by_sector['S'+arc]
            arc_loads = np.array([obhl.dictionary[kk] for kk in arc_circuits])
            avg_load = np.nanmean(arc_loads, axis=0)

            avg_loads.append(avg_load)
            avg_varnames.append('S%s_QBS_AVG_ARC.POSST'%arc)

    # Build temporary object to compute arc averages
    obhl_store = tm.AlignedTimberData(timestamps=obraw.timestamps,
            data=np.array(qbs_recalc + avg_loads).T,
            variables=(calibration.circuits + avg_varnames))

    other = {}
    other['issues'] = issues

    return obhl_store, other

h5_storage.store_qbs(filln=filln, qbs_ob=obhl_store, use_dP=with_P_drop)
