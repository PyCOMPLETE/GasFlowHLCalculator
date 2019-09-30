import numpy as np

import LHCMeasurementTools.TimberManager as tm
import LHCMeasurementTools.LHC_Heatloads as HL

import heatload_recalc as hlr

from instrumented_cells_config import instrumented_cells_config

def recalc_multiple_circuits(raw_data_object, calibration,
        circuit_selection, with_P_drop):

    if circuit_selection == 'full_lhc':
        circuits = calibration.circuits
    elif circuit_selection == 'all_instrumented':
        circuits = sorted(instrumented_cells_config.keys())
    else:
        raise ValueError('Not implemented!')

    obraw = raw_data_object

    qbs_recalc = []
    issues = []
    instrum_cell_recalc_dict = {}
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

        if circuit_selection == 'all_instrumented':
            instrum_cell_config = instrumented_cells_config[circuit]

            (n_channels_circuits, magnet_lengths_circuits,
                    in_sensor_names_circuits, out_sensor_names_circuits
                    ) = hlr.extract_info_from_instrum_config_dict(
                    config_dict=instrum_cell_config)

            t_out_magnets_circuits = [[obraw.dictionary[vv] for vv in
                out_sensor_names_circuits[ii]] for ii in [0, 1]]
            t_in_magnets_circuits = [[obraw.dictionary[vv] for vv in
                in_sensor_names_circuits[ii]] for ii in [0, 1]]

            qbs_magnets_circuits, other_instr = \
                hlr.compute_heat_loads_instrumented_cell(
                    mass_flow = other['mass_flow'], p1=p1,
                    t_in_magnets_circuits=t_in_magnets_circuits,
                    t_out_magnets_circuits=t_out_magnets_circuits,
                    magnet_lengths_circuits=magnet_lengths_circuits,
                    n_channels_circuits=n_channels_circuits,
                    channel_radius=cell_calib['channel_radius'],
                    channel_roughness=cell_calib['roughness'],
                    dp_toll = 0.001, n_iter_max=200)

            magnet_beam_circuits = [
                    instrum_cell_config['circuit_%s_beam'%cc]
                                    for cc in ['a', 'b']]

            dict_output = hlr.build_instrumented_hl_dict(
                config_dict=instrum_cell_config, circuit=circuit,
                 qbs_magnets_circuits=qbs_magnets_circuits)

            instrum_cell_recalc_dict.update(dict_output)

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

    instrum_varnames = sorted(instrum_cell_recalc_dict.keys())
    instrum_qbs_recalc = [instrum_cell_recalc_dict[kk]
                            for kk in instrum_varnames]

    obhl_store = tm.AlignedTimberData(timestamps=obraw.timestamps,
            data=np.array(qbs_recalc + avg_loads + instrum_qbs_recalc).T,
            variables=(calibration.circuits + avg_varnames + instrum_varnames))

    other = {}
    other['issues'] = issues

    return obhl_store, other

