import numpy as np
import LHCMeasurementTools.TimberManager as tm

from calibration_config import calibration_config
from calibration import Calibration, CalibrationManager
from h5_storage import H5_storage
import heatload_recalc as hlr

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

if compute_instrumented:
    obraw = h5_storage.load_special_data_file(filln=filln)
else:
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

if compute_instrumented:
    from instrumented_cells_config import instrumented_cells_config

    instrum_cell_config = instrumented_cells_config[circuit]

    (n_channels_circuits, magnet_lengths_circuits, in_sensor_names_circuits,
            out_sensor_names_circuits) = \
                hlr.extract_info_from_instrum_config_dict(
                    config_dict=instrum_cell_config)

    T_out_magnets_circuits = [[obraw.dictionary[vv] for vv in
            out_sensor_names_circuits[ii]] for ii in [0, 1]]
    T_in_magnets_circuits = [[obraw.dictionary[vv] for vv in
            in_sensor_names_circuits[ii]] for ii in [0, 1]]

    Qbs_magnets_circuits, other_instr = hlr.compute_heat_loads_instrumented_cell(
        mass_flow = other['mass_flow'], P1=P1,
        T_in_magnets_circuits=T_in_magnets_circuits,
        T_out_magnets_circuits=T_out_magnets_circuits,
        magnet_lengths_circuits=magnet_lengths_circuits,
        n_channels_circuits=n_channels_circuits,
        channel_radius=cell_calib['channel_radius'],
        channel_roughness=cell_calib['roughness'],
        dp_toll = 0.001, N_iter_max=200)

    magnet_beam_circuits = [instrum_cell_config['circuit_%s_beam'%cc]
                                    for cc in ['A', 'B']]

    dict_output = hlr.build_instrumented_hl_dict(
            config_dict=instrum_cell_config, circuit=circuit,
            Qbs_magnets_circuits=Qbs_magnets_circuits)

    # Some plots
    import matplotlib.pyplot as plt
    plt.close('all')

    for i_mag, name_mag in enumerate(instrum_cell_config['magnet_names']):
        fig = plt.figure(i_mag+1)
        ax = fig.add_subplot(111)

        nn = circuit.replace('.POSST', '_%s.POSST'%name_mag)
        nnb1 = circuit.replace('.POSST', '_%sB1.POSST'%name_mag)
        nnb2 = circuit.replace('.POSST', '_%sB2.POSST'%name_mag)

        ax.plot(dict_output[nn], color='k')
        ax.plot(dict_output[nnb1], color='b')
        ax.plot(dict_output[nnb2], color='r')

        fig.suptitle(name_mag)
    plt.show()

