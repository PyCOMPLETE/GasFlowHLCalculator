import numpy as np

from h5_storage import H5_storage
import Helium_properties as hp
from valve_LT import valve_LT
from Pressure_drop import pd_factory

cell_description = {
    'P1': 'QRLAA_29L2_PT961.POSST',
    'P4': 'QRLAA_29L2_PT991.POSST',
    'T1': 'QRLAA_29L2_TT961.POSST',
    'T3': 'QRLAB_31L2_TT943.POSST',

    'T2': 'LBARB_31L2_TT843.POSST',

    'CV': 'QRLAB_31L2_CV943.POSST',
    'EH': 'LBARB_31L2_EH843.POSST',

    'channel_radius': 3.7e-3/2.,
    'roughness': 1e-5,

    'n_channels_tot': 4,

    'length': 53.,

    'circuit_A_sensors': ['QBQI_31L2_TT825.POSST', 'QBBI_A31L2_TT824.POSST',
                          'QBBI_B31L2_TT826.POSST', 'QQBI_31L2_TT824.POSST',
                          'QBQI_32L2_TT825.POSST'],
    'circuit_B_sensors': ['QBQI_31L2_TT825.POSST', 'QBBI_A31L2_TT826.POSST',
                          'QBBI_B31L2_TT824.POSST', 'QQBI_31L2_TT826.POSST',
                          'QBQI_32L2_TT825.POSST'],

    'n_channels_circuit_A': 2,
    'n_channels_circuit_B': 2,


    'magnet_sequence': ['Q1', 'D2', 'D3', 'D4'],
    'magnet_lengths': [5.9, 15.7, 15.7, 15.7],

    'b1_circuit': ['A', 'B', 'A', 'B'],
    'b2_circuit': ['B', 'A', 'B', 'A'],
}

cell_calibration = {
 'R': 38.,
 'Qs': 8.,
 'Kv': 0.39
 }

h5_storage = H5_storage(h5_dir='/eos/user/l/lhcecld/heatload_data_storage')

obraw = h5_storage.load_special_data_file(filln=6737)

T1 = obraw.dictionary[cell_description['T1']]
T3 = obraw.dictionary[cell_description['T3']]
P1 = obraw.dictionary[cell_description['P1']]
P4 = obraw.dictionary[cell_description['P4']]
CV= obraw.dictionary[cell_description['CV']]
EH = obraw.dictionary[cell_description['EH']]

# T2 = obraw.dictionary[cell_description['T2']]

# Evaluate enthalpy at circuit entrance
H1 = hp.interp_P_T_hPT(P1, T1)

# We initially neglect the pressure drop
P3_0 = P1.copy()

# Evaluate enthalpy at circuit exit 
H3_0 = hp.interp_P_T_hPT(P3_0, T3)

# Evaluate density at valve
rho = hp.interp_P_T_DPT(P3_0, T3)

# Evaluate the mass flow from the valve characteristics
m_L_0 = valve_LT(pin=P3_0, pout=P4, rho=rho, kv=cell_calibration['Kv'],
        u=CV, R=cell_calibration['R'])

# Estimate pressure drop DP = P1 - P3
pressure_drop = pd_factory(D=2*cell_description['channel_radius'],
                           rug=cell_description['roughness'])

P3 = P3_0.copy()
P3_list = []

DP_prev = 0.
N_iter_max = 100
scale_correction = 0.3

mask_iter = P3 > P4

for i_iter in range(N_iter_max):

    P3_prev_iter = P3[mask_iter].copy()
    DP_prev_iter = P1[mask_iter] - P3_prev_iter

    m_L_iter = valve_LT(pin=P3[mask_iter], pout=P4[mask_iter], rho=rho[mask_iter],
            kv=cell_calibration['Kv'], u=CV[mask_iter], R=cell_calibration['R'])

    #H3_iter = H3_0[mask_iter]
    H3_iter = hp.interp_P_T_hPT(P3[mask_iter], T3[mask_iter])
    H2_iter = (m_L_iter * H1[mask_iter] + EH[mask_iter]) / m_L_iter
    H_ave_iter = 0.5*(H2_iter + H3_iter)
    rho_DP_iter = hp.interp_P_H_DPH(P3[mask_iter], H_ave_iter)
    mu_iter = hp.inperp_P_H_mu(P3[mask_iter], H_ave_iter)

    #T_ave = (T2 + T3)/2.
    #rho_DP = hp.interp_P_T_DPT(P3, T_ave)
    #mu = hp.interp_P_T_mu(P3, T_ave)

    DP_new_iter = pressure_drop(m=m_L_iter/cell_description['n_channels_tot'],
                    L=cell_description['length'], mu=mu_iter , rho=rho_DP_iter)

    DP_iter = DP_prev_iter + scale_correction * (DP_new_iter - DP_prev_iter)

    P3_iter = P1[mask_iter] - DP_iter

    # Identify negative (P3-P4)
    mask_negative_iter = P3_iter < P4[mask_iter]

    # Stop iteration for negative values
    mask_iter[mask_iter][mask_negative_iter] = False

    # Update only positive (P3-P4)  
    P3[mask_iter] = P3_iter[~mask_negative_iter]

    # Stop iteration for point where convergence is found
    mask_iter[mask_iter] = (np.abs((P3_iter[~mask_negative_iter] \
            - P3_prev_iter[~mask_negative_iter])\
             / P3_prev_iter[~mask_negative_iter]) > 0.001)

    if np.sum(mask_iter) == 0:
        break

    P3_list.append(P3.copy())


P3_list = np.array(P3_list)

# Re-evaluate mass flow
m_L = valve_LT(pin=P3, pout=P4, rho=rho, kv=cell_calibration['Kv'],
        u=CV, R=cell_calibration['R'])

# Re-evaluate H3
H3 = hp.interp_P_T_hPT(P3, T3)

# Compute heat load
Q_bs = m_L * (H3 - H1) - cell_calibration['Qs'] - EH

#####################################################################
# Compute share between the two beam screens

n_channels_circuits = [cell_description['n_channels_circuit_%s'%cc]
                            for cc in ['A', 'B']]
magnet_lengths_circuits = [cell_description['magnet_lengths']
                                for _ in ['A', 'B']]

out_sensor_names_circuits = [cell_description['circuit_%s_sensors'%cc][1:]
                                for cc in ['A', 'B']]

in_sensor_names_circuits = [cell_description['circuit_%s_sensors'%cc][:-1]
                                for cc in ['A', 'B']]

T_out_magnets_circuits = [[obraw.dictionary[vv] for vv in
        out_sensor_names_circuits[ii]] for ii in [0, 1]]

T_in_magnets_circuits = [[obraw.dictionary[vv] for vv in
        in_sensor_names_circuits[ii]] for ii in [0, 1]]


frac_flow_list = []
dp_diff_list = []
frac_flow = 0.5 + 0*Q_bs

for i_iter in xrange(N_iter_max):
    mL_circuits = [m_L * frac_flow, m_L * (1. - frac_flow)]
    dp_circuits = []
    for i_circ, mL_circuit in enumerate(mL_circuits):

        T_out_magnets = T_in_magnets_circuits[i_circ]
        magnet_lengths = magnet_lengths_circuits[i_circ]
        n_channels_circuit = n_channels_circuits[i_circ]

        rho_out_magnets = [
                hp.interp_P_T_DPT(P1, ttout) for ttout in T_out_magnets]
        mu_out_magnets = [
                hp.interp_P_T_mu(P1, ttout) for ttout in T_out_magnets]

        dp_magnets = [pressure_drop(m=mL_circuit/n_channels_circuit, # This is not there in Benjamin's implementation!!!!
                            L=ll, mu=mumu , rho=rhorho)
                                for ll, mumu, rhorho in zip(magnet_lengths,
                                        mu_out_magnets, rho_out_magnets)]

        dp_circuit = np.sum(np.array(dp_magnets), axis=0)

        dp_circuits.append(dp_circuit.copy())

    frac_flow *= (1 + 0.05*(dp_circuits[1] - dp_circuits[0])
            /(dp_circuits[1] + dp_circuits[0]))

    dp_diff_list.append(dp_circuits[0] - dp_circuits[1])
    frac_flow_list.append(frac_flow.copy())


# Final mass flow sharing
mL_circuits = [m_L * frac_flow, m_L * (1. - frac_flow)]

# Qbs for individual beam screens
