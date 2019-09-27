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

    'circuit_1_sensors': ['QBQI_31L2_TT825.POSST', 'QBBI_A31L2_TT824.POSST',
                          'QBBI_B31L2_TT826.POSST', 'QQBI_31L2_TT824.POSST',
                          'QBQI_32L2_TT825.POSST'],
    'circuit_2_sensors': ['QBQI_31L2_TT825.POSST', 'QBBI_A31L2_TT826.POSST',
                          'QBBI_B31L2_TT824.POSST', 'QQBI_31L2_TT826.POSST',
                          'QBQI_32L2_TT825.POSST'],

    'magnet_sequence': ['Q1', 'D2', 'D3', 'D4'],

    'b1_circuit': [1,2,1,2],
    'b2_circuit': [2,1,2,1],
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

T2 = obraw.dictionary[cell_description['T2']]

# Evaluate density at circuit entrance
rho = hp.interp_P_T_DPT(P1, T1)

# Evaluate enthalpy at circuit entrance
H1 = hp.interp_P_T_hPT(P1, T1)

# We initially neglect the pressure drop
P3_0 = P1

# Evaluate enthalpy at circuit exit 
H3_0 = hp.interp_P_T_hPT(P3_0, T3)

# Evaluate the mass flow from the valve characteristics
m_L_0 = valve_LT(pin=P3_0, pout=P4, rho=rho, kv=cell_calibration['Kv'],
        u=CV, R=cell_calibration['R'])

# Estimate pressure drop DP = P1 - P3
pressure_drop = pd_factory(D=2*cell_description['channel_radius'],
                           rug=cell_description['roughness'])

P3 = P3_0.copy()
m_L = m_L_0.copy()
H3 = H3_0.copy() # In fact it is not updated in the calculation

P3_list = []

for _ in range(30):

    m_L =valve_LT(pin=P3, pout=P4, rho=rho, kv=cell_calibration['Kv'],
        u=CV, R=cell_calibration['R'])

    #H3 = hp.interp_P_T_hPT(P3, T3)
    H2 = (m_L * H1 + EH) / m_L
    H_ave = 0.5*(H2 + H3)
    rho_DP = hp.interp_P_H_DPH(P3, H_ave)
    mu = hp.inperp_P_H_mu(P3, H_ave)

    #T_ave = (T2 + T3)/2.
    #rho_DP = hp.interp_P_T_DPT(P3, T_ave)
    #mu = hp.interp_P_T_mu(P3, T_ave)

    DP = pressure_drop(m=m_L/cell_description['n_channels_tot'],
            L=cell_description['length'], mu=mu , rho=rho_DP)
    P3 = P1 - DP

    P3_list.append(P3.copy())

P3_list = np.array(P3_list)

