import numpy as np
import hepak as hep
from pycryo import pressureDrop
from pycryo import valve_m

from . import Helium_properties as hp
from .valve_LT import valve_LT
from .Pressure_drop import pd_factory


def compute_heat_load(P1, T1, T3, P4, CV1,CV2, EH, Qs, Kvmax, R, u0,
        cell_length, n_channels, channel_radius, channel_roughness,
        with_P_drop=True, N_iter_max=10, scale_correction=0.3,
        iter_toll=1e-3):

    # create all vectors
    n = len(P1)
    counter= np.zeros(n)
    P3_temp = np.zeros(n)
    m_L = np.zeros(n)
    H2=np.zeros(n)
    H3=np.zeros(n)
    rho3=np.zeros(n)
    g3=np.zeros(n)
    DP=np.zeros(n)
    mu_avg=np.zeros(n)
    rho_avg = np.zeros(n)
    DP_old = 0
    list_issues = []
    
    #linear valve or EQ valve ?
    if R == 0.0:
        Type = "LIN"
    else:
        Type = "EQ"
        
    
    #any zeros values in pressure/temperatures
    if len(P1[P1<=1.0]):
        list_issues.append("problem P1 null values")
    if len(P4[P4<=1.0]):
        list_issues.append("problem P4 null values")
    if len(T1[T1<=1.0]):
        list_issues.append("problem T1 null values")
    if len(T3[T3<=1.0]) or np.std(T3) == 0:
        list_issues.append("problem T3 null values")
    if (len(CV1[CV1<=5.0])):
        list_issues.append("problem CV94x closed")
    
    #replace zero values by default values if any
    P1[P1<=1.0] = 3.0
    P4[P4<1.0] = 1.2
    T1[T1<=1.0] = 5.0
    T3[T3<=1.0] = 20
    
    P3 = P1.copy()
    
    #Inlet enthalpy (line C)
    H1 = hep.calculate(9,1,P1*1e5,2,T1)
    
    #Pressure drop calculation
    if with_P_drop:
        for i in range(n):
            while(np.abs(P3_temp[i]-P3[i])/P3[i] > iter_toll and counter[i] < N_iter_max):
                if P3[i] < P4[i] or np.isnan(P3[i]):
                    list_issues.append("Pb in dP(P3<P4) at i="+str(i))
                    P3[i] = P1[i]
                    break
                elif P3_temp[i] >0:
                    P3[i]=P3_temp[i]

                rho3[i] = hep.calculate(3,1,P3[i]*1e5,2,T3[i])
                g3[i] = hep.calculate(16,1,P3[i]*1e5,2,T3[i])
                H3[i] = hep.calculate(9,1,P3[i]*1e5,2,T3[i]) 
                
                m_L[i] =valve_m(pin=P3[i],pout=P4[i], rho=rho3[i],g=g3[i],Kv0=Kvmax, R=R, Pos=CV1[i],u0=u0,phase="gas",Type=Type)
                
                #if second valve, add it
                if CV2 is not None and CV2[i] != 0.0:
                    m_L[i] =m_L[i] + valve_m(pin=P3[i],pout=P4[i], rho=rho3[i],g=g3[i],Kv0=Kvmax, R=R, Pos=CV2[i],u0=u0,phase="gas",Type=Type)

                if m_L[i] > 1e-5:
                    H2[i] = (m_L[i] * H1[i] + EH[i]+Qs) / m_L[i]
                    rho_avg[i] = hep.calculate(3,1,(P1[i]+P3[i])*1e5/2,9, (H2[i]+H3[i])/2)
                    mu_avg[i] = hep.calculate(25,1,(P1[i]+P3[i])*1e5/2,9, (H2[i]+H3[i])/2)
                    DP[i] = pressureDrop(m=m_L[i]/n_channels, mu=mu_avg[i] , rho=rho_avg[i],D=2*channel_radius,L=cell_length,rug=channel_roughness)
                else:
                    H2[i] = H3[i]
                    rho_avg[i] =0.0
                    mu_avg[i] = 0.0
                    DP[i] = 0.0

                DP_new = DP_old + (DP[i]-DP_old)*scale_correction   #smooth convergence
                P3_temp[i] = P1[i]-DP_new
                DP_old = DP_new
                counter[i] = counter[i]+1
                if counter[i] == N_iter_max and i != 0:
                    list_issues.append("max iteration (%i) reached at sample %i"%(N_iter_max,i))
                 
        P3 = P3_temp
    
    # Re-evaluate pressure, enthalpy, density before the valve for final massflow calculation
    P3[P3<=1.0] = 3.0
    H3 = hep.calculate(9,1,P3*1e5,2,T3)
    rho3 = hep.calculate(3,1,P3*1e5,2,T3)
    g3 = hep.calculate(16,1,P3*1e5,2, T3)
    m_L = valve_m(pin=P3,pout=P4, rho=rho3,g=g3,Kv0=Kvmax, R=R, Pos=CV1,u0=u0,phase="gas",Type=Type)

    #if second valve, add it
    if CV2 is not None and CV2.any() != 0.0:
        m_L = m_L + valve_m(pin=P3,pout=P4, rho=rho3,g=g3,Kv0=Kvmax, R=R, Pos=CV2,u0=u0,phase="gas",Type=Type)

    # Compute heat load
    Q_bs = m_L * (H3 - H1) - Qs - EH

    # # list all issues and remove invalid data
    # mask_invalid_data = (T3==0) | (CV1==0)
    # N_invalid = np.sum(mask_invalid_data)
    # if N_invalid > 0:
    #     Q_bs[mask_invalid_data] = np.nan
    #     list_issues.append('P3 or CV is 0 in %d points of %d.'%(N_invalid,len(mask_invalid_data)))

    # mask_P1T1P4is0 = (P4==0) | (T1==0) | (P1==0)
    # N_P1T1P4is0 = np.sum(mask_P1T1P4is0)
    # if N_P1T1P4is0 > 0:
    #     list_issues.append('P1, P4 or T1 is 0. in %d points of %d.'%(N_P1T1P4is0,len(mask_P1T1P4is0)))
    
    other = {
        'P3': P3,
        'mass_flow': m_L,
        'n_iter': counter,
        'issues': list_issues
        }

    return Q_bs,other

# Old function (possibly to be used for Run 1 data)
# def compute_heat_load(P1, T1, T3, P4, CV, EH, Qs_calib, Kv_calib, R_calib,
#         cell_length, n_channels, channel_radius, channel_roughness,
#         with_P_drop=True, N_iter_max=100, scale_correction=0.3,
#         iter_toll=1e-3):

#     # Evaluate enthalpy at circuit entrance
#     H1 = hp.interp_P_T_hPT(P1, T1)

#     # We initially neglect the pressure drop
#     P3_0 = P1.copy()

#     # Evaluate enthalpy at circuit exit 
#     H3_0 = hp.interp_P_T_hPT(P3_0, T3)

#     # Evaluate density at valve
#     rho = hp.interp_P_T_DPT(P3_0, T3)

#     # Evaluate the mass flow from the valve characteristics
#     m_L_0 = valve_LT(pin=P3_0, pout=P4, rho=rho, kv=Kv_calib,
#             u=CV, R=R_calib)

#     # Estimate pressure drop DP = P1 - P3
#     pressure_drop = pd_factory(D=2*channel_radius,
#                                rug=channel_roughness)

#     P3 = P3_0.copy()
#     P3_list = []
#     mask_negative = P4 > P1
#     n_iter = np.nan
#     if with_P_drop:
#         DP_prev = 0.
#         mask_iter = P3 > P4
#         n_iter = np.zeros_like(P4, dtype=np.int)
#         for i_iter in range(N_iter_max):

#             if np.sum(mask_iter) == 0:
#                 break

#             P3_prev_iter = P3[mask_iter].copy()
#             DP_prev_iter = P1[mask_iter] - P3_prev_iter

#             m_L_iter = valve_LT(pin=P3[mask_iter],
#                     pout=P4[mask_iter], rho=rho[mask_iter],
#                     kv=Kv_calib, u=CV[mask_iter], R=R_calib)

#             #H3_iter = H3_0[mask_iter]
#             H3_iter = hp.interp_P_T_hPT(P3[mask_iter], T3[mask_iter])
#             H2_iter = (m_L_iter * H1[mask_iter] + EH[mask_iter]) / m_L_iter
#             H_ave_iter = 0.5*(H2_iter + H3_iter)
#             rho_DP_iter = hp.interp_P_H_DPH(P3[mask_iter], H_ave_iter)
#             mu_iter = hp.inperp_P_H_mu(P3[mask_iter], H_ave_iter)

#             #T_ave = (T2 + T3)/2.
#             #rho_DP = hp.interp_P_T_DPT(P3, T_ave)
#             #mu = hp.interp_P_T_mu(P3, T_ave)

#             DP_new_iter = pressure_drop(m=m_L_iter/n_channels,
#                             L=cell_length, mu=mu_iter , rho=rho_DP_iter)

#             DP_iter = DP_prev_iter + scale_correction * (DP_new_iter - DP_prev_iter)

#             P3_iter = P1[mask_iter] - DP_iter

#             # Identify negative (P3-P4)
#             mask_negative_iter = P3_iter < P4[mask_iter]

#             # Stop iteration for negative values
#             mask_negative[mask_iter] = mask_negative[mask_iter] | mask_negative_iter
#             mask_iter[mask_iter] = (mask_iter[mask_iter]) & (~mask_negative_iter)

#             # Update only positive (P3-P4)  
#             P3[mask_iter] = P3_iter[~mask_negative_iter]

#             # Stop iteration for point where convergence is found
#             mask_iter[mask_iter] = (np.abs((P3_iter[~mask_negative_iter] \
#                     - P3_prev_iter[~mask_negative_iter])\
#                      / P3_prev_iter[~mask_negative_iter]) > iter_toll)

#             P3_list.append(P3.copy())

#             n_iter[mask_iter] = i_iter + 1

#         P3_list = np.array(P3_list)

#     # Re-evaluate mass flow
#     m_L = valve_LT(pin=P3, pout=P4, rho=rho, kv=Kv_calib,
#             u=CV, R=R_calib)

#     # Re-evaluate H3
#     H3 = hp.interp_P_T_hPT(P3, T3)

#     # Compute heat load
#     Q_bs = m_L * (H3 - H1) - Qs_calib - EH

#     list_issues = []

#     # Remove invalid data
#     mask_invalid_data = (T3==0) | (CV==0)
#     N_invalid = np.sum(mask_invalid_data)
#     if N_invalid > 0:
#         Q_bs[mask_invalid_data] = np.nan
#         list_issues.append('Invalid data in %d points of %d.'%(N_invalid,
#             len(mask_invalid_data)))

#     mask_P1T1P4is0 = (P4==0) | (T1==0) | (P1==0)
#     N_P1T1P4is0 = np.sum(mask_P1T1P4is0)
#     if N_P1T1P4is0 > 0:
#         list_issues.append('P1, P4 or T1 is 0. in %d points of %d.'%(N_P1T1P4is0,
#             len(mask_P1T1P4is0)))

#     N_negative = np.sum(mask_negative)
#     if N_negative > 0:
#         list_issues.append('Negative pressure drop in %d points of %d.'%(
#             N_negative, len(mask_negative)))

#     if with_P_drop:
#         N_convergence = np.sum(mask_iter)
#         if N_convergence > 0:
#             list_issues.append('No convergence in P drop for %d points of %d.'%(
#                 N_convergence, len(mask_iter)))


#     other = {
#         'H1': H1,
#         'H3': H3,
#         'P3': P3,
#         'mass_flow': m_L,
#         'P3_list': P3_list,
#         'n_iter': n_iter,
#         'issues': list_issues
#         }

#     return Q_bs, other

def compute_heat_loads_instrumented_cell(mass_flow, P1,
        T_in_magnets_circuits, T_out_magnets_circuits,
        magnet_lengths_circuits, n_channels_circuits,
        channel_radius, channel_roughness,
        N_iter_max=200, dp_toll = 0.001):

    m_L = mass_flow

    pressure_drop = pd_factory(D=2*channel_radius,
                           rug=channel_roughness)

    frac_flow_list = []
    dp_diff_list = []
    frac_flow = 0.5 + 0 * mass_flow

    for i_iter in range(N_iter_max):
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

        dp_diff = dp_circuits[0] - dp_circuits[1]
        dp_diff_list.append(dp_diff)
        frac_flow_list.append(frac_flow.copy())

        if np.sum(dp_diff>dp_toll) == 0:
            break

    # Final mass flow sharing
    mL_circuits = [m_L * frac_flow, m_L * (1. - frac_flow)]
    Qbs_magnets_circuits = [[], []]
    # Qbs for individual beam screens
    for i_circ in [0, 1]:
        mL_c = mL_circuits[i_circ]
        magnet_lengths_c = magnet_lengths_circuits[i_circ]
        T_in_magnets_c = T_in_magnets_circuits[i_circ]
        T_out_magnets_c = T_out_magnets_circuits[i_circ]

        for i_mag, lmag in enumerate(magnet_lengths_c):
            Hin = hp.interp_P_T_hPT(P1, T_in_magnets_c[i_mag])
            Hout = hp.interp_P_T_hPT(P1, T_out_magnets_c[i_mag])

            Qbs_mag = mL_c * (Hout - Hin)
            Qbs_magnets_circuits[i_circ].append(Qbs_mag)

    other = {
            'mass_flow_circuits': mL_circuits
        }

    return Qbs_magnets_circuits, other


def extract_info_from_instrum_config_dict(config_dict):

    n_channels_circuits = [config_dict['n_channels_circuit_%s'%cc]
                                for cc in ['A', 'B']]
    magnet_lengths_circuits = [config_dict['magnet_lengths']
                                    for _ in ['A', 'B']]
    out_sensor_names_circuits = [config_dict['circuit_%s_sensors'%cc][1:]
                                    for cc in ['A', 'B']]
    in_sensor_names_circuits = [config_dict['circuit_%s_sensors'%cc][:-1]
                                    for cc in ['A', 'B']]

    return (n_channels_circuits, magnet_lengths_circuits, in_sensor_names_circuits,
            out_sensor_names_circuits)


def build_instrumented_hl_dict(config_dict, circuit, Qbs_magnets_circuits):
    dict_output = {}

    magnet_beam_circuits = [config_dict['circuit_%s_beam'%cc]
                                    for cc in ['A', 'B']]

    magnet_names = config_dict['magnet_names']
    for i_circ in [0, 1]:
        magnets_beam_c = magnet_beam_circuits[i_circ]

        for i_mag, name_mag in enumerate(magnet_names):
            dict_output[circuit.split('.POSST')[0]
                   +'_%sB%s.POSST'%(name_mag, magnets_beam_c[i_mag])]=\
            Qbs_magnets_circuits[i_circ][i_mag]

    for name_mag in magnet_names:
        dict_output[circuit.split('.POSST')[0]
            +'_%s.POSST'%(name_mag)] = \
             dict_output[circuit.split('.POSST')[0]
                   +'_%sB%s.POSST'%(name_mag, 1)] +\
             dict_output[circuit.split('.POSST')[0]
                   +'_%sB%s.POSST'%(name_mag, 2)]

    # Hide last magnet
    name_last_magnet = magnet_names[-1]
    for kk in list(dict_output.keys()):
        for bb in [1,2]:
            if '_%sB%d'%(name_last_magnet, bb) in kk:
                dict_output[kk] *= 0.

    return dict_output



