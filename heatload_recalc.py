import numpy as np

import Helium_properties as hp
from valve_LT import valve_LT
from Pressure_drop import pd_factory


# Functions
def compute_heat_load(P1, T1, T3, P4, CV, EH, Qs_calib, Kv_calib, R_calib,
        cell_length, n_channels, channel_radius, channel_roughness,
        with_P_drop=True, N_iter_max=100, scale_correction=0.3,
        iter_toll=1e-3):

    # Evaluate enthalpy at circuit entrance
    H1 = hp.interp_P_T_hPT(P1, T1)

    # We initially neglect the pressure drop
    P3_0 = P1.copy()

    # Evaluate enthalpy at circuit exit 
    H3_0 = hp.interp_P_T_hPT(P3_0, T3)

    # Evaluate density at valve
    rho = hp.interp_P_T_DPT(P3_0, T3)

    # Evaluate the mass flow from the valve characteristics
    m_L_0 = valve_LT(pin=P3_0, pout=P4, rho=rho, kv=Kv_calib,
            u=CV, R=R_calib)

    # Estimate pressure drop DP = P1 - P3
    pressure_drop = pd_factory(D=2*channel_radius,
                               rug=channel_roughness)

    P3 = P3_0.copy()
    P3_list = []
    if with_P_drop:
        DP_prev = 0.
        mask_iter = P3 > P4
        mask_negative = P4 > P1
        n_iter = np.zeros_like(P4, dtype=np.int)
        for i_iter in range(N_iter_max):

            if np.sum(mask_iter) == 0:
                break

            P3_prev_iter = P3[mask_iter].copy()
            DP_prev_iter = P1[mask_iter] - P3_prev_iter

            m_L_iter = valve_LT(pin=P3[mask_iter],
                    pout=P4[mask_iter], rho=rho[mask_iter],
                    kv=Kv_calib, u=CV[mask_iter], R=R_calib)

            #H3_iter = H3_0[mask_iter]
            H3_iter = hp.interp_P_T_hPT(P3[mask_iter], T3[mask_iter])
            H2_iter = (m_L_iter * H1[mask_iter] + EH[mask_iter]) / m_L_iter
            H_ave_iter = 0.5*(H2_iter + H3_iter)
            rho_DP_iter = hp.interp_P_H_DPH(P3[mask_iter], H_ave_iter)
            mu_iter = hp.inperp_P_H_mu(P3[mask_iter], H_ave_iter)

            #T_ave = (T2 + T3)/2.
            #rho_DP = hp.interp_P_T_DPT(P3, T_ave)
            #mu = hp.interp_P_T_mu(P3, T_ave)

            DP_new_iter = pressure_drop(m=m_L_iter/n_channels,
                            L=cell_length, mu=mu_iter , rho=rho_DP_iter)

            DP_iter = DP_prev_iter + scale_correction * (DP_new_iter - DP_prev_iter)

            P3_iter = P1[mask_iter] - DP_iter

            # Identify negative (P3-P4)
            mask_negative_iter = P3_iter < P4[mask_iter]

            # Stop iteration for negative values
            mask_negative[mask_iter] = mask_negative[mask_iter] | mask_negative_iter
            mask_iter[mask_iter] = (mask_iter[mask_iter]) & (~mask_negative_iter)

            # Update only positive (P3-P4)  
            P3[mask_iter] = P3_iter[~mask_negative_iter]

            # Stop iteration for point where convergence is found
            mask_iter[mask_iter] = (np.abs((P3_iter[~mask_negative_iter] \
                    - P3_prev_iter[~mask_negative_iter])\
                     / P3_prev_iter[~mask_negative_iter]) > iter_toll)

            P3_list.append(P3.copy())

            n_iter[mask_iter] = i_iter + 1

        P3_list = np.array(P3_list)

    # Re-evaluate mass flow
    m_L = valve_LT(pin=P3, pout=P4, rho=rho, kv=Kv_calib,
            u=CV, R=R_calib)

    # Re-evaluate H3
    H3 = hp.interp_P_T_hPT(P3, T3)

    # Compute heat load
    Q_bs = m_L * (H3 - H1) - Qs_calib - EH

    list_issues = []

    # Remove invalid data
    mask_invalid_data = (P1 == 0) | (P4==0) | (T1==0) | (T3==0) | (CV==0)
    N_invalid = np.sum(mask_invalid_data)
    if N_invalid > 0:
        Q_bs[mask_invalid_data] = np.nan
        list_issues.append('Invalid data in %d points of %d.'%(N_invalid,
            len(mask_invalid_data)))

    N_negative = np.sum(mask_negative)
    if N_negative > 0:
        Q_bs[mask_negative] = np.nan
        list_issues.append('Negative pressure drop in %d points of %d.'%(
            N_negative, len(mask_negative)))

    N_convergence = np.sum(mask_iter)
    if N_convergence > 0:
        Q_bs[mask_iter] = np.nan
        list_issues.append('No convergence in P drop for %d points of %d.'%(
            N_convergence, len(mask_iter)))


    other = {
        'H1': H1,
        'H3': H3,
        'P3': P3,
        'mass_flow': m_L,
        'P3_list': P3_list,
        'n_iter': n_iter,
        'issues': list_issues
        }

    return Q_bs, other

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



