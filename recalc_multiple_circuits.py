import numpy as np

import LHCMeasurementTools.TimberManager as tm
import LHCMeasurementTools.LHC_Heatloads as HL
import hepak as hep
import pycryo
from . import heatload_recalc as hlr

from .instrumented_cells_config import instrumented_cells_config

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
                print(('Circuit %d/%d'%(ii, len(circuits))))
        else:
            print((ii, circuit))

        cell_calib = calibration.get_circuit(circuit)

        T1 = obraw.dictionary[cell_calib['T1']]
        T3 = obraw.dictionary[cell_calib['T3']]
        P1 = obraw.dictionary[cell_calib['P1']]
        P4 = obraw.dictionary[cell_calib['P4']]
        CV1 = obraw.dictionary[cell_calib['CV1']]
        CV2 = obraw.dictionary.get(cell_calib['CV2'], None)
        EH = obraw.dictionary[cell_calib['EH']]

        # T2 = obraw.dictionary[cell_calib['T2']]

        Q_bs, other = hlr.compute_heat_load(P1=P1, T1=T1, T3=T3, P4=P4, CV1=CV1, CV2=CV2, EH=EH,
                Qs=cell_calib['Qs_calib'],
                Kvmax=cell_calib['Kv_calib'],
                R=cell_calib['R_calib'],
                u0=cell_calib['u0_calib'],
                cell_length=cell_calib['length'],
                n_channels=cell_calib['n_channels_tot'],
                channel_radius=cell_calib['channel_radius'],
                channel_roughness=cell_calib['roughness'],
                with_P_drop=with_P_drop, N_iter_max=100, scale_correction=0.3,
                iter_toll=1e-3)

        qbs_recalc.append(Q_bs)

        if len(other['issues'])>0:
            print(('Issues found for circuit %s:'%circuit))
            print(('\n'.join(other['issues'])))
            issues.append([circuit, other['issues']])

        if circuit_selection == 'all_instrumented':
            instrum_cell_config = instrumented_cells_config[circuit]

            # (n_channels_circuits, magnet_lengths_circuits,
            #         in_sensor_names_circuits, out_sensor_names_circuits
            #         ) = hlr.extract_info_from_instrum_config_dict(
            #         config_dict=instrum_cell_config)

            C1TT = [obraw.dictionary[vv] for vv in instrum_cell_config['circuit_A_sensors']]
            C2TT = [obraw.dictionary[vv] for vv in instrum_cell_config['circuit_B_sensors']]
            C1FE = obraw.dictionary[instrum_cell_config["C1-FE"]]/1000
            C2FE = obraw.dictionary[instrum_cell_config["C2-FE"]]/1000
            faults1 = []
            faults2 = []
            for ii in range(len(C1TT)):
                if np.all(C1TT[ii] == 0):
                    print(f'Faulty temperature sensor: ', instrum_cell_config['circuit_A_sensors'][ii]) 
                    faults1.append(ii)
                if np.all(C2TT[ii] == 0):
                    print(f'Faulty temperature sensor: ', instrum_cell_config['circuit_B_sensors'][ii]) 
                    faults2.append(ii)
            
            hc_type = int("QBS943" in circuit) # 0 = QDDD (X47) / 1 = DDDQ (X43)
            # don't use coriolis yet
            coriolis=True
            for coriolis in [True, False]:
                if coriolis:
                    #pass
                    FT1  = C1FE
                    FT2  = C2FE
                    mass_flow = FT1 + FT2
                    H1 = hep.calculate(9,1,P1*1e5,2,T1)
                    H2 = (mass_flow * H1 + EH + cell_calib['Qs_calib']) / mass_flow
                    H3 = hep.calculate(9,1,P1*1e5,2,T3)
                    # rho_avg = hep.calculate(3,1,P1*1e5/2,9, (H2+H3)/2) ## BUG
                    # mu_avg = hep.calculate(25,1,P1*1e5/2,9, (H2+H3)/2) ## BUG
                    rho_avg = hep.calculate(3,1,P1*1e5,9, (H2+H3)/2)
                    mu_avg = hep.calculate(25,1,P1*1e5,9, (H2+H3)/2)
                    P3 = P1 - pycryo.pressureDrop(m=mass_flow/cell_calib['n_channels_tot'],
                                                  mu=mu_avg, rho=rho_avg, D=2*cell_calib['channel_radius'],
                                                  L=cell_calib['length'], rug=cell_calib['roughness'])
                    H3 = hep.calculate(9,1,P3*1e5,2,T3)#re-evaluate H3 with P3
                    # # QBS_temp = (FT1+FT2) * (H3 - H1) - QS - EH

                    #### pressure drop
                    #####

                else:
                    mass_flow = other['mass_flow']
                    P3 = other['P3']
                    x = hlr.massflowDistrib(mass_flow, C1TT, C2TT, P1, 
                                            hc_type, 2*cell_calib["channel_radius"], cell_calib["roughness"])
                    FT1 = mass_flow * x
                    FT2 = mass_flow * (1. - x)

                total_length = sum(instrumented_cells_config[circuit]['magnet_lengths'])
                dP = (P1-P3)/total_length

                qbs_magnets_circuits = [[],[]]
                ## only if run3
                run3 = True
                # from cryo
                # LD = 15
                # LQ = 8
                # Lmag1 = LQ ##BUG
                Lmag1 = instrumented_cells_config[circuit]['magnet_lengths'][0] 
                LD = instrumented_cells_config[circuit]['magnet_lengths'][1]
                #First magnet 
                j=0
                qbs_mag = hlr.hlr_mag(C1TT[j],C1TT[j],C1TT[j+1],P1,P1-dP*Lmag1,FT1)
                if j in faults1 or j+1 in faults1: qbs_mag *= 0
                qbs_magnets_circuits[0].append( qbs_mag ) #beam1
                qbs_mag = hlr.hlr_mag(C2TT[j],C2TT[j],C2TT[j+1],P1,P1-dP*Lmag1,FT2)
                # qbs_mag = hlr.hlr_mag(C2TT[j],C2TT[j],C2TT[j+1],P1,P1-dP+Lmag1,FT2) ## BUG
                if j in faults2 or j+1 in faults2: qbs_mag *= 0
                qbs_magnets_circuits[1].append( qbs_mag ) #beam2
                #Second magnet
                j=1
                qbs_mag = hlr.hlr_mag(C2TT[j],C2TT[j],C2TT[j+1],P1-dP*Lmag1,P1-dP*(Lmag1+LD),FT2)
                if j in faults2 or j+1 in faults2: qbs_mag *= 0
                qbs_magnets_circuits[1].append( qbs_mag ) #beam1
                qbs_mag = hlr.hlr_mag(C1TT[j],C1TT[j],C1TT[j+1],P1-dP*Lmag1,P1-dP*(Lmag1+LD),FT1)
                if j in faults1 or j+1 in faults1: qbs_mag *= 0
                qbs_magnets_circuits[0].append( qbs_mag ) #beam2
                #Third magnet  
                j=2
                qbs_mag = hlr.hlr_mag(C1TT[j],C1TT[j],C1TT[j+1],P1-dP*(Lmag1+LD),P1-dP*(Lmag1+2*LD),FT1)
                if j in faults1 or j+1 in faults1: qbs_mag *= 0
                qbs_magnets_circuits[0].append( qbs_mag ) #beam1
                qbs_mag = hlr.hlr_mag(C2TT[j],C2TT[j],C2TT[j+1],P1-dP*(Lmag1+LD),P1-dP*(Lmag1+2*LD),FT2)
                if j in faults2 or j+1 in faults2: qbs_mag *= 0
                qbs_magnets_circuits[1].append( qbs_mag ) #beam2
                #Fourth magnet
                j=3
                if run3:
                    qbs_mag = hlr.hlr_mag(C2TT[j],C2TT[j],C2TT[j+1],P1-dP*(Lmag1+2*LD),P3,FT2)
                    if j in faults2 or j+1 in faults2: qbs_mag *= 0
                    qbs_magnets_circuits[1].append( qbs_mag ) #beam1
                    qbs_mag = hlr.hlr_mag(C1TT[j],C1TT[j],C1TT[j+1],P1-dP*(Lmag1+2*LD),P3,FT1)
                    if j in faults1 or j+1 in faults1: qbs_mag *= 0
                    qbs_magnets_circuits[0].append( qbs_mag ) #beam2
                else: # only one temperature sensor so average between two??
                    ### UNTESTED
                    qbs_mag = hlr.hlr_mag(C2TT[j],C1TT[j],C2TT[j+1],P1-dP*(Lmag1+2*LD),P3,FT2)
                    if j in faults1 or j in faults2 or j+1 in faults2: qbs_mag *= 0
                    qbs_magnets_circuits[1].append( qbs_mag ) #beam1
                    qbs_mag = hlr.hlr_mag(C1TT[j],C2TT[j],C1TT[j+1],P1-dP*(Lmag1+2*LD),P3,FT1)
                    if j in faults1 or j in faults2 or j+1 in faults1: qbs_mag *= 0
                    qbs_magnets_circuits[0].append( qbs_mag ) #beam2

                #old run2
                ## t_out_magnets_circuits = [[obraw.dictionary[vv] for vv in
                ##     out_sensor_names_circuits[ii]] for ii in [0, 1]]
                ## t_in_magnets_circuits = [[obraw.dictionary[vv] for vv in
                ##     in_sensor_names_circuits[ii]] for ii in [0, 1]]

                ## qbs_magnets_circuits, other_instr = \
                ##     hlr.compute_heat_loads_instrumented_cell(
                ##         mass_flow = other['mass_flow'], P1=P1,
                ##         T_in_magnets_circuits=t_in_magnets_circuits,
                ##         T_out_magnets_circuits=t_out_magnets_circuits,
                ##         magnet_lengths_circuits=magnet_lengths_circuits,
                ##         n_channels_circuits=n_channels_circuits,
                ##         channel_radius=cell_calib['channel_radius'],
                ##         channel_roughness=cell_calib['roughness'],
                ##         dp_toll = 0.001, N_iter_max=200)

                ## unused
                ## magnet_beam_circuits = [
                ##         instrum_cell_config['circuit_%s_beam'%cc]
                ##                         for cc in ['A', 'A']]
            
                # ## rename magnets with _COR when using coriolis
                # if coriolis:
                #     instrum_cell_config_cor = copy.deepcopy(instrum_cell_config)
                #     magnet_names = instrum_cell_config_cor['magnet_names']
                #     for ii in range(len(magnet_names)):
                #         magnet_names[ii] = f"{magnet_names[ii]}_COR"
                #     instrum_cell_config_cor['magnet_names'] = magnet_names
                # # else:

                dict_output = hlr.build_instrumented_hl_dict(
                    config_dict=instrum_cell_config, circuit=circuit,
                     Qbs_magnets_circuits=qbs_magnets_circuits, coriolis=coriolis)
                instrum_cell_recalc_dict.update(dict_output)
                # print(dict_output.keys())

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
            variables=(circuits + avg_varnames + instrum_varnames))

    other = {}
    other['issues'] = issues

    return obhl_store, other

