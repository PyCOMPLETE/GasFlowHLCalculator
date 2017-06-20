# Note: There are different naming conventions for the special cells.
# For exampsle 12R4 or 13R4.
import numpy as np

from Helium_properties import interp_P_T_hPT, interp_P_T_DPT
from data_S45_details import *
from compute_QBS_magnet import QbsMagnetCalculator
from valve_LT import valve_LT

zeros = lambda *x: np.zeros(shape=(x), dtype=float)

cell_list = ['13R4', '33L5', '13L5', '31L2']
cell_list_old = ['13R4', '33L5', '13L5']

def mass_flow(atd):
    n_tt = len(atd.timestamps)
    n_list = len(TT94x_list)

    m_L = zeros(n_tt, n_list) # mass flow
    T1 = zeros(n_tt, n_list) #TT961
    T3 = zeros(n_tt, n_list) #TT94x
    CV = zeros(n_tt, n_list) #CV94x
    EH = zeros(n_tt, n_list) #EH84x
    P1 = zeros(n_tt, n_list) #PT961
    P4 = zeros(n_tt, n_list) #PT991
    hC = zeros(n_tt, n_list) #enthalpy of line C
    h3 = zeros(n_tt, n_list) #enthalpy of BS output
    ro = zeros(n_tt, n_list) #density
    Qbs = zeros(n_tt, n_list)

    arr_list = [T1, T3, CV, EH, P1, P4]
    name_list = [TT961_list, TT94x_list, CV94x_list, EH84x_list, PT961_list, PT991_list]
    varlist = list(atd.variables)

    for names, arr in zip(name_list, arr_list):
        for ii, var in enumerate(names):
            index = varlist.index(var)
            arr[:,ii] = atd.data[:,index]

    for i in xrange(n_list):
        for j in xrange(n_tt):
            hC[j,i] = interp_P_T_hPT(P1[j,i],T1[j,i])
            h3[j,i] = interp_P_T_hPT(P1[j,i],T3[j,i])
            ro[j,i] = interp_P_T_DPT(P1[j,i],T3[j,i])
        m_L[:,i] = valve_LT(P1[:,i], P4[:,i], ro[:,i], Kv_list[i], CV[:,i], R_list[i])
        Qbs[:,i] = m_L[:,i]*(h3[:,i] - hC[:,i]) - Qs_list[i] - EH[:,i]

    Compute_QBS_magnet = QbsMagnetCalculator(interp_P_T_hPT, atd, P1, m_L).Compute_QBS_magnet
    return Compute_QBS_magnet, Qbs

def make_dict(Compute_QBS_magnet, Qbs, atd, new_cell):
    qbs_special = {}
    qbs_special['timestamps'] = atd.timestamps
    if new_cell:
        qbs_special['cells'] = cell_list
    else:
        qbs_special['cells'] = cell_list_old

    #compute each magnet QBS
    QBS_Q1_12R4 = Compute_QBS_magnet(0,Q1_Tin_12R4,Q1_Tout_12R4)
    QBS_D2_12R4 = Compute_QBS_magnet(0,D2_Tin_12R4,D2_Tout_12R4)
    QBS_D3_12R4 = Compute_QBS_magnet(0,D3_Tin_12R4,D3_Tout_12R4)
    QBS_D4_12R4 = Compute_QBS_magnet(0,D4_Tin_12R4,D4_Tout_12R4)
    QBS_12R4_sum = QBS_Q1_12R4 + QBS_D2_12R4 + QBS_D3_12R4 + QBS_D4_12R4
    # Be careful of naming conventions for cells!
    qbs_special['13L5'] = {
            'Q1': QBS_Q1_12R4,
            'D2': QBS_D2_12R4,
            'D3': QBS_D3_12R4,
            'D4': QBS_D4_12R4,
            'Sum': QBS_12R4_sum,
            'qbs': Qbs[:,0],
            }

    # This is the cell with the faulty temperature sensor
    QBS_Q1_32R4 = Compute_QBS_magnet(1,Q1_Tin_32R4,Q1_Tout_32R4)
    QBS_D2_32R4 = Compute_QBS_magnet(1,D2_Tin_32R4,D2_Tout_32R4)
    QBS_D3_32R4 = Compute_QBS_magnet(1,D3_Tin_32R4,D3_Tout_32R4)
    QBS_D4_32R4 = Compute_QBS_magnet(1,D4_Tin_32R4,D4_Tout_32R4)
    QBS_32R4_sum = QBS_Q1_32R4 + QBS_D2_32R4 + QBS_D3_32R4 + QBS_D4_32R4
    # Be careful of naming conventions for cells!
    qbs_special['33L5'] = {
            'Q1': QBS_Q1_32R4,
            'D2': QBS_D2_32R4,
            'D3': QBS_D3_32R4,
            'D4': QBS_D4_32R4,
            'Sum': QBS_32R4_sum,
            'qbs': Qbs[:,1],
            }

    # This is the reversed cell
    QBS_Q1_13L5 = Compute_QBS_magnet(2,Q1_Tin_13L5,Q1_Tout_13L5)
    QBS_D2_13L5 = Compute_QBS_magnet(2,D2_Tin_13L5,D2_Tout_13L5)
    QBS_D3_13L5 = Compute_QBS_magnet(2,D3_Tin_13L5,D3_Tout_13L5)
    QBS_D4_13L5 = Compute_QBS_magnet(2,D4_Tin_13L5,D4_Tout_13L5)
    QBS_13L5_sum = QBS_Q1_13L5 + QBS_D2_13L5 + QBS_D3_13L5 + QBS_D4_13L5
    # Be careful of naming conventions for cells!
    qbs_special['13R4'] = {
            'Q1': QBS_Q1_13L5,
            'D2': QBS_D2_13L5,
            'D3': QBS_D3_13L5,
            'D4': QBS_D4_13L5,
            'Sum': QBS_13L5_sum,
            'Qbs': Qbs[:,2],
            }

    # This is the new cell
    if new_cell:
        QBS_Q1_32L2 = Compute_QBS_magnet(3,Q1_Tin_32L2,Q1_Tout_32L2)
        QBS_D2_32L2 = Compute_QBS_magnet(3,D2_Tin_32L2,D2_Tout_32L2)
        QBS_D3_32L2 = Compute_QBS_magnet(3,D3_Tin_32L2,D3_Tout_32L2)
        QBS_D4_32L2 = Compute_QBS_magnet(3,D4_Tin_32L2,D4_Tout_32L2)
        QBS_32L2_sum = QBS_Q1_32L2 + QBS_D2_32L2 + QBS_D3_32L2 + QBS_D4_32L2
        # Be careful of naming conventions for cells!
        qbs_special['31L2'] = {
                'Q1': QBS_Q1_32L2,
                'D2': QBS_D2_32L2,
                'D3': QBS_D3_32L2,
                'D4': QBS_D4_32L2,
                'Sum': QBS_32L2_sum,
                'Qbs': Qbs[:,3],
                }

    return qbs_special

def make_dict_separate(Compute_QBS_magnet, Qbs, atd, qbs_special):
    raise ValueError('This does not currently work')
    top, bot = 0, 1 #826 and 824 temp sensors

    #compute each magnet QBS

    QBS_Q1_12R4_b2 = Compute_QBS_magnet(0,Q1_Tin_12R4,Q1_Tout_12R4[top])
    QBS_Q1_12R4_b1 = Compute_QBS_magnet(0,Q1_Tin_12R4,Q1_Tout_12R4[bot])
    QBS_D2_12R4_b2 = Compute_QBS_magnet(0,D2_Tin_12R4[bot],D2_Tout_12R4[top])
    QBS_D2_12R4_b1 = Compute_QBS_magnet(0,D2_Tin_12R4[top],D2_Tout_12R4[bot])
    QBS_D3_12R4_b2 = Compute_QBS_magnet(0,D3_Tin_12R4[bot],D3_Tout_12R4[top])
    QBS_D3_12R4_b1 = Compute_QBS_magnet(0,D3_Tin_12R4[top],D3_Tout_12R4[bot])

    # Be careful of naming conventions for cells!
    qbs_special['13L5'].update({
            'Q1_b2': QBS_Q1_12R4_b2,
            'Q1_b1': QBS_Q1_12R4_b1,
            'D2_b2': QBS_D2_12R4_b2,
            'D2_b1': QBS_D2_12R4_b1,
            'D3_b1': QBS_D3_12R4_b1,
            'D3_b2': QBS_D3_12R4_b2,
            })

    # This is the cell with the faulty temperature sensor
    QBS_Q1_32R4_b2 = Compute_QBS_magnet(1,Q1_Tin_32R4,Q1_Tout_32R4[top])
    QBS_Q1_32R4_b1 = Compute_QBS_magnet(1,Q1_Tin_32R4,Q1_Tout_32R4[bot])
    QBS_D2_32R4_b2 = Compute_QBS_magnet(1,D2_Tin_32R4[bot],D2_Tout_32R4[top])
    QBS_D2_32R4_b1 = Compute_QBS_magnet(1,D2_Tin_32R4[top],D2_Tout_32R4[bot])
    QBS_D3_32R4_b1 = Compute_QBS_magnet(1,D3_Tin_32R4[top],D3_Tout_32R4)
    # Be careful of naming conventions for cells!
    qbs_special['33L5'].update({
            'Q1_b1': QBS_Q1_32R4_b1,
            'Q1_b2': QBS_Q1_32R4_b2,
            'D2_b1': QBS_D2_32R4_b1,
            'D2_b2': QBS_D2_32R4_b2,
            'D3_b1': QBS_D3_32R4_b1,
            })

    # This is the reversed cell
    QBS_D4_13L5_b1 = Compute_QBS_magnet(2,D4_Tin_13L5,D4_Tout_13L5[top])
    QBS_D4_13L5_b2 = Compute_QBS_magnet(2,D4_Tin_13L5,D4_Tout_13L5[bot])
    QBS_D3_13L5_b1 = Compute_QBS_magnet(2,D3_Tin_13L5[bot],D3_Tout_13L5[top])
    QBS_D3_13L5_b2 = Compute_QBS_magnet(2,D3_Tin_13L5[top],D3_Tout_13L5[bot])
    QBS_D2_13L5_b1 = Compute_QBS_magnet(2,D2_Tin_13L5[bot],D2_Tout_13L5[top])
    QBS_D2_13L5_b2 = Compute_QBS_magnet(2,D2_Tin_13L5[top],D2_Tout_13L5[bot])
    # Be careful of naming conventions for cells!
    qbs_special['13R4'].update({
            'D2_b1': QBS_D2_13L5_b1,
            'D2_b2': QBS_D2_13L5_b2,
            'D3_b1': QBS_D3_13L5_b1,
            'D3_b2': QBS_D3_13L5_b2,
            'D4_b1': QBS_D4_13L5_b1,
            'D4_b2': QBS_D4_13L5_b2,
            })

    return qbs_special

def compute_qbs_special(atd, new_cell, separate=False):
    Compute_QBS_magnet, Qbs = mass_flow(atd)
    special_dict = make_dict(Compute_QBS_magnet, Qbs, atd, new_cell)
    if not separate:
        return special_dict
    else:
        return make_dict_separate(Compute_QBS_magnet, Qbs, atd, special_dict)


if __name__ == '__main__':
    import os
    import matplotlib.pyplot as plt
    import LHCMeasurementTools.TimberManager as tm
    from LHCMeasurementTools.SetOfHomogeneousVariables import SetOfHomogeneousNumericVariables as shnv
    plt.close('all')
    dt_seconds = 30
    filename = os.path.dirname(os.path.abspath(__file__)) + '/TIMBER_DATA_special_5030.csv'
    tv = tm.parse_timber_file(filename)
    atd = shnv(tv.keys(), tv).aligned_object(dt_seconds)
    qbs_special = compute_qbs_special(atd)
    tt = (qbs_special['timestamps'] - qbs_special['timestamps'][0]) / 3600.

    fig = plt.figure()
    for ctr, cell in enumerate(qbs_special['cells']):
        sp = plt.subplot(2,2,ctr+1)
        hl_dict = qbs_special[cell]
        for key in hl_dict:
            sp.plot(tt,hl_dict[key], label=key)
        if ctr == 1:
            sp.legend(bbox_to_anchor=(1.1,1))
        sp.set_xlabel('Time [h]')
        sp.set_ylabel('Qdbs [W]')
        sp.set_title(cell)
    plt.show()
