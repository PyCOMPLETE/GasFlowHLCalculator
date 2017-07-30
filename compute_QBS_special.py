# Note: There are different naming conventions for the special cells.
# For exampsle 12R4 or 13R4.
from __future__ import division
import numpy as np

from Helium_properties import interp_P_T_hPT, interp_P_T_DPT
from data_S45_details import cell_timber_vars_dict
from compute_QBS_magnet import QbsMagnetCalculator
from valve_LT import valve_LT

zeros = lambda *x: np.zeros(shape=x, dtype=float)

cell_list = ['13R4', '33L5', '13L5', '31L2']
cell_list_old = ['13R4', '33L5', '13L5']

# Also the order of magnets in a cell with non reversed gasflow
magnet_ids = ['Q1', 'D2', 'D3', 'D4']

def mass_flow(atd, new_cell):
    n_tt = len(atd.timestamps)

    # Account for missing cell
    if new_cell:
        cells = cell_list
    else:
        cells = cell_list_old
    n_cells = len(cells)

    m_L = zeros(n_tt, n_cells) # mass flow
    T1 = zeros(n_tt, n_cells) #TT961
    T3 = zeros(n_tt, n_cells) #TT94x
    CV = zeros(n_tt, n_cells) #CV94x
    EH = zeros(n_tt, n_cells) #EH84x
    P1 = zeros(n_tt, n_cells) #PT961
    P4 = zeros(n_tt, n_cells) #PT991
    hC = zeros(n_tt, n_cells) #enthalpy of line C
    h3 = zeros(n_tt, n_cells) #enthalpy of BS output
    ro = zeros(n_tt, n_cells) #density
    Qbs = zeros(n_tt, n_cells)

    arr_list = [T1, T3, CV, EH, P1, P4]
    name_list = ['TT961', 'TT94x', 'CV94x', 'EH84x', 'PT961', 'PT991']
    varlist = list(atd.variables)

    for name, arr in zip(name_list, arr_list):
        for ii, cell in enumerate(cells):
            var = cell_timber_vars_dict[cell][name]
            index = varlist.index(var)
            arr[:,ii] = atd.data[:,index]

    for i, cell in enumerate(cells):
        Kv = cell_timber_vars_dict[cell]['Kv']
        Qs = cell_timber_vars_dict[cell]['Qs']
        R = cell_timber_vars_dict[cell]['R']
        for j in xrange(n_tt):
            hC[j,i] = interp_P_T_hPT(P1[j,i],T1[j,i])
            h3[j,i] = interp_P_T_hPT(P1[j,i],T3[j,i])
            ro[j,i] = interp_P_T_DPT(P1[j,i],T3[j,i])
            m_L[:,i] = valve_LT(P1[:,i], P4[:,i], ro[:,i], Kv, CV[:,i], R)
        Qbs[:,i] = m_L[:,i]*(h3[:,i] - hC[:,i]) - Qs - EH[:,i]

    Compute_QBS_magnet = QbsMagnetCalculator(interp_P_T_hPT, atd, P1, m_L, cell_list).Compute_QBS_magnet
    return Compute_QBS_magnet, Qbs

def make_dict(Compute_QBS_magnet, Qbs, atd, new_cell):
    qbs_special = {}
    qbs_special['timestamps'] = atd.timestamps
    if new_cell:
        cells = cell_list
    else:
        cells = cell_list_old
    qbs_special['cells'] = cells

    for cell_ctr, cell in enumerate(cells):
        qbs_special[cell] = {}
        sum_magnet_hl = 0
        for magnet_id in magnet_ids:
            magnet_hl = Compute_QBS_magnet(cell, cell_timber_vars_dict[cell][magnet_id]['Tin'], cell_timber_vars_dict[cell][magnet_id]['Tout'])
            sum_magnet_hl += magnet_hl
            qbs_special[cell][magnet_id] = magnet_hl
        qbs_special[cell]['Sum'] = sum_magnet_hl
        qbs_special[cell]['qbs'] = Qbs[:,cell_ctr]

    return qbs_special

def make_dict_separate(Compute_QBS_magnet, Qbs, atd, new_cell):
    qbs_special = {}
    qbs_special['timestamps'] = atd.timestamps
    if new_cell:
        cells = cell_list
    else:
        cells = cell_list_old
    qbs_special['cells'] = cells

    for cell_ctr, cell in enumerate(cells):
        qbs_special[cell] = {}
        sum_magnet_hl = 0
        for magnet_id in magnet_ids:
            qbs_special[cell][magnet_id] = {}
            for beam_number in (1,2):
                Tin, Tout = _get_Tin_Tout(cell, magnet_id, beam_number)
                magnet_hl = Compute_QBS_magnet(cell, Tin, Tout)/2.
                sum_magnet_hl += magnet_hl
                qbs_special[cell][magnet_id][beam_number] = magnet_hl
        qbs_special[cell]['Sum'] = sum_magnet_hl
        qbs_special[cell]['qbs'] = Qbs[:,cell_ctr]
    return qbs_special


def _get_Tin_Tout(cell, magnet, beam):

    # for cell with the gas flow in normal direction
    dict_beam_Tin = {
        1: 'TT826',
        2: 'TT824',
    }

    dict_beam_Tout = {
        1: 'TT824',
        2: 'TT826',
    }

    dd = cell_timber_vars_dict[cell]
    list_Tin = dd[magnet]['Tin']
    list_Tout = dd[magnet]['Tout']
    if len(list_Tin) == 1:
        Tin = list_Tin[0]
    else:
        identifier = dict_beam_Tin[beam]
        Tin = filter(lambda x: identifier in x, list_Tin)[0]

    if len(list_Tout) == 1:
        Tout = list_Tout[0]
    else:
        identifier = dict_beam_Tout[beam]
        Tout = filter(lambda x: identifier in x, list_Tout)[0]

    return Tin, Tout

def compute_qbs_special(atd, new_cell, separate=False):
    Compute_QBS_magnet, Qbs = mass_flow(atd, new_cell)
    if not separate:
        return make_dict(Compute_QBS_magnet, Qbs, atd, new_cell)
    else:
       return make_dict_separate(Compute_QBS_magnet, Qbs, atd, new_cell)
    raise ValueError('Currently not implemented')


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
