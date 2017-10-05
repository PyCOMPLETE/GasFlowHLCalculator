# Note: There are different naming conventions for the special cells.
# For example 12R4 or 13R4.
from __future__ import division
import numpy as np

import LHCMeasurementTools.TimberManager as tm

from Helium_properties import interp_P_T_hPT, interp_P_T_DPT
from data_S45_details import cell_timber_vars_dict, cell_list, cell_list_pre_EYETS16
import compute_QBS_magnet as cqm
from valve_LT import valve_LT

zeros = lambda *x: np.zeros(shape=x, dtype=float)

# Also the order of magnets in a cell with non reversed gasflow
magnet_ids = ['Q1', 'D2', 'D3', 'D4']

def mass_flow(atd, new_cell):
    """
    Returns at object which contains the necessary information on the cell mass flow and pressure.
    """
    n_tt = len(atd.timestamps)

    # Account for missing cell
    if new_cell:
        cells = cell_list
    else:
        cells = cell_list_pre_EYETS16
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

    mag_calc = cqm.QbsMagnetCalculator(interp_P_T_hPT, atd, P1, m_L, cell_list)
    return mag_calc, Qbs

def make_dict(mag_calc, Qbs, atd, new_cell):
    qbs_special = {}
    qbs_special['timestamps'] = atd.timestamps
    if new_cell:
        cells = cell_list
    else:
        cells = cell_list_pre_EYETS16
    qbs_special['cells'] = cells

    for cell_ctr, cell in enumerate(cells):
        qbs_special[cell] = {}
        sum_all_magnet_hl = 0
        for magnet_id in magnet_ids:
            magnet_hl = mag_calc.Compute_QBS_magnet(cell, cell_timber_vars_dict[cell][magnet_id]['Tin'], cell_timber_vars_dict[cell][magnet_id]['Tout'])
            sum_all_magnet_hl += magnet_hl
            qbs_special[cell][magnet_id] = magnet_hl
        qbs_special[cell]['Sum'] = sum_all_magnet_hl
        qbs_special[cell]['Qbs'] = Qbs[:,cell_ctr]

    return qbs_special

def make_dict_separate(mag_calc, Qbs, atd, new_cell):
    qbs_special = {}
    qbs_special['timestamps'] = atd.timestamps
    if new_cell:
        cells = cell_list
    else:
        cells = cell_list_pre_EYETS16
    qbs_special['cells'] = cells

    for cell_ctr, cell in enumerate(cells):
        qbs_special[cell] = {}
        sum_all_magnet_hl = np.zeros_like(atd.timestamps)
        for magnet_id in magnet_ids:
            qbs_special[cell][magnet_id] = {}
            sum_magnet_hl = np.zeros_like(atd.timestamps)
            for beam_number in (1,2):
                Tin, Tout = _get_Tin_Tout(cell, magnet_id, beam_number)
                try:
                    magnet_hl = mag_calc.Compute_QBS_magnet_single(cell, Tin, Tout)/2.
                except cqm.InvalidDataError:
                    magnet_hl = np.empty_like(atd.timestamps)*np.nan
                else:
                    sum_all_magnet_hl += magnet_hl
                    sum_magnet_hl += magnet_hl

                qbs_special[cell][magnet_id][beam_number] = magnet_hl
            if np.all(sum_magnet_hl == 0):
                sum_magnet_hl = np.empty_like(sum_magnet_hl)*np.nan
            qbs_special[cell][magnet_id]['Sum'] = sum_magnet_hl

        if np.all(sum_all_magnet_hl == 0):
            sum_magnet_hl = np.empty_like(sum_all_magnet_hl)*np.nan
        qbs_special[cell]['Sum'] = sum_all_magnet_hl
        qbs_special[cell]['Qbs'] = Qbs[:,cell_ctr]
    return qbs_special


def _get_Tin_Tout(cell, magnet, beam):
    """
    Returns the variable names of temperature sensors for a given cell, magnet and beam.
    """
    # For cell with normal gas flow
    dict_beam_Tin_normal = {
            1: 'TT826',
            2: 'TT824',
        }

    dict_beam_Tout_normal = {
            1: 'TT824',
            2: 'TT826',
        }

    # reversed
    dict_beam_Tin_reversed = dict_beam_Tout_normal
    dict_beam_Tout_reversed = dict_beam_Tin_normal


    if cell == '13L5':
        dict_beam_Tin = dict_beam_Tin_reversed
        dict_beam_Tout = dict_beam_Tout_reversed
    else:
        dict_beam_Tin = dict_beam_Tin_normal
        dict_beam_Tout = dict_beam_Tout_normal


    dd = cell_timber_vars_dict[cell]
    list_Tin = dd[magnet]['Tin']
    list_Tout = dd[magnet]['Tout']
    if len(list_Tin) != 1:
        identifier = dict_beam_Tin[beam]
        list_Tin = filter(lambda x: identifier in x, list_Tin)

    if len(list_Tout) != 1:
        identifier = dict_beam_Tout[beam]
        list_Tout = filter(lambda x: identifier in x, list_Tout)

    assert len(list_Tout) == 1
    assert len(list_Tin) == 1

    return list_Tin[0], list_Tout[0]

def compute_qbs_special(atd, new_cell, separate=True, aligned=False):
    mag_calc, Qbs = mass_flow(atd, new_cell)

    if separate:
        dict_ = make_dict_separate(mag_calc, Qbs, atd, new_cell)
        if aligned:
            return dict_to_aligned_separate(dict_)
        else:
            return dict_
    else:
        dict_ = make_dict(mag_calc, Qbs, atd, new_cell)
        if aligned:
            return dict_to_aligned(dict_)
        else:
            return dict_


def dict_to_aligned(dict_):
    timestamps = dict_['timestamps']
    variables = []
    data = []
    for cell in dict_['cells']:
        dd = dict_[cell]
        for key, arr in dd.iteritems():
            main_key = cell + '_' + key
            variables.append(main_key)
            data.append(arr)

    data_arr = np.array(data).T
    return tm.AlignedTimberData(timestamps, data_arr, np.array(variables))

def dict_to_aligned_separate(dict_):
    timestamps = dict_['timestamps']
    variables = []
    data = []
    for cell in dict_['cells']:
        dd = dict_[cell]
        for magnet_id in magnet_ids:
            mag_dict = dd[magnet_id]
            for beam_id, arr in mag_dict.iteritems():
                if beam_id == 'Sum':
                    main_key = '_'.join([cell, magnet_id])
                else:
                    main_key = '_'.join([cell, magnet_id, str(beam_id)])
                variables.append(main_key)
                data.append(arr)
        for key in ('Sum', 'Qbs'):
            main_key = '_'.join([cell, key])
            variables.append(main_key)
            data.append(dd[key])

    data_arr = np.array(data, dtype=float).T
    return tm.AlignedTimberData(timestamps, data_arr, np.array(variables))


def aligned_to_dict(qbs_ob):
    output = {}
    output['timestamps'] = qbs_ob.timestamps
    output['cells'] = cell_list
    for cell in cell_list:
        output[cell] = dd = {}
        for key in qbs_ob.variables:
            if cell in key:
                subkey = key.split(cell+'_')[1]
                dd[subkey] = qbs_ob.dictionary[key]
    return output

def aligned_to_dict_separate(qbs_ob):
    output = {}
    output['timestamps'] = qbs_ob.timestamps
    cells = set()

    # reverse dict to aligned seperate
    for var, arr in qbs_ob.dictionary.iteritems():
        split_var = var.split('_')
        if len(split_var) == 3:
            cell, magnet_id, beam = split_var

            if cell not in output:
                output[cell] = {}
            if magnet_id not in output[cell]:
                output[cell][magnet_id] = {}

            output[cell][magnet_id][int(beam)] = arr
        elif len(split_var) == 2:
            cell, key = split_var
            if cell not in output:
                output[cell] = {}

            if key in ('Qbs', 'Sum'):
                output[cell][key] = arr
            elif key in magnet_ids:
                if key not  in output[cell]:
                    output[cell][key] = {}

                output[cell][key]['Sum'] = arr
            else:
                raise ValueError(key)
            cells.add(cell)

    if set(cell_list) == cells:
        cells = cell_list
    elif set(cell_list_pre_EYETS16) == cells:
        cells = cell_list_pre_EYETS16

    output['cells'] = cells
    return output

