from __future__ import division

import numpy as np
from scipy.interpolate import interp2d

import LHCMeasurementTools.TimberManager as tm

from Helium_properties import P, T, D_PT, h_PT, mu_PT
from valve_LT import valve_LT, valve_LT_arr
from Pressure_drop import Pressure_drop
from data_qbs import Data_qbs, data_qbs
import h5_storage

zeros = lambda *x: np.zeros(shape=(x), dtype=float)
dq_latest = data_qbs
max_iterations = 5

interp_P_T_hPT = interp2d(P,T,h_PT)
interp_P_T_DPT = interp2d(P,T,D_PT)
interp_P_T_mu = interp2d(P,T,mu_PT)

def value_getter_factory(missing_variables, corrected_variables, atd_ob, nan_arr, nan_variables):
    def get_value(arr, names, correct_zero):
        correct_first = False
        for cell_ctr, name in enumerate(names):
            try:
                data = atd_ob.dictionary[name]
            except KeyError:
                missing_variables.add(name)
            else:
                arr[:,cell_ctr] = data
                if arr[0,cell_ctr] == 0:
                    if correct_zero:
                        if cell_ctr > 0:
                            arr[:,cell_ctr] = arr[:,cell_ctr-1]
                            corrected_variables.add(name)
                        else:
                            correct_first = True
                    else:
                        nan_variables.add(name)
                        nan_arr[cell_ctr] = True
        return correct_first
    return get_value
        self.nan_arr = np.zeros_like(atd_ob.data, dtype=np.bool)
        self.

def compute_qbs(atd_ob, use_dP, version=h5_storage.version, strict=True):

    if version != h5_storage.version:
        dq = Data_qbs(version)
    else:
        dq = dq_latest

    Ncell = len(dq.Cell_list) # number of half cells
    Nvalue = len(atd_ob.timestamps) # number of data points in time

    T1 = zeros(Nvalue, Ncell) #TT961
    T3 = zeros(Nvalue, Ncell) #TT94x
    CV = zeros(Nvalue, Ncell) #CV94x;
    EH = zeros(Nvalue, Ncell) #EH84x;
    P1 = zeros(Nvalue, Ncell) #PT961;
    P4 = zeros(Nvalue, Ncell) #PT991;
    T2 = zeros(Nvalue, Ncell) #TT84x

    var_data_dict = {
        'T1': {
            'data': T1,
            'vars': dq.TT961_list,
            'correct': True,
        },
        'T3': {
            'data': T3,
            'vars': dq.TT94x_list,
            'correct': False,
        },
        'CV': {
            'data': CV,
            'vars': dq.CV94x_list,
            'correct': False,
        },
        'EH': {
            'data': EH,
            'vars': dq.EH84x_list,
            'correct': False,
        },
        'P1': {
            'data': P1,
            'vars': dq.PT961_list,
            'correct': True,
        },
        'P4': {
            'data': P4,
            'vars': dq.PT991_list,
            'correct': True,
        },
        'T2': {
            'data': T2,
            'vars': dq.TT84x_list,
            'correct': False,
        },
    }

    missing_variables = set()
    corrected_variables = set()
    nan_variables = set()
    nan_arr = np.zeros(Ncell, dtype=np.bool)
    get_value = value_getter_factory(missing_variables, corrected_variables, atd_ob, nan_arr, nan_variables)

    for dd in var_data_dict.itervalues():
        correct_first = get_value(dd['data'], dd['vars'], dd['correct'])
        if correct_first:
            dd['data'][:,0] = dd['data'][:,-1]

    if missing_variables:
        print('Missing variables:', missing_variables)
        if strict:
            raise ValueError('There have been missing variables!')
    if corrected_variables:
        print('Corrected variables:', corrected_variables)
    if nan_variables:
        print('Nan variables:', nan_variables)

    ro          = zeros(Nvalue,Ncell)    #density
    m_L         = zeros(Nvalue,Ncell)    #massflow
    qbs         = zeros(Nvalue,Ncell)    #BS heat load per half-cell
    hC          = zeros(Nvalue,Ncell)    #enthalpy of line C
    h3          = zeros(Nvalue,Ncell)    #enthalpy of BS output
    counter_int = np.zeros(shape=(Nvalue,Ncell),dtype=int)  #internal counter

    #interp_P_T_g = interp2d(P,T,gamma_PT) # not used

    P3 = np.copy(P1) #intermediate pressure ater the beam screen and before the valve
    rug = dq.rug
    Radius = dq.Radius
    R_list = dq.R_list
    Kv_list = dq.Kv_list
    nc_list = dq.nc_list
    L_list = dq.L_list
    Qs_list = dq.Qs_list

    for i, isnan in enumerate(nan_arr):
        if isnan:
            continue
        for j in xrange(Nvalue):
            hC[j,i] = interp_P_T_hPT(P1[j,i],T1[j,i])
            h3[j,i] = interp_P_T_hPT(P1[j,i],T3[j,i])
            ro[j,i] = interp_P_T_DPT(P1[j,i],T3[j,i])
            #gamma[j,i] = interp_P_T_g(P1[j,i],T3[j,i]) # not used

        #compute the intermediate pressure P3 by iteration
            if use_dP:
                P3[j,i], counter_int[j,i] = compute_P3(P1[j,i], P4[j,i], T2[j,i], T3[j,i], Kv_list[i], CV[j,i], R_list[i], ro[j,i], nc_list[i], Radius, rug, L_list[i])

        m_L[:,i] = valve_LT_arr(P3[:,i],P4[:,i],ro[:,i],Kv_list[i],CV[:,i],R_list[i])
        qbs[:,i] = m_L[:,i]*(h3[:,i]-hC[:,i])-Qs_list[i]-EH[:,i]

    n_max_iterations = np.sum(counter_int == max_iterations)
    if n_max_iterations != 0:
        print('Warning: Maximum iterations were reached %i times!' % n_max_iterations)

    # Nan values
    qbs[:,nan_arr] = np.nan
    return tm.AlignedTimberData(atd_ob.timestamps, qbs, dq.Cell_list)

def compute_P3(P1, P4, T2, T3, Kv, CV, R, ro, nc, Radius, rug, L):
    P3_temp = 0
    P3 = P1
    counter_int = 0
    T_avg = (T2 + T3)/2.
    #iterative loop to compute the pressure drop with error of 1% or until max iteration
    while abs((P3_temp-P3)/P3) > 0.01 and counter_int < max_iterations:
        counter_int += 1
        #protect against P3 < P4
        if P3 < P4:
            P3 = P1
            break
        elif P3_temp != 0:
            P3 = P3_temp
        m_L = valve_LT(P3,P4,ro,Kv,CV,R)

        ro_dP = interp_P_T_DPT(P3,T_avg)[0]
        mu = interp_P_T_mu(P3,T_avg)[0]
        dP = Pressure_drop(m_L/nc,2*Radius, L, mu, ro_dP, rug)
        P3_temp = P1 - dP

    if P3 < P4:
        P3 = P1
    else:
        P3 = P3_temp
    return P3, counter_int

