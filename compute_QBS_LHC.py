from __future__ import division
import sys
import os
import h5py

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d

if '..' not in sys.path: sys.path.append('..')
import LHCMeasurementTools.TimberManager as tm

from Helium_properties import P, T, D_PT, h_PT, mu_PT
from valve_LT import valve_LT
from Pressure_drop import Pressure_drop
from data_qbs import Data_qbs, data_qbs
import h5_storage

zeros = lambda *x: np.zeros(shape=(x), dtype=float)
dq_latest = data_qbs

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
                                        
    if use_dP:
        T2= np.zeros(shape=(Nvalue,Ncell)) ##TT84x

    arr_list = [T1, T3, CV, EH, P1, P4]
    name_list = [dq.TT961_list, dq.TT94x_list, dq.CV94x_list, dq.EH84x_list, dq.PT961_list, dq.PT991_list]
    correct_zeros = [True, False, False, False, True, True]
    if use_dP:
        arr_list.append(T2)
        name_list.append(dq.TT84x_list)
        correct_zeros.append(False)

    variables = list(atd_ob.variables)
    missing_variables = []
    corrected_variables = []

    for ctr, (arr, names, correct_zero) in enumerate(zip(arr_list, name_list, correct_zeros)):
        for i in xrange(Ncell):
            try:
                j = variables.index(names[i])
            except ValueError as e:
                missing_variables.append(names[i])
                arr[:,i] = arr[:,i-1]
            else:
                arr[:,i] = atd_ob.data[:,j]
                if correct_zero and arr[0,i] == 0 and i > 0:
                    arr[:,i] = arr[:,i-1]
                    corrected_variables.append(names[i])
    if missing_variables:
        print('Missing variables:', missing_variables)
        if strict:
            raise ValueError('There have been missing variables!')
    if corrected_variables:
        print('Corrected variables:', corrected_variables)

    ro          = zeros(Nvalue,Ncell)    #density
    P3_temp     = zeros(Nvalue,Ncell)    #temporary intermediate pressure ater the beam screen and before the valve
    ro_dP       = zeros(Nvalue,Ncell)    #Average density in beam screen
    dP          = zeros(Nvalue,Ncell)    #Pressure drop in beam screen
    gamma       = zeros(Nvalue,Ncell)    #ratio of heat capacities
    mu          = zeros(Nvalue,Ncell)    #viscosity
    m_L         = zeros(Nvalue,Ncell)    #massflow
    qbs         = zeros(Nvalue,Ncell)    #BS heat load per half-cell
    hC          = zeros(Nvalue,Ncell)    #enthalpy of line C
    h3          = zeros(Nvalue,Ncell)    #enthalpy of BS output
    counter_int = np.zeros(shape=(Nvalue,Ncell),dtype=int)  #internal counter
     
    interp_P_T_hPT = interp2d(P,T,h_PT)
    interp_P_T_DPT = interp2d(P,T,D_PT)
    interp_P_T_mu = interp2d(P,T,mu_PT)
    #interp_P_T_g = interp2d(P,T,gamma_PT) # not used

    P3 = np.copy(P1) #intermediate pressure ater the beam screen and before the valve

    rug = dq.rug
    Radius = dq.Radius
    R_list = dq.R_list
    Kv_list = dq.Kv_list
    nc_list = dq.nc_list
    L_list = dq.L_list
    Qs_list = dq.Qs_list

    max_iterations = 5
    for i in xrange(Ncell):
        for j in xrange(Nvalue):
            hC[j,i] = interp_P_T_hPT(P1[j,i],T1[j,i])
            h3[j,i] = interp_P_T_hPT(P1[j,i],T3[j,i])
            ro[j,i] = interp_P_T_DPT(P1[j,i],T3[j,i])
            #gamma[j,i] = interp_P_T_g(P1[j,i],T3[j,i]) # not used

        #compute the intermediate pressure P3 by iteration
        if use_dP:
            for j in xrange(Nvalue):    #loop along the data values
                while abs((P3_temp[j,i]-P3[j,i])/P3[j,i]) > 0.01 and counter_int[j,i] < max_iterations:
                    #iterative loop to compute the pressure drop with error of 1% or 5 max iteration
                    counter_int[j,i] += 1
                    #protect against P3 < P4
                    if P3[j,i] < P4[j,i]:
                        P3[j,i] = P1[j,i]
                        break
                    elif P3_temp[j,i] != 0:
                        P3[j,i] = P3_temp[j,i]
                    m_L[j,i] = valve_LT(P3[j,i],P4[j,i],ro[j,i],gamma[j,i],Kv_list[i],CV[j,i],R_list[i])
                    ro_dP[j,i] = interp_P_T_DPT(P3[j,i],(T2[j,i]+T3[j,i])/2.)
                    mu[j,i] = interp_P_T_mu(P3[j,i],(T2[j,i]+T3[j,i])/2.)
                    dP[j,i] = Pressure_drop(m_L[j,i]/nc_list[i],2*Radius, L_list[i], mu[j,i], ro_dP[j,i], rug)
                    P3_temp[j,i] = P1[j,i] - dP[j,i];
            
                if P3[j,i] < P4[j,i]:
                    P3[j,i] = P1[j,i]
                else:
                    P3[j,i] = P3_temp[j,i]

        m_L[:,i] = valve_LT(P3[:,i],P4[:,i],ro[:,i],gamma[:,i],Kv_list[i],CV[:,i],R_list[i])
        qbs[:,i] = m_L[:,i]*(h3[:,i]-hC[:,i])-Qs_list[i]-EH[:,i]

    # Protect against P1 == 0
    mask_0 = P1 == 0
    qbs[mask_0] = 0

    n_max_iterations = np.sum(counter_int == max_iterations)
    if n_max_iterations != 0:
        print('Warning: Maximum iterations were reached %i times!' % n_max_iterations)

    return tm.AlignedTimberData(atd_ob.timestamps, qbs, dq.Cell_list)

