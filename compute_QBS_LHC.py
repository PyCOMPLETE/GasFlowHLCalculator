from __future__ import division

import numpy as np

import LHCMeasurementTools.TimberManager as tm

from Helium_properties import interp_P_T_hPT, interp_P_T_DPT, interp_P_T_mu
from valve_LT import valve_LT, valve_LT_arr
from Pressure_drop import Pressure_drop
from data_qbs import Data_qbs, data_qbs
from var_getter import VarGetter
import h5_storage

zeros = lambda *x: np.zeros(shape=(x), dtype=float)
max_iterations = 5


class HeatLoadComputer(VarGetter):

    def __init__(self, atd_ob, dq, strict=True, use_dP=True):
        self.use_dP = use_dP
        self.max_iter_cells = []
        self.fail_iter_cells = []

        super(HeatLoadComputer,self).__init__(atd_ob, dq, strict)

        self.computed_values = {}
        self._compute_H()
        self.compute_heat_load()

    def _compute_H(self):
        """
        Computes h3, hC, ro
        """
        P1 = self.data_dict['P1']
        T1 = self.data_dict['T1']
        for key, interpolator in zip(
            ('hC', 'h3', 'ro'),
            (interp_P_T_hPT, interp_P_T_hPT, interp_P_T_DPT)
        ):
            data = zeros(self.Nvalue, self.Ncell)
            for i, isnan in enumerate(self.nan_arr):
                if isnan:
                    data[:,i] = np.nan
                else:
                    for j in xrange(self.Nvalue):
                        data[j,i] = interpolator(P1[j,i], T1[j,i])
            self.computed_values[key] = data

    def _compute_P3(self):
        """
        Iterative computation of P3.
        """
        P3_arr  = zeros(self.Nvalue, self.Ncell)

        P1 = self.data_dict['P1']
        T2 = self.data_dict['T2']
        T3 = self.data_dict['T3']
        P4 = self.data_dict['P4']
        CV = self.data_dict['CV']

        ro = self.computed_values['ro']

        Kv_list = self.dq.Kv_list
        R_list  = self.dq.R_list
        nc_list = self.dq.nc_list
        Radius  = self.dq.Radius
        L_list  = self.dq.L_liist
        rug     = self.dq.rug
        cells   = self.dq.Cell_list

        fails = []
        max_iter = []

        for i, isnan in enumerate(self.nan_arr):
            if isnan:
                P3_arr[:,i] = np.nan
                continue

            Kv = Kv_list[i]
            R  = R_list[i]
            nc = nc_list[i]
            L  = L_list[i]

            for j in xrange(self.Nvalue):
                P3_temp = 0
                counter_int = 0
                P3 = P1[j,i]
                T_avg = (T2[j,i] + T3[j,i])/2.

                #iterative loop to compute the pressure drop with error of 1% or until max iteration
                while abs((P3_temp-P3)/P3) > 0.01 and counter_int < max_iterations:
                    counter_int += 1
                    #protect against P3 < P4
                    if P3 < P4[j,i]:
                        break
                    elif P3_temp != 0:
                        P3 = P3_temp
                    m_L = valve_LT(P3, P4[j,i], ro[j,i], Kv, CV[j,i] ,R)
                    ro_dP = interp_P_T_DPT(P3,T_avg)[0]
                    mu = interp_P_T_mu(P3,T_avg)[0]
                    dP = Pressure_drop(m_L/nc,2*Radius, L, mu, ro_dP, rug)
                    P3_temp = P1 - dP

                if P3 < P4[j,i]:
                    P3_arr[j,i] = P1[j,i]
                    fails.append(cells[i])
                else:
                    P3_arr[j,i] = P3_temp

                if counter_int == max_iterations:
                    max_iter.append(cells[i])

        self.computed_values['P3'] = P3_arr
        self.max_iter_cells = max_iter
        self.fail_iter_cells = fails

    def compute_heat_load(self):

        Qs_list = self.dq.Qs_list
        Kv_list = self.dq.Kv_list
        R_list  = self.dq.R_list

        EH = self.data_dict['EH']
        P4 = self.data_dict['P4']
        CV = self.data_dict['CV']

        h3 = self.computed_values['h3']
        hC = self.computed_values['hC']
        ro = self.computed_values['ro']

        if self.use_dP:
            self._compute_P3()
            P3 = self.computed_values['P3']
        else:
            P3 = self.computed_values['P1']


        m_L = zeros(self.Nvalue, self.Ncell)
        qbs     = zeros(self.Nvalue, self.Ncell)
        for i, isnan in enumerate(self.nan_arr):
            if isnan:
                m_L[:,i] = np.nan
                qbs[:,i] = np.nan
            else:
                m_L[:,i] = valve_LT_arr(P3[:,i],P4[:,i],ro[:,i],Kv_list[i],CV[:,i],R_list[i])
                qbs[:,i] = m_L[:,i]*(h3[:,i]-hC[:,i])-Qs_list[i]-EH[:,i]

        self.computed_values['m_L'] = m_L
        self.computed_values['qbs'] = qbs

    def report(self):
        super(HeatLoadComputer, self).report()

        for list_, title_str in zip(
            (self.fail_iter_cells, self.max_iter_cells),
            ('Iterative approach failed for:\n', 'Maximum iterations for:\n')
        ):
            if list_:
                string = title_str
                for cell in list_:
                    string += '\tcell'
                print(string)

def compute_qbs(atd_ob, use_dP, version=h5_storage.version, strict=True):

    if version != h5_storage.version:
        dq = Data_qbs(version)
    else:
        dq = data_qbs

    Ncell = len(dq.Cell_list) # number of half cells
    Nvalue = len(atd_ob.timestamps) # number of data points in time

    var_getter = HeatLoadComputer(atd_ob, dq, strict)
    data_dict = var_getter.data_dict
    nan_arr = var_getter.nan_arr
    var_getter.report()

    T1          = data_dict['T1']                       #TT961
    T3          = data_dict['T3']                       #TT94x
    CV          = data_dict['CV']                       #CV94x;
    EH          = data_dict['EH']                       #EH84x;
    P1          = data_dict['P1']                       #PT961;
    P4          = data_dict['P4']                       #PT991;
    T2          = data_dict['T2']                       #TT84x

    ro          = zeros(Nvalue, Ncell)                  #density
    m_L         = zeros(Nvalue, Ncell)                  #massflow
    qbs         = zeros(Nvalue, Ncell)                  #BS heat load per half-cell
    hC          = zeros(Nvalue, Ncell)                  #enthalpy of line C
    h3          = zeros(Nvalue, Ncell)                  #enthalpy of BS output

    counter_int = np.zeros((Nvalue, Ncell), dtype=int)  #internal counter
    P3 = np.copy(P1) #intermediate pressure ater the beam screen and before the valve

    rug     = dq.rug
    Radius  = dq.Radius
    R_list  = dq.R_list
    Kv_list = dq.Kv_list
    nc_list = dq.nc_list
    L_list  = dq.L_list
    Qs_list = dq.Qs_list

    for i, isnan in enumerate(nan_arr):
        if isnan:
            continue
        for j in xrange(Nvalue):
            hC[j,i] = interp_P_T_hPT(P1[j,i],T1[j,i])
            h3[j,i] = interp_P_T_hPT(P1[j,i],T3[j,i])
            ro[j,i] = interp_P_T_DPT(P1[j,i],T3[j,i])

            if use_dP:
                #compute the intermediate pressure P3 by iteration
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

