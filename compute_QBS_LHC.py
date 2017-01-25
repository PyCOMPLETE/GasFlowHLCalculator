from __future__ import division, print_function
import sys
import numpy as np
import LHCMeasurementTools.TimberManager as tm

import h5_storage
from Helium_properties import interp_P_T_hPT, interp_P_T_DPT, interp_P_T_mu
from valve_LT import valve_LT, valve_LT_arr
from Pressure_drop import Pressure_drop
from data_qbs import Data_qbs, data_qbs
from var_getter import VarGetter

zeros = lambda *x: np.zeros(shape=(x), dtype=float)
max_iterations = 5

class HeatLoadComputer(VarGetter):
    __doc__ = VarGetter.__doc__ + '\t-use_dP: Improved, iterative calculation. Default: True'

    def __init__(self, atd_ob, dq, strict=True, report=False, use_dP=True):
        super(HeatLoadComputer,self).__init__(atd_ob, dq, strict, report=False)

        self.use_dP = use_dP

        self.computed_values = {}
        self._compute_H()
        if self.use_dP:
            self._compute_P3()
        self._compute_heat_load()
        if report:
            self.report()

        self.qbs_atd = tm.AlignedTimberData(atd_ob.timestamps, self.computed_values['qbs'], dq.Cell_list)

    def _compute_H(self):
        """
        Computes h3, hC, ro
        """
        P1 = self.data_dict['P1']
        T1 = self.data_dict['T1']
        T3 = self.data_dict['T3']

        hC = zeros(self.Nvalue, self.Ncell)
        h3 = zeros(self.Nvalue, self.Ncell)
        ro = zeros(self.Nvalue, self.Ncell)

        for i, isnan in enumerate(self.nan_arr):
            if isnan:
                hC[:,i] = np.nan
                h3[:,i] = np.nan
                ro[:,i] = np.nan
            else:
                for j in xrange(self.Nvalue):
                    hC[j,i] = interp_P_T_hPT(P1[j,i],T1[j,i])
                    h3[j,i] = interp_P_T_hPT(P1[j,i],T3[j,i])
                    ro[j,i] = interp_P_T_DPT(P1[j,i],T3[j,i])

        self.computed_values['hC'] = hC
        self.computed_values['h3'] = h3
        self.computed_values['ro'] = ro

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
        dq      = self.dq
        Radius  = dq.Radius
        rug     = dq.rug

        for i, isnan in enumerate(self.nan_arr):
            if isnan:
                P3_arr[:,i] = np.nan
                continue

            Kv = dq.Kv_list[i]
            R  = dq.R_list[i]
            nc = dq.nc_list[i]
            L  = dq.L_list[i]

            T_avg = (T2 + T3)/2.

            for j in xrange(self.Nvalue):
                P3_temp = 0
                counter_int = 0
                P3 = P1[j,i]

                #iterative loop to compute the pressure drop with error of 1% or until max iteration
                while abs((P3_temp-P3)/P3) > 0.01 and counter_int < max_iterations:
                    counter_int += 1
                    #protect against P3 < P4
                    if P3 < P4[j,i]:
                        break
                    elif P3_temp != 0:
                        P3 = P3_temp
                    m_L = valve_LT(P3, P4[j,i], ro[j,i], Kv, CV[j,i] ,R)
                    ro_dP = interp_P_T_DPT(P3,T_avg[j,i])[0]
                    mu = interp_P_T_mu(P3,T_avg[j,i])[0]
                    dP = Pressure_drop(m_L/nc,2*Radius, L, mu, ro_dP, rug)
                    P3_temp = P1[j,i] - dP

                if P3 < P4[j,i]:
                    P3_arr[j,i] = P1[j,i]
                    self._insert_to_problem_cells(i, '', 'no_convergence')
                else:
                    P3_arr[j,i] = P3_temp

                if counter_int == max_iterations:
                    self._insert_to_problem_cells(i, j, 'max_iter')

        self.computed_values['P3'] = P3_arr

    def _compute_heat_load(self):

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
            P3 = self.computed_values['P3']
        else:
            P3 = self.data_dict['P1']

        m_L = zeros(self.Nvalue, self.Ncell)
        qbs     = zeros(self.Nvalue, self.Ncell)
        for i, isnan in enumerate(self.nan_arr):
            if isnan:
                m_L[:,i] = np.nan
                qbs[:,i] = np.nan
            else:
                m_L[:,i] = valve_LT_arr(P3[:,i], P4[:,i], ro[:,i], Kv_list[i], CV[:,i], R_list[i])
                qbs[:,i] = m_L[:,i]*(h3[:,i]-hC[:,i])-Qs_list[i]-EH[:,i]

                if np.any(m_L[:,i] < 0):
                    self._insert_to_problem_cells(i, 'm_L', 'negative')

        self.computed_values['m_L'] = m_L
        self.computed_values['qbs'] = qbs

def compute_qbs(atd_ob, use_dP, version=h5_storage.version, strict=True):

    if version != h5_storage.version:
        print((version, h5_storage.version))
        dq = Data_qbs(version)
    else:
        dq = data_qbs

    hl_comp = HeatLoadComputer(atd_ob, dq, strict=strict, use_dP=use_dP)
    hl_comp.report()

    return hl_comp.qbs_atd

