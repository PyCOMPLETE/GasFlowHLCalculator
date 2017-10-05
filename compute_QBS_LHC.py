from __future__ import division, print_function
import numpy as np
import LHCMeasurementTools.TimberManager as tm

import h5_storage
from Helium_properties import interp_P_T_hPT, interp_P_T_DPT, interp_P_T_mu
from valve_LT import valve_LT
from Pressure_drop import pd_factory, calc_re
from config_qbs import Config_qbs, config_qbs

zeros = lambda *x: np.zeros(shape=(x), dtype=float)

class HeatLoadComputer(object):
    """
    This class can be used to relate the raw timber data for a given cell or all cells.
    Parameters:
        -atd_ob: Timber data
        -version: Which config_qbs version to use
        -strict: Raise error if there are missing variables. Default: True
        -details: Print details of report() method.
        -use_dP: Use pressure drop. This takes significantly longer.
        -compute_Re: Compute reynold's number. (Not used any longer).
        -only_raw_data: Only load the pressures, temperatures etc but do not compute anything.
    """

    max_iterations = 5 # For pressure drop

    def __init__(self, atd_ob, version=h5_storage.version, strict=True, details=False, use_dP=True, compute_Re=False, only_raw_data=False):

        # Initialization
        if version == h5_storage.version:
            cq = config_qbs
        else:
            cq = Config_qbs(version)

        self.atd_ob  = atd_ob
        self.cq      = cq
        self.use_dP = use_dP
        self.strict  = strict

        self.Ncell = len(cq.Cell_list)
        self.Nvalue = len(atd_ob.timestamps)

        self.problem_cells = {}
        self.missing_variables = []
        self.nan_arr = np.zeros(len(cq.Cell_list), dtype=np.bool)

        # Read input variables: pressures temperatures, electrical heater etc.
        self._store_all_cell_data()

        # Optionally stop here if no computation is desired
        if only_raw_data:
            self.report(details=details)
            return

        # Main computation
        self.computed_values = {}
        self._compute_ro(use_P3=False)

        if self.use_dP:
            self._compute_P3()
            self._compute_ro(use_P3=True)
        self._compute_H()
        self._compute_heat_load()

        if compute_Re:
            self._compute_Re()

        # Some sanity checks. Then report failing cells to stdout.
        self.assure()
        self.report(details=details)

        # Create object that holds data
        self.qbs_atd = tm.AlignedTimberData(atd_ob.timestamps, self.computed_values['qbs'], cq.Cell_list)

    def _store_all_cell_data(self):
        """
        Create a data_dict that contains all raw data.
        """

        cq = self.cq

        # correct: If no data available, copy it from previous cell.
        # negative: This data may be negative.
        self.var_data_dict = {
            'T1': {
                'vars': cq.TT961_list,
                'correct': True,
                'negative': False,
            },
            'T3': {
                'vars': cq.TT94x_list,
                'correct': False,
                'negative': False,
            },
            'CV': {
                'vars': cq.CV94x_list,
                'correct': False,
                'negative': False,
            },
            'EH': {
                'vars': cq.EH84x_list,
                'correct': False,
                'negative': True,
            },
            'P1': {
                'vars': cq.PT961_list,
                'correct': True,
                'negative': False,
            },
            'P4': {
                'vars': cq.PT991_list,
                'correct': True,
                'negative': False,
            },
            'T2': {
                'vars': cq.TT84x_list,
                'correct': False,
                'negative': False,
            },
        }
        data_dict = {}
        for key, dd in self.var_data_dict.iteritems():
            arr = np.zeros((self.Nvalue, self.Ncell), dtype=float)
            negative_allowed = dd['negative']
            can_be_corrected = dd['correct']

            correct_first = False
            for cell_ctr, var_name in enumerate(dd['vars']):
                try:
                    data = self.atd_ob.dictionary[var_name]
                except KeyError:
                    self.missing_variables.append(var_name)
                    continue

                arr[:,cell_ctr] = data
                if (negative_allowed and np.all(arr[:,cell_ctr] == 0)) or (not negative_allowed and np.all(arr[:,cell_ctr] <= 0)):
                    if can_be_corrected:
                        if cell_ctr != 0:
                            arr[:,cell_ctr] = arr[:,cell_ctr-1]
                            self._insert_to_problem_cells(cell_ctr, var_name, 'corrected')
                        else:
                            correct_first = True
                    else:
                        self._insert_to_problem_cells(cell_ctr, var_name, 'no_data')
                        self.nan_arr[cell_ctr] = True
                elif not negative_allowed and np.any(arr[:,cell_ctr] <= 0):
                    self._insert_to_problem_cells(cell_ctr, var_name, 'negative')

            if correct_first:
                arr[:,0] = arr[:,-1]

            data_dict[key] = arr

            if self.missing_variables:
                print('Warning! Some variables are missing!')
                print(self.missing_variables)
                if self.strict:
                    raise ValueError('Missing variables!')
        self.data_dict = data_dict


    def _compute_H(self):
        """
        Computes h3, hC
        """
        P1 = self.data_dict['P1']
        T1 = self.data_dict['T1']
        T3 = self.data_dict['T3']
        if self.use_dP:
            P3 = self.computed_values['P3']
        else:
            P3 = P1

        hC = zeros(self.Nvalue, self.Ncell)
        h3 = zeros(self.Nvalue, self.Ncell)

        for i, isnan in enumerate(self.nan_arr):
            if isnan:
                hC[:,i] = np.nan
                h3[:,i] = np.nan
            else:
                for j in xrange(self.Nvalue):
                    hC[j,i] = interp_P_T_hPT(P1[j,i], T1[j,i])
                    h3[j,i] = interp_P_T_hPT(P3[j,i], T3[j,i])

        self.computed_values['hC'] = hC
        self.computed_values['h3'] = h3

    def _compute_ro(self, use_P3):
        """
        Computes helium density.
        """
        if use_P3:
            P = self.computed_values['P3']
        else:
            P = self.data_dict['P1']
        T3 = self.data_dict['T3']
        ro = zeros(self.Nvalue, self.Ncell)
        for i, isnan in enumerate(self.nan_arr):
            if isnan:
                ro[:,i] = np.nan
            else:
                for j in xrange(self.Nvalue):
                    ro[j,i] = interp_P_T_DPT(P[j,i],T3[j,i])

        self.computed_values['ro'] = ro

    def _compute_Re(self):
        """
        Computes Reynold's number.
        """
        m_L = self.computed_values['m_L']
        D   = 2*self.cq.Radius
        T_avg = (self.data_dict['T2'] + self.data_dict['T3'])/2.
        if self.use_dP:
            P3 = self.computed_values['P3']
        else:
            P3 = self.data_dict['P1']

        mu = np.zeros_like(m_L)
        for i in xrange(self.Ncell):
            for j in xrange(self.Nvalue):
                mu[j,i] = interp_P_T_mu(P3[j,i], T_avg[j,i])

        self.computed_values['Re'] = calc_re(D, m_L, mu)

    def _compute_P3(self):
        """
        Iterative computation of P3 (pressure drop).
        """
        P3_arr  = zeros(self.Nvalue, self.Ncell)

        P1 = self.data_dict['P1']
        T2 = self.data_dict['T2']
        T3 = self.data_dict['T3']
        P4 = self.data_dict['P4']
        CV = self.data_dict['CV']
        ro = self.computed_values['ro']
        cq      = self.cq
        Pressure_drop = pd_factory(2*cq.Radius, cq.rug)

        for i, isnan in enumerate(self.nan_arr):
            if isnan:
                P3_arr[:,i] = np.nan
                continue

            Kv = cq.Kv_list[i]
            R  = cq.R_list[i]
            nc = cq.nc_list[i]
            L  = cq.L_list[i]

            T_avg = (T2 + T3)/2.

            for j in xrange(self.Nvalue):
                P3_temp = 0
                iterations = 0
                P3 = P1[j,i]

                #iterative loop to compute the pressure drop with error of 1% or until max iteration
                while abs((P3_temp-P3)/P3) > 0.01 and iterations < self.max_iterations:
                    iterations += 1
                    #protect against P3 < P4
                    if P3 < P4[j,i]:
                        break
                    elif P3_temp != 0:
                        P3 = P3_temp

                    m_L         = valve_LT(P3, P4[j,i], ro[j,i], Kv, CV[j,i] ,R)
                    ro_dP       = interp_P_T_DPT(P3,T_avg[j,i])[0]
                    mu          = interp_P_T_mu(P3,T_avg[j,i])[0]
                    dP          = Pressure_drop(m_L/nc, L, mu, ro_dP)
                    P3_temp     = P1[j,i] - dP

                if P3 < P4[j,i]:
                    P3_arr[j,i] = P1[j,i]
                    self._insert_to_problem_cells(i, j, 'no_convergence')
                else:
                    P3_arr[j,i] = P3_temp

                if iterations == self.max_iterations:
                    self._insert_to_problem_cells(i, j, 'max_iter')

        self.computed_values['P3'] = P3_arr

    def _compute_heat_load(self):
        """
        Final step: mass flow and heat load.
        """

        Qs_list = self.cq.Qs_list
        Kv_list = self.cq.Kv_list
        R_list  = self.cq.R_list

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
        qbs = zeros(self.Nvalue, self.Ncell)
        for i, isnan in enumerate(self.nan_arr):
            if isnan:
                m_L[:,i] = np.nan
                qbs[:,i] = np.nan
            else:
                m_L[:,i] = valve_LT(P3[:,i], P4[:,i], ro[:,i], Kv_list[i], CV[:,i], R_list[i])
                qbs[:,i] = m_L[:,i]*(h3[:,i]-hC[:,i])-Qs_list[i]-EH[:,i]

        self.computed_values['m_L'] = m_L
        self.computed_values['qbs'] = qbs

    def get_single_cell_data(self, cell):
        """
        Returns recomputed and raw data for one single cell.
        """
        output_dict = {}
        candidate_cells = filter(lambda x: cell in x, list(self.cq.Cell_list))
        if len(candidate_cells) == 0:
            raise ValueError('Cell %s not found!' % cell)
        elif len(candidate_cells) != 1:
            raise ValueError('Too many cells with %s: %s' % (cell, candidate_cells))
        else:
            index = list(self.cq.Cell_list).index(candidate_cells[0])

        for key, arr in self.computed_values.iteritems():
            output_dict[key] = arr[:,index]
        for key, arr in self.data_dict.iteritems():
            output_dict[key] = arr[:,index]

        return output_dict

    def assure(self):
        """
        Check conditions:
        P1 > P4
        m_L > 0
        """

        P1 = self.data_dict['P1']
        P4 = self.data_dict['P4']
        for cell_ctr, isnan in enumerate(self.nan_arr):
            if not isnan and np.any(P1[:,cell_ctr] < P4[:,cell_ctr]):
                self._insert_to_problem_cells(cell_ctr, 'P1 < P4', 'failed_checks')

        m_L = self.computed_values['m_L']
        for cell_ctr, isnan in enumerate(self.nan_arr):
            if not isnan:
                if np.any(m_L[:,cell_ctr] < 0):
                    self._insert_to_problem_cells(cell_ctr, 'm_L', 'negative')

    def report(self, details=False):
        """
        Print out the missing and corrected variables as well as the nan cells.
        """

        for type_, dd in self.problem_cells.iteritems():
            print('%i problems of type %s' % (len(dd), type_))
            if details:
                for cell, subdict in dd.iteritems():
                    print('%s in S%s of type %s: %s' % (cell, subdict['sector'], subdict['type'], subdict['list']))

    def _insert_to_problem_cells(self, cell_ctr, var, type_):
        """
        Utility to store what kind of problems and for which cells they occur.
        """

        problem_cells = self.problem_cells
        cell = self.cq.Cell_list[cell_ctr]
        if type_ not in problem_cells:
            problem_cells[type_] = {}
        if cell not in problem_cells[type_]:
            problem_cells[type_][cell] = {
                'sector': self.cq.Sector_list[cell_ctr],
                'type': self.cq.Type_list[cell_ctr],
                'list': set(),
            }
        problem_cells[type_][cell]['list'].add(var)

# Main interface of this file
def compute_qbs(atd_ob, use_dP, version=h5_storage.version, strict=True, details=False):
    hl_comp = HeatLoadComputer(atd_ob, version=version, strict=strict, use_dP=use_dP, details=details)
    return hl_comp.qbs_atd

