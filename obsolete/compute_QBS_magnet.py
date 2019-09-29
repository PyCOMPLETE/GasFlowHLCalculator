from __future__ import division
import numpy as np

zeros = lambda *x: np.zeros(x, dtype=float)

class InvalidDataError(ValueError):
    pass

class QbsMagnetCalculator(object):
    def __init__(self, interp_P_T_hPT, atd, P1, m_L, cell_list):
        self.P1 = P1
        self.m_L = m_L
        self.variables = list(atd.variables)
        self.n_tt = len(atd.timestamps)
        self.data = atd.data
        self.interp_P_T_hPT = interp_P_T_hPT
        self.variables_set = set(self.variables)
        self.cell_list = cell_list
        self.data_dictionary = atd.dictionary

    def Compute_QBS_magnet(self, cell, Tin_list, Tout_list):
        print('Warning! Obsolete function Compute_QBS_magnet is called. Call Compute_QBS_magnet_single instead!')
        cell_index = self.cell_list.index(cell)
        if not isinstance(Tin_list, list):
            Tin_list = [Tin_list]
        if not isinstance(Tout_list, list):
            Tout_list = [Tout_list]

        Tin_list = list(set(Tin_list).intersection(self.variables_set))
        Tout_list = list(set(Tout_list).intersection(self.variables_set))

        n_tt = self.n_tt
        n_in = len(Tin_list)
        n_out = len(Tout_list)

        Tin= zeros(n_tt,n_in)
        Tout= zeros(n_tt,n_out)
        hin= zeros(n_tt)
        hout = zeros(n_tt)

        variables = self.variables
        data = self.data
        P1 = self.P1
        m_L = self.m_L
        interp_P_T_hPT = self.interp_P_T_hPT

        #build vector of data
        for tt_var_list, tt_data_arr in [(Tin_list, Tin),(Tout_list, Tout)]:
            for i, tt_var in enumerate(tt_var_list):
                index = variables.index(tt_var)
                if np.all(data[:,index] == 0):
                    tt_data_arr[:,i] = np.nan
                else:
                    tt_data_arr[:,i] = data[:,index]


        #compute average if 2 sensors
        if n_in > 1:
            Tin_avg = np.nanmean(Tin,axis=1)
        else:
            Tin_avg = Tin

        if n_out > 1:
            Tout_avg = np.nanmean(Tout, axis=1)
        else:
            Tout_avg = Tout

        #Beam screen load calculation
        for i in xrange(n_tt):
            hin[i] = interp_P_T_hPT(P1[i,cell_index],Tin_avg[i])
            try:
                hout[i] = interp_P_T_hPT(P1[i,cell_index],Tout_avg[i])
            except:
                import pdb ; pdb.set_trace()
        QBS_mag = m_L[:,cell_index]*(hout-hin)

        return QBS_mag

    def Compute_QBS_magnet_single(self, cell, Tin_var, Tout_var):
        cell_index = self.cell_list.index(cell)
        P1 = self.P1
        n_tt = self.n_tt
        data_dictionary = self.data_dictionary

        hin= zeros(n_tt)
        hout = zeros(n_tt)

        if Tin_var in data_dictionary and Tout_var in data_dictionary:
            Tin = data_dictionary[Tin_var]
            Tout = data_dictionary[Tout_var]
        else:
            raise InvalidDataError('Key not found')

        if np.all(Tin == 0) or np.all(Tout == 0):
            raise InvalidDataError('Data is 0')

        for i in xrange(n_tt):
            hin[i] = self.interp_P_T_hPT(P1[i,cell_index],Tin[i])
            hout[i] = self.interp_P_T_hPT(P1[i,cell_index],Tout[i])
        QBS_mag = self.m_L[:,cell_index]*(hout-hin)

        return QBS_mag

