from __future__ import division
import numpy as np

zeros = lambda *x: np.zeros((x), dtype=float)

class QbsMagnetCalculator(object):
    def __init__(self, interp_P_T_hPT, atd, P1, m_L):
        self.P1 = P1
        self.m_L = m_L
        self.variables = list(atd.variables)
        self.n_tt = len(atd.timestamps)
        self.data = atd.data
        self.interp_P_T_hPT = interp_P_T_hPT

    def Compute_QBS_magnet(self, n, Tin_list, Tout_list):
        if not isinstance(Tin_list, list):
            Tin_list = [Tin_list]
        if not isinstance(Tout_list, list):
            Tout_list = [Tout_list]

        n_tt = self.n_tt
        n_in = len(Tin_list)
        n_out = len(Tout_list)

        Tin= zeros(n_tt,n_in)
        Tout= zeros(n_tt,n_out)
        Tin_avg= zeros(n_tt)
        Tout_avg= zeros(n_tt)
        hin= zeros(n_tt)
        hout = zeros(n_tt)
        QBS_mag= zeros(n_tt)

        variables = self.variables
        data = self.data
        P1 = self.P1
        m_L = self.m_L
        interp_P_T_hPT = self.interp_P_T_hPT

        #build vector of data
        for i in xrange(n_in):
            index = variables.index(Tin_list[i])
            Tin[:,i] = data[:,index]

        for i in xrange(n_out):
            index = variables.index(Tout_list[i])
            Tout[:,i] = data[:,index]

        #compute average if 2 sensors
        if n_in > 1:
            Tin_avg = np.mean(Tin,axis=1)
        else:
            Tin_avg = Tin

        if n_out > 1:
            Tout_avg = np.mean(Tout, axis=1)
        else:
            Tout_avg = Tout

        #Beam screen load calculation
        for i in xrange(n_tt):
            hin[i] = interp_P_T_hPT(P1[i,n],Tin_avg[i])
            hout[i] = interp_P_T_hPT(P1[i,n],Tout_avg[i])
        QBS_mag = m_L[:,n]*(hout-hin)

        return QBS_mag
