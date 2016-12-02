from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d

sys.path.append('..')
from LHCMeasurementTools import TimberManager as tm

from Helium_properties import *
from data_QBS_LHC import *
from valve_LT import valve_LT
from Pressure_drop import Pressure_drop

def compute_QBS_LHC(aligned_timber_data, use_dP, return_qbs=False):

    atd_ob = aligned_timber_data
    
    Ncell = len(TT94x_list) # number of half cells
    Nvalue = atd_ob.data.shape[0] # number of data points
    Nlist = atd_ob.data.shape[1] # Number of variables

    T1 = np.zeros(shape=(Nvalue, Ncell)) #TT961
    T3 = np.zeros(shape=(Nvalue, Ncell)) #TT94x
    CV = np.zeros(shape=(Nvalue, Ncell)) #CV94x;
    EH = np.zeros(shape=(Nvalue, Ncell)) #EH84x;
    P1 = np.zeros(shape=(Nvalue, Ncell)) #PT961;
    P4 = np.zeros(shape=(Nvalue, Ncell)) #PT991;
                                        
    if use_dP:
        T2= np.zeros(shape=(Nvalue,Ncell)) ##TT84x

    arr_list = [T1, T3, CV, EH, P1, P4]
    name_list = [TT961_list, TT94x_list, CV94x_list, EH84x_list, PT961_list, PT991_list]
    if use_dP:
        arr_list.append(T2)
        name_list.append(TT84x_list)

    for i in xrange(Ncell):
        for arr, names in zip(arr_list, name_list):
            try:
                j = atd_ob.variables.index(names[i])
            except ValueError as e:
                print('Warning: %s' % e)
                arr[:,i] = arr[:,i-1]
            else:
                arr[:,i] = atd_ob.data[:,j]
                if arr[0,i] == 0 and i > 0:
                    arr[:,i] = arr[:,i-1]

    zeros = lambda *x: np.zeros(shape=(x), dtype=float)

    ro= zeros(Nvalue,Ncell)            #density
    P3_temp= zeros (Nvalue,Ncell)      #temporary intermediate pressure ater the beam screen and before the valve
    ro_dP= zeros(Nvalue,Ncell)         #Average density in beam screen
    dP= zeros(Nvalue,Ncell)            #Pressure drop in beam screen
    gamma= zeros(Nvalue,Ncell)         #ratio of heat capacities
    mu= zeros(Nvalue,Ncell)            #viscosity
    m_L= zeros(Nvalue,Ncell)           #massflow
    Qbs = zeros(Nvalue,Ncell)          #BS heat load per half-cell
    QBS_ARC_AVG = zeros(Nvalue,8)      #BS Average values per ARC
    hC= zeros(Nvalue,Ncell)            #enthalpy of line C
    h3 = zeros(Nvalue,Ncell)           #enthalpy of BS output
    counter_int= np.zeros(shape=(Nvalue,Ncell),dtype=int)  #internal counter
     
    interp_P_T_hPT = interp2d(P,T,h_PT)
    interp_P_T_DPT = interp2d(P,T,D_PT)
    interp_P_T_g = interp2d(P,T,gamma_PT)
    interp_P_T_mu = interp2d(P,T,mu_PT)

    P3 = np.copy(P1) #intermediate pressure ater the beam screen and before the valve
    for i in xrange(Ncell):
        hC[:,i] = np.diag(interp_P_T_hPT(P1[:,i],T1[:,i]))
        h3[:,i] = np.diag(interp_P_T_hPT(P1[:,i],T3[:,i]))
        ro[:,i] = np.diag(interp_P_T_DPT(P1[:,i],T3[:,i]))
        gamma[:,i] = np.diag(interp_P_T_g(P1[:,i],T3[:,i]))

        #compute the intermediate pressure P3 by iteration
        if use_dP:
            for j in xrange(Nvalue):    #loop along the data values
                while abs((P3_temp[j,i]-P3[j,i])/P3[j,i]) > 0.01 and counter_int[j,i] < 5:
                    #iterative loop to compute the pressure drop with error of 1% or 5 max iteration
                    counter_int[j,i] += 1
                    #protect against P3 < P4
                    if P3[j,i] < P4[j,i]:
                        P3[j,i] = P1[j,i]
                        break
                    elif P3_temp[j,i] != 0:
                        P3[j,i] = P3_temp[j,i]
                    m_L[j,i] = valve_LT(P3[j,i],P4[j,i],ro[j,i],gamma[j,i],Kv_list[i],CV[j,i],R_list[i])
                    ro_dP[j,i] = interp_P_T_DPT(P3[j,i],(T2[j,i]+T3[j,i])/2)
                    mu[j,i] = interp_P_T_mu(P3[j,i],(T2[j,i]+T3[j,i])/2)
                    dP[j,i] = Pressure_drop(m_L[j,i]/nc_list[i],2*Radius, L_list[i], mu[j,i], ro_dP[j,i], rug)
                    P3_temp[j,i] = P1[j,i] - dP[j,i];
            
                if P3[j,i] < P4[j,i]:
                    P3[j,i] = P1[j,i]
                else:
                    P3[j,i] = P3_temp[j,i]

        m_L[:,i] = valve_LT(P3[:,i],P4[:,i],ro[:,i],gamma[:,i],Kv_list[i],CV[:,i],R_list[i])
        Qbs[:,i] = m_L[:,i]*(h3[:,i]-hC[:,i])-Qs_list[i]-EH[:,i]

    #remove NaN values
    mask_nan = np.isnan(Qbs)
    mask_not_nan = np.logical_not(mask_nan)
    print('There are %i nan out of %i values in Qbs!' % (np.sum(mask_nan), np.sum(mask_not_nan)))
    Qbs[mask_nan] = 0.

    #compute average per ARC
    for k in xrange(8):
        first = arc_index[k,0]
        last = arc_index[k,1]
        for i in xrange(Nvalue):
            QBS_ARC_AVG[i,k] = np.mean(Qbs[i,first:last])
    arc_list = ['ARC12','ARC23','ARC34','ARC45','ARC56','ARC67','ARC78','ARC81']

    if return_qbs:
        return QBS_ARC_AVG, arc_list, Qbs
    else:
        return QBS_ARC_AVG, arc_list

if __name__ == '__main__':
    show_plot = True
    use_dP = True
    filename = './TIMBER_DATA_Fill5416_LHCBEAMSCREEN_TT84x_injec.csv' #Select the timber file you want to extract
    atd_ob = tm.parse_aligned_csv_file(filename)
    atd_ob.timestamps -= atd_ob.timestamps[0]

    QBS_ARC_AVG, arc_list, Qbs = compute_QBS_LHC(atd_ob, use_dP, return_qbs=True)

    if show_plot:
        plt.close('all')
        t = atd_ob.timestamps

        tend = t[-1]/3600.
        plt.figure()
        plt.subplot(2,1,1)
        plt.plot(t/3600,QBS_ARC_AVG)
        plt.xlabel('time (hr)')
        plt.ylabel('Qdbs (W)')
        plt.title('Average beam screen heat load per ARC')
        #axis([0 tend 0 250])
        plt.legend(arc_list)
        plt.subplot(2,1,2)
        plt.plot(t/3600,Qbs)
        plt.xlabel('time (hr)')
        plt.ylabel('Qdbs (W)')
        plt.title('Beam screen heat loads over LHC')
        #axis([0 tend 0 250]);

        plt.show()
