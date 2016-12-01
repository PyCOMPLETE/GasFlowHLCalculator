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

use_dP = True
show_plot = True
filename = './TIMBER_DATA_Fill5416_LHCBEAMSCREEN_TT84x_injec.csv' #Select the timber file you want to extract

if __name__ != '__main__':
    show_plot = False

timber_data = tm.parse_aligned_csv_file(filename)
timber_data.timestamps -= timber_data.timestamps[0]

Nlist = len(TT94x_list) # number of half cells # 485
Nvalue = timber_data.data.shape[0] # number of data points # 120
Ncell = timber_data.data.shape[1] # 2498
print(Nlist, Nvalue)

T1 = np.zeros(shape=(Nvalue, Nlist)) #TT961
T3 = np.zeros(shape=(Nvalue, Nlist)) #TT94x
CV = np.zeros(shape=(Nvalue, Nlist)) #CV94x;
EH = np.zeros(shape=(Nvalue, Nlist)) #EH84x;
P1 = np.zeros(shape=(Nvalue, Nlist)) #PT961;
P4 = np.zeros(shape=(Nvalue, Nlist)) #PT991;
                                    
if use_dP:
    T2= np.zeros(shape=(Nvalue,Nlist)) ##TT84x

#print(Nvalue, Nlist, len(TT961_list), timber_data.data.shape)

for i in xrange(Nlist):
    j = timber_data.variables.index(TT961_list[i])
    T1[:,i] = timber_data.data[:,j]
    if T1[0,i] == 0 and i > 0:
        T1[:,i] = T1[:,i-1]

    if use_dP:
        try:
            j = timber_data.variables.index(TT84x_list[i])
        except ValueError as e:
            print('Warning: %s' % e)
            T2[:,i] = T2[:,i-1]
        else:
            T2[:,i] = timber_data.data[:,j]

    try:
        j = timber_data.variables.index(TT94x_list[i])
    except ValueError as e:
        print('Warning: %s' % e)
        T3[:,i] = T3[:,i-1]
    else:
        T3[:,i] = timber_data.data[:,j]

    try:
        j = timber_data.variables.index(CV94x_list[i])
    except ValueError as e:
        print('Warning: %s' % e)
        CV[:,i] = CV[:,i-1]
    else:
        CV[:,i] = timber_data.data[:,j]

    try:
        j = timber_data.variables.index(EH84x_list[i])
    except ValueError as e:
        print('Warning: %s' % e)
        EH[:,i] = EH[:,i-1]
    else:
        EH[:,i] = timber_data.data[:,j]

    j = timber_data.variables.index(PT961_list[i])
    P1[:,i] = timber_data.data[:,j]

    j = timber_data.variables.index(PT991_list[i])
    P4[:,i] = timber_data.data[:,j]
    if P4[0,i] == 0 and i > 0:
        P4[:,i] = P4[:,i-1]

if False:
    for ctr, l in enumerate((T1, T2, T3, CV, EH, P4)):
        mask_0 = l == 0
        mask_1 = l != 0
        print(np.sum(mask_0), np.sum(mask_1))

zeros = lambda *x: np.zeros(shape=(x), dtype=float)

ro= zeros(Nvalue,Nlist)            #density
P3_temp= zeros (Nvalue,Nlist)      #temporary intermediate pressure ater the beam screen and before the valve
ro_dP= zeros(Nvalue,Nlist)         #Average density in beam screen
dP= zeros(Nvalue,Nlist)            #Pressure drop in beam screen
gamma= zeros(Nvalue,Nlist)         #ratio of heat capacities
mu= zeros(Nvalue,Nlist)            #viscosity
m_L= zeros(Nvalue,Nlist)           #massflow
Qbs = zeros(Nvalue,Nlist)          #BS heat load per half-cell
QBS_ARC_AVG = zeros(Nvalue,8)      #BS Average values per ARC
hC= zeros(Nvalue,Nlist)            #enthalpy of line C
h3 = zeros(Nvalue,Nlist)           #enthalpy of BS output
counter_int= np.zeros(shape=(Nvalue,Nlist),dtype=int)  #internal counter
 
interp_P_T_hPT = interp2d(P,T,h_PT)
interp_P_T_DPT = interp2d(P,T,D_PT)
interp_P_T_g = interp2d(P,T,gamma_PT)
interp_P_T_mu = interp2d(P,T,mu_PT)

P3 = np.copy(P1) #intermediate pressure ater the beam screen and before the valve
for i in xrange(Nlist):
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
if show_plot:
    plt.close('all')
    t = timber_data.timestamps

    tend = t[-1]/3600.
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(t/3600,QBS_ARC_AVG)
    plt.xlabel('time (hr)')
    plt.ylabel('Qdbs (W)')
    plt.title('Average beam screen heat load per ARC')
    #axis([0 tend 0 250])
    plt.legend(['ARC12','ARC23','ARC34','ARC45','ARC56','ARC67','ARC78','ARC81'])
    plt.subplot(2,1,2)
    plt.plot(t/3600,Qbs)
    plt.xlabel('time (hr)')
    plt.ylabel('Qdbs (W)')
    plt.title('Beam screen heat loads over LHC')
    #axis([0 tend 0 250]);

    plt.show()
