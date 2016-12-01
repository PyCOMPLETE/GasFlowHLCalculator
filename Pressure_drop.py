from __future__ import division
import numpy as np
def Pressure_drop(m, D, L, mu, rho, rug):

    R = D/2;
    S = 3.1416*R**2;
    G = m/S;
      
    Re = D*G/mu;
    fr = 1./((2.*np.log(D/(2*rug))+1.74))**2

    if Re < 3000:
        fl = 64./Re
    elif Re < 100000:
        fl = 0.3164/Re**0.25
    else:
        fl = 0.221/Re**0.237+0.0032

    if fl > fr:
        ff = fl
    else:
        ff = fr

    dP = (ff*m**2 * 1./(rho*S**2)*L/D)/1e5

    return dP
