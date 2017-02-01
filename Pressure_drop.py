from __future__ import division
import numpy as np
def Pressure_drop(m, D, L, mu, rho, rug):

    R = D/2;
    S = 3.1416*R**2;
    G = m/S;

    Re = D*G/mu;
    fr = 1./((2.*np.log(D/(2*rug))+1.74))**2

    if not hasattr(Re, '__len__'):
        # For single float values
        if Re < 3e3:
            fl = 64./Re
        elif Re < 1e5:
            fl = 0.3164/Re**0.25
        else:
            fl = 0.221/Re**0.237+0.0032

        ff = max(fl,fr)
    else:
        # For arrays
        fl = np.select([Re < 3e3, Re < 1e5, True], [64./Re, 0.3164/Re**.25, 0.221/Re**0.237+0.0032])
        ff = np.where(fl > fr, fl, fr)

    dP = (ff*m**2 * 1./(rho*S**2)*L/D)/1e5

    return dP

if __name__ == '__main__':
    r = lambda : np.random.random(100)*1e6
    n_args = 6
    function = Pressure_drop
    args = [r() for i in xrange(n_args)]

    single = np.zeros_like(args[0])
    for ctr, arg in enumerate(zip(*args)):
        single[ctr] = function(*arg)
    assert np.all(single == function(*args))

    if False:
        ln_Re = np.linspace(np.log(1e3), np.log(1e6), num=1000)
        Re = np.exp(ln_Re)
        fl = np.select([Re < 3e3, Re < 1e5, True], [64./Re, 0.3164/Re**.25, 0.221/Re**0.237+0.0032])
        import matplotlib.pyplot as plt
        plt.close('all')
        fig = plt.figure()
        fig.set_facecolor('w')
        plt.plot(ln_Re, fl)
        plt.grid(True)
        plt.ylabel('fl', fontsize=18)
        plt.xlabel('log Re', fontsize=18)
        plt.show()

