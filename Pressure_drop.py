from __future__ import division
import numpy as np

take_statistics = True

if take_statistics:
    fl_arr = []
    fr_val = []

def calc_fl(Re):
    return np.select([Re < 3e3, Re < 1e5, True], [64/Re, 0.3164/Re**0.25, 0.221/Re**0.237+0.0032])

def pd_factory(D, rug):
    fr = 1./((2.*np.log(D/(2*rug))+1.74))**2
    if take_statistics:
        fr_val.append(fr)
    R = D/2;
    S = 3.1416*R**2;
    def Pressure_drop(m, L, mu, rho):
        if take_statistics:
            G = m/S;
            Re = D*G/mu
            fl_arr.append(calc_fl(Re))

        return fr*m**2*L/(D*rho*S**2)*1e-5

    return Pressure_drop

if __name__ == '__main__':
    r = lambda : np.random.random(100)*1e6
    n_args = 6
    def function(*args):
        return pd_factory(*args[:2])(*args[2:])
    args = [r() for i in xrange(n_args)]

    single = np.zeros_like(args[0])
    for ctr, arg in enumerate(zip(*args)):
        single[ctr] = function(*arg)
    assert np.all(single == function(*args))

    if True:
        ln_Re = np.linspace(np.log(1e3), np.log(1e6), num=1000)
        Re = np.exp(ln_Re)
        fl = np.select([Re < 3e3, Re < 1e5, True], [64./Re +.021419, 0.3164/Re**.25, 0.221/Re**0.237+0.0032])
        import matplotlib.pyplot as plt
        plt.close('all')
        fig = plt.figure()
        fig.set_facecolor('w')
        plt.plot(ln_Re, fl)
        plt.grid(True)
        plt.ylabel('fl', fontsize=18)
        plt.xlabel('log Re', fontsize=18)
        plt.show()

