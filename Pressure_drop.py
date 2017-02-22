from __future__ import division
import numpy as np

def calc_fl(Re):
    return np.select([Re < 3e3, Re < 1e5, True], [64/Re, 0.3164/Re**0.25, 0.221/Re**0.237+0.0032])

def calc_fr(D, rug):
    return 1./((2.*np.log(D/(2*rug))+1.74))**2

def calc_re(D, m, mu):
    R = D/2.
    S = np.pi*R**2
    G = m/S;
    Re = D*G/mu
    return Re

def pd_factory(D, rug):
    fr = calc_fr(D,rug)
    R = D/2.
    S = np.pi*R**2
    def Pressure_drop(m, L, mu, rho):
        G = m/S;
        Re = D*G/mu
        fl = calc_fl(Re)
        ff = np.where(fl>fr, fl, fr)
        return ff*m**2*L/(D*rho*S**2)*1e-5

    return Pressure_drop

if __name__ == '__main__':
    from config_qbs import config_qbs

    r = lambda : 1e6**np.random.random(100)
    n_args = 6
    def function(*args):
        return pd_factory(*args[:2])(*args[2:])
    args = [r() for i in xrange(n_args)]

    single = np.zeros_like(args[0])
    for ctr, arg in enumerate(zip(*args)):
        single[ctr] = function(*arg)
    assert np.all(single == function(*args))

