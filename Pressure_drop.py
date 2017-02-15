from __future__ import division
import numpy as np

def calc_fl(Re):
    return np.select([Re < 3e3, Re < 1e5, True], [64/Re, 0.3164/Re**0.25, 0.221/Re**0.237+0.0032])

def calc_fr(D, rug):
    return 1./((2.*np.log(D/(2*rug))+1.74))**2

def pd_factory(D, rug):
    fr = calc_fr(D,rug)
    R = D/2;
    S = np.pi*R**2
    def Pressure_drop(m, L, mu, rho):
        G = m/S;
        Re = D*G/mu
        fl = calc_fl(Re)
        ff = np.where(fl>fr, fl, fr)
        return ff*m**2*L/(D*rho*S**2)*1e-5

    return Pressure_drop

if __name__ == '__main__':
    import argparse
    import matplotlib.pyplot as plt
    from config_qbs import config_qbs

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', type=str, help='Path where to save the figure', default=None)
    arguments = parser.parse_args()

    r = lambda : np.random.random(100)*1e6
    n_args = 6
    def function(*args):
        return pd_factory(*args[:2])(*args[2:])
    args = [r() for i in xrange(n_args)]

    single = np.zeros_like(args[0])
    for ctr, arg in enumerate(zip(*args)):
        single[ctr] = function(*arg)
    assert np.all(single == function(*args))

    ln_Re = np.linspace(np.log(1e3), np.log(1e7), num=1000)
    Re = np.exp(ln_Re)
    fl = calc_fl(Re)

    try:
        from RcParams import init_pyplot
        init_pyplot()
    except ImportError:
        ticksize = 18
        plt.rcParams['font.size'] = ticksize
        plt.rcParams['ytick.labelsize'] = ticksize
        plt.rcParams['xtick.labelsize'] = ticksize

    plt.close('all')

    fig = plt.figure()
    fig.set_facecolor('w')

    sp = plt.subplot(1,1,1)
    sp.plot(Re, fl, lw=2, label='f_L')
    #plt.title("Friction factor depending on Reynold's number")
    sp.grid(True)
    sp.set_ylabel('fl', fontsize=18)
    sp.set_xlabel('Re', fontsize=18)
    sp.set_xscale('log')
    for xx in [3e3, 1e5]:
        sp.axvline(xx, color='black', lw=1)
    for xx in [4.89e3, 4.389e5]:
        sp.axvline(xx, color='red', lw=1)
    fr = calc_fr(config_qbs.Radius, config_qbs.rug)
    sp.axhline(fr, color='green', lw=2, label='f_R')

    sp.legend(loc=1)

    if arguments.o is None:
        plt.show()
    else:
        fig.savefig(arguments.o, dpi=200)

