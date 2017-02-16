from Pressure_drop import *
from config_qbs import config_qbs
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', type=str, help='Path where to save the figure', default=None)
parser.add_argument('--noshow', help='Do not make a plot', action='store_true')
arguments = parser.parse_args()

if not arguments.noshow:
    import matplotlib.pyplot as plt
    plt.close('all')

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

    fig = plt.figure()
    fig.set_facecolor('w')

    sp = plt.subplot(1,1,1)
    sp.plot(Re, fl, lw=2, label='f_L')
    sp.grid(True)
    sp.set_ylabel('fl', fontsize=18)
    sp.set_xlabel('Re', fontsize=18)
    sp.set_xscale('log')
    for xx in [3e3, 1e5]:
        sp.axvline(xx, color='black', lw=2)
    for xx in [4.89e3, 4.389e5]:
        sp.axvline(xx, color='red', lw=1)
    fr = calc_fr(config_qbs.Radius, config_qbs.rug)
    sp.axhline(fr, color='black', lw=2, label='f_R')

    sp.legend(loc=1)

    if arguments.o is None:
        plt.show()
    else:
        fig.savefig(arguments.o, dpi=200)
