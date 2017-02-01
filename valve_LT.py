from __future__ import division
import numpy as np

gamma = 5./3.

### Erika Oedlund Virtual instrumentation - Samson method
A = 2.*gamma / (gamma-1.)
B = 2./gamma
C = (gamma-1.)/gamma
###

def valve_LT(pin, pout, rho, kv, u, R):

    #valve equation from Laurent Tavian
    #u = u-10;              #old formulation with constant pre-constraint.
    u = u-10.*(1-u/100.);   #new formulation with variable pre-constraint.

    if hasattr(pin, '__len__'):
        x = np.where(pin > pout, pout/pin, pin/pout)
        K = np.where(x <= .42, 1., np.sqrt(A * x**B * (1.-x**C)))
    else:
        if pin > pout:
            x = pout/pin
        else:
            x = pin/pout

        if x <= 0.42:
            K = 1.
        else:
            K = np.sqrt(A * x**B * (1.-x**C))

    m_dot = np.sign(pin-pout) * K*1.25e-5*np.sqrt(rho*1e5*pin) * kv/R* np.exp(u/100.*np.log(R))
    return m_dot

if __name__ == '__main__':
    print A, B, C
    r = lambda : np.random.random(1000)
    n_args = 6
    function = valve_LT
    args = [r() for i in xrange(n_args)]

    single = np.zeros_like(args[0])
    for ctr, arg in enumerate(zip(*args)):
        single[ctr] = function(*arg)
    assert np.all(single == function(*args))
