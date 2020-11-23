
import numpy as np


### Erika Oedlund Virtual instrumentation - Samson method
#gamma = 5./3.
#A = 2.*gamma / (gamma-1.)
A = 9.57
B = 1.2
C = 0.4
###

def K_(x):
    return np.sqrt(A * x**B * (1.-x**C))

def valve_LT(pin, pout, rho, kv, u, R):

    #valve equation from Laurent Tavian
    #u = u-10;              #old formulation with constant pre-constraint.
    u = u-10.*(1-u/100.);   #new formulation with variable pre-constraint.

    if hasattr(pin, '__len__'):
        x = np.where(pout < pin, pout/pin, pin/pout)
        K = np.where(x <= .42, 1., K_(x))
    else:
        if pout < pin:
            x = pout/pin
        else:
            x = pin/pout
        if x <= 0.42:
            K = 1.
        else:
            K = K_(x)

    m_dot = np.sign(pin-pout) * K*1.25e-5*np.sqrt(rho*1e5*pin) * kv * R**(u/100.-1.)
    return m_dot

if __name__ == '__main__':
    r = lambda : np.random.random(1000)
    n_args = 6
    function = valve_LT
    args = [r() for i in range(n_args)]

    single = np.zeros_like(args[0])
    for ctr, arg in enumerate(zip(*args)):
        single[ctr] = function(*arg)
    assert np.all(single == function(*args))

    def K(x):
        K = np.where(x <= .42, 1., np.sqrt(A * x**B * (1.-x**C)))
        return K

    print(A, B, C)
    print(K(np.array([ .41, .42, .43])))
