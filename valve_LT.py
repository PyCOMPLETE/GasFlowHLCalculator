from __future__ import division
import numpy as np

def valve_LT(pin,pout,rho,gamma,kv,u,R):
    #valve equation from Laurent Tavian
    #u = u-10;              #old formulation with constant pre-constraint.
    u = u-10.*(1-u/100);   #new formulation with variable pre-constraint.
    try:
        N = len(pin)
    except:
        return valve_LT_single(pin, pout, rho, gamma, kv, u, R)
    K = np.zeros(N);
    x = np.zeros(N);

    for i in xrange(N):
        if pin[i] > pout[i]:
            x[i] = pout[i]/pin[i]
        else:
            x[i] = pin[i]/pout[i]

        if x[i] <= 0.42:
            K[i] = 1
        else:
            K[i] = np.sqrt(9.57 * x[i]**1.2 * (1-x[i]**0.4))

    m = np.sign(pin-pout) * K*1.25e-5*np.sqrt(rho*1e5*pin) * kv/R* np.exp(u/100*np.log(R))
    return m

def valve_LT_single(pin,pout,rho,gamma,kv,u,R):
    #valve equation from Laurent Tavian
    #u = u-10;              #old formulation with constant pre-constraint.
    u = u-10.*(1-u/100);   #new formulation with variable pre-constraint.

    if pin > pout:
        x = pout/pin
    else:
        x = pin/pout

    if x <= 0.42:
        K = 1
    else:
        K = np.sqrt(9.57 * x**1.2 * (1-x**0.4))

    m = np.sign(pin-pout) * K*1.25e-5*np.sqrt(rho*1e5*pin) * kv/R* np.exp(u/100*np.log(R))
    return m


