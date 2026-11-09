"""Grid-refinement stability check of the reported numbers (rerun solve.py logic at 2x density)."""
import numpy as np
from scipy.interpolate import CubicSpline

tau0, gamma0, n = 80.0, 0.002, 3.5
J_target = 0.120
cc = gamma0 / tau0**n
alpha_w = 1.0/(n+1.0); a_exp = n/(n+1.0); a = n/(n+1.0)

def run(nlin, ngeo, tmax, kJ):
    tg = np.concatenate([np.linspace(0.0, 1.0, nlin), np.geomspace(1.0, tmax, ngeo)[1:]])
    th = np.arctan(tg) + np.arctan(n*tg)
    T = a*((1+tg**2)/(1+n**2*tg**2))**(1.0/(2*(n+1)))
    Sr = T*tg/np.sqrt(1+tg**2); St = T/np.sqrt(1+tg**2)
    cp = 1/np.sqrt(1+tg**2); dth = 1/(1+tg**2) + n/(1+n**2*tg**2)
    Ft = np.concatenate([[0.0], np.cumsum(0.5*(T[:-1]**n*cp[:-1]*dth[:-1] + T[1:]**n*cp[1:]*dth[1:])*np.diff(tg))])
    o = np.argsort(th)
    Sr_s, St_s, T_s, F_s = (CubicSpline(th[o], v[o]) for v in (Sr, St, T, Ft))
    thq = np.linspace(-np.pi, np.pi, kJ); tu = np.abs(thq)
    Srv, Stv, Tv = Sr_s(tu), St_s(tu), T_s(tu)
    def JR(R, A):
        tr = A*R**(a_exp-1)*Srv; tt = A*R**(a_exp-1)*Stv*np.sign(thq)
        gm = cc*A**n*R**(alpha_w-1)*Tv**n
        gr = gm*Srv/Tv; gth = gm*Stv/Tv*np.sign(thq)
        W = (n/(n+1.0))*(A*R**(a_exp-1)*Tv)*gm
        return np.trapezoid(W*np.cos(thq) - tr*(gr*np.cos(thq)-gth*np.sin(thq)), thq)*R
    J1 = JR(0.05, 1.0); A = (J_target/J1)**(1/(n+1)); D = cc*A**n
    out = [2*D*0.02**alpha_w*float(F_s(np.pi)), 2*D*0.08**alpha_w*float(F_s(np.pi))]
    t60 = 0.0
    from scipy.optimize import brentq
    t60 = brentq(lambda t: np.arctan(t)+np.arctan(n*t)-np.pi/3, 0.1, 0.4, xtol=1e-15)
    T60 = a*((1+t60**2)/(1+n**2*t60**2))**(1/(2*(n+1)))
    Sr60 = T60*t60/np.sqrt(1+t60**2); St60 = T60/np.sqrt(1+t60**2)
    out += [A*0.02**(a_exp-1)*Sr60, A*0.02**(a_exp-1)*St60]
    for Tth in (a, float(T_s(np.pi/2)), a*n**(-1/(n+1))):
        out.append((D*Tth**n/0.01)**(1/(1-alpha_w)))
    return out

r1 = run(4001, 8001, 1e6, 200001)
r2 = run(16001, 32001, 1e9, 400001)
names = ["dw(0.02) [mm]", "dw(0.08) [mm]", "tau_rz(60deg) [MPa]", "tau_thz(60deg) [MPa]",
         "r(g=0.01, 0deg) [mm]", "r(g=0.01, 90deg) [mm]", "r(g=0.01, 180deg) [mm]"]
print(f"{'quantity':28s} {'coarse':>16s} {'fine':>16s} {'rel.diff':>10s}")
for nm, v1, v2 in zip(names, r1, r2):
    print(f"{nm:28s} {v1:16.9g} {v2:16.9g} {abs(v2-v1)/abs(v2):10.2e}")
