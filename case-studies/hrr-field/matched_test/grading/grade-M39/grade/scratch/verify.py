import numpy as np, sys, os
sys.path.insert(0, "knowledge/checks")
from crack_tip_checks import check_field, format_report

# effective near-tip material constants (from solve.py)
E=70000.0; nu=0.30; sy=10.0; J=0.8; Czz=-0.005
from scipy.optimize import brentq
P_of = lambda T: 2*nu/E + 1/T - 1/E
Czz_of = lambda T: -P_of(T)*np.sqrt(J*T)/np.sqrt(2*np.pi)
T = brentq(lambda T: Czz_of(T)-Czz, 1, E-1e-6)
D = (1+nu)/E + 1.5*(1/T-1/E); B=(1-2*nu)/E
nu_eff=(D-B)/(2*D+B)
print("T=%.6f nu_eff=%.6f"%(T,nu_eff))

# mode I angular stress functions (r^1/2 * sigma), arbitrary normalization
th = np.linspace(0, np.pi, 4001)
c, s = np.cos(th/2), np.sin(th/2)
S_rr = c*(1+s*s)
S_tt = c*(1-s*s)
S_rt = s*c*c

rep = check_field(th, S_rr, S_tt, S_rt, p=0.5, beta=np.pi, plane="stress",
                  law="linear", nu=nu_eff, mode="I", leading=True, h=1e-5, h2=1e-3)
print(format_report(rep))
