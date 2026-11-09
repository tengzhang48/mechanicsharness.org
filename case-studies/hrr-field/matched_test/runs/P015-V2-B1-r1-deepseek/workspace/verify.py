import numpy as np, sys
sys.path.insert(0,'knowledge/checks')
import solve as S
from crack_tip_checks import check_field, format_report

th = np.linspace(0, np.pi, 4001)
pref = S.K/np.sqrt(2*np.pi)
frr,ftt,frt = S.f_polar(th)
rep = check_field(th, pref*frr, pref*ftt, pref*frt, p=0.5, beta=np.pi,
                  plane="stress", law="linear", nu=S.nu2, mode="I", leading=True)
print(format_report(rep))
