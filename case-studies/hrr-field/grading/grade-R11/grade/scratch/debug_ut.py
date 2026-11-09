import sys, numpy as np
sys.path.insert(0, "knowledge")
import hrr

n = 6.5
p = 1.0 / (n + 1.0)
s = hrr.singular_exponent(n)

pshoot, sol_ivp = hrr._shoot(True, n)
M = hrr._quadratic_coeff(True, n, s)

def fields(theta):
    z = np.abs(theta)
    Y = sol_ivp.sol(z)
    sgn = np.sign(theta)
    Y = Y * np.vstack([sgn**0, sgn**1, sgn**2, sgn**3])
    f, fp, fpp, fppp = Y
    sr = s * f + fpp
    st = s * (s - 1.0) * f
    srt = (1.0 - s) * fp
    se = np.sqrt(0.75 * (sr - st) ** 2 + 3.0 * srt ** 2)
    return f, fp, fpp, fppp, sr, st, srt, se

thf = np.linspace(0.0, np.pi, 200001)
sef = fields(thf)[7]
c = 1.0 / sef.max()

th = np.linspace(-np.pi, np.pi, 32001)
f, fp, fpp, fppp, sr, st, srt, se = [q * c for q in fields(th)]
rho = se ** (n - 1.0)
ur = rho * 0.75 * (sr - st) / p
utp = rho * 0.75 * (st - sr) - ur
ut = np.concatenate([[0.0], np.cumsum((utp[1:] + utp[:-1]) / 2 * np.diff(th))])

print("max|ur| =", np.abs(ur).max(), " max|utp| =", np.abs(utp).max())
print("max|ut| =", np.abs(ut).max())
print("max|sr-st| =", np.abs(sr - st).max(), " max se =", se.max())
print("se at 0, pi/2, pi:", se[0], se[len(th)//2], se[-1])

# shear route with careful derivative of ur (spectral-ish: high-order FD)
def deriv(y, x):
    return np.gradient(y, x, edge_order=2)
urp = deriv(ur, th)
ut_shear = (urp - 3.0 * rho * srt) / (1.0 - p)
diff = np.abs(ut - ut_shear)
i = np.argmax(diff)
print("max diff", diff[i], "at theta =", th[i], " ut =", ut[i], " ut_shear =", ut_shear[i])
# exclude endpoints
dm = diff[1:-1]; im = np.argmax(dm)
print("max interior diff", dm[im], "at theta =", th[im+1], " ut =", ut[im+1], " ut_shear =", ut_shear[im+1])
# where is ut big?
j = np.argmax(np.abs(ut))
print("ut max at theta =", th[j], "ut =", ut[j], "ut_shear =", ut_shear[j])
# check the eps-theta route against code's ut formula (shear) on interior grid
print("diff stats interior: mean", dm.mean(), " p99", np.percentile(dm, 99))
