"""Follow-up probes: 90deg field identity, opening on T-stress, leading+p<=0, mode II 90, scale=0."""
import sys
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline

sys.path.insert(0, str(Path(__file__).resolve().parent / "crack_tip_asymptotics" / "checks"))
from crack_tip_checks import check_field, williams_exponent, format_report


def failed(rep):
    return {r["check"] for r in rep if r["passed"] is False}


def by_name(rep):
    return {r["check"]: r for r in rep}


# 90-degree field used in tests: srr=4 sin^2, stt=4 cos^2, srt=4 sin cos
# Is this Airy-derived Williams p=0 mode I?
beta = np.pi / 2
th = np.linspace(-beta, beta, 4001)
srr = 4 * np.sin(th) ** 2
stt = 4 * np.cos(th) ** 2
srt = 4 * np.sin(th) * np.cos(th)
print("faces stt,srt at +-beta", stt[0], srt[0], stt[-1], srt[-1])
# Cartesian: sxx = srr c^2 + stt s^2 - 2 srt s c
c, s = np.cos(th), np.sin(th)
sxx = srr * c**2 + stt * s**2 - 2 * srt * s * c
syy = srr * s**2 + stt * c**2 + 2 * srt * s * c
sxy = (srr - stt) * s * c + srt * (c**2 - s**2)
print("sxx min/max", sxx.min(), sxx.max())
print("syy min/max", syy.min(), syy.max())
print("sxy min/max", sxy.min(), sxy.max())
# This is sigma_yy = 4, sigma_xx=sigma_xy=0: uniform tension parallel to the faces? 
# wedge occupies |theta|<pi/2, faces at theta=+-pi/2 which are the y-axis (x=0).
# Traction on face theta=pi/2 (the +y axis? wait theta=pi/2 is +y, outward normal of material |th|<pi/2
# at +pi/2 is +e_theta which is -e_x. Traction-free means sigma_tt=sigma_rt=0.
# Uniform sigma_yy=4 is tension along the bisector, faces are vertical (x=0), so traction on x=0 is sigma_xx, sigma_xy.
# Yes that's the standard non-singular mode I field for a half-plane.

# Williams p=0 mode I Airy at beta=pi/2
p = 0.0
lam = 1.0
sA = lam + 1  # 2
# B = -cos((lam+1)beta)/cos((lam-1)beta) = -cos(pi)/cos(0) = -(-1)/1 = 1
B = -np.cos((lam + 1) * beta) / np.cos((lam - 1) * beta)
print("B", B)
F = np.cos((lam + 1) * th) + B * np.cos((lam - 1) * th)
print("F at 0, +-pi/2", F[th.size//2], F[0], F[-1])
# F = cos(2th) + 1 = 2 cos^2 th
# Sigma_tt = sA*(sA-1)*F = 2*1*F = 2(1+cos 2th)= 4 cos^2 th  MATCHES if amplitude 1
# Sigma_rt = -(sA-1) F' ; F' = -2 sin 2th; so -1 * (-2 sin 2th)= 2 sin 2th = 4 sin th cos th MATCHES
print("test field matches Williams p=0 mode I Airy with F=1+cos2th")

# T-stress opening: U_t ~ 0, passed because < 0 is True for -7e-16
# hydrostatic opening: U_t=0, FAIL because 0 < 0 is False
# Is U_t=0 opening or not? For hydrostatic, faces don't move apart in theta.

# leading=True T-stress: non-singular report AND crack exponent FAIL
# TOPIC says leading=False skips exponent checks; for p<=0 it still applies crack exponent if crack.

# Does 90 wedge with leading=True apply Williams p=0 and pass? Yes.

# Mode II 90 deg: p=-1? That's lambda=2?
print("mode II 90 p", williams_exponent(np.pi/2, "II"))
print("mode II 179deg/2 p", williams_exponent(np.radians(179)/2, "II"))
print("mode II 181deg/2 p", williams_exponent(np.radians(181)/2, "II"))
print("mode I 179", williams_exponent(np.radians(179)/2, "I"))
print("mode I 181", williams_exponent(np.radians(181)/2, "I"))

# table "none singular" for mode II at 180 and 240: p is negative, not a crash.
# At 180, p=-1 is a bounded (actually vanishing as r) field? p=-1 means stress ~ r^{+1}, regular.
# "none singular" is correct if they mean p>0. Returning p=-1 is honest.
# At exactly 257.45, crash.

# scale=0 field
th = np.linspace(-np.pi, np.pi, 401)
zeros = np.zeros_like(th)
try:
    rep = check_field(th, zeros, zeros, zeros, p=0.5, law="linear")
    print("zero field failed", failed(rep))
except Exception as e:
    print("zero field ERROR", type(e).__name__, e)

# opening check on T-stress: ut ~ 0. Should it apply? Faces of a crack with T-stress
# u_theta on faces for T-stress: T-stress is sigma_xx=T, displacement is linear in r, Poisson.
# For plane strain, u_y on y=0 (crack plane) is 0 by symmetry; on faces theta=pi which is the
# negative x-axis, e_theta = -e_y (theta=pi: e_r=-e_x, e_theta=-e_y).
# T-stress: uniaxial sigma_xx. On the crack faces (negative x), the faces are along x, so
# they remain closed in the sense of no opening? T-stress doesn't open the crack.
# The check passed with ut ~ -1e-15. Strictly 0 is FAIL (hydrostatic). Fragile.

# Mixed: leading True T-stress also runs opening and PASSES due to roundoff sign.

# field-scale with tiny field
th, a, b, c = np.linspace(-np.pi, np.pi, 4001), None, None, None
k = 1e-16 / np.sqrt(2 * np.pi)
a = k * (1.25 * np.cos(th / 2) - 0.25 * np.cos(1.5 * th))
b = k * (0.75 * np.cos(th / 2) + 0.25 * np.cos(1.5 * th))
c = k * (0.25 * np.sin(th / 2) + 0.25 * np.sin(1.5 * th))
rep = check_field(th, a, b, c, p=0.5, law="linear")
print("tiny Williams failed", failed(rep))
for name in ["equilibrium (angular ODEs, declared p)",
             "equilibrium (Cartesian, finite differences)",
             "compatibility (Saint-Venant, strains from the declared law)"]:
    print(name, by_name(rep)[name]["value"], "pass", by_name(rep)[name]["passed"])

# Cartesian denom uses max(max(mags), scale) - if mag terms cancel, scale saves it. Good.

# power-law n=1 as linear incompressible
# skip if no hrr

# Does leading=False skip Williams on a 90 wedge? 
th = np.linspace(-np.pi / 2, np.pi / 2, 4001)
rep = check_field(th, 4 * np.sin(th)**2, 4 * np.cos(th)**2, 4 * np.sin(th)*np.cos(th),
                  p=0.0, law="linear", beta=np.pi/2, leading=False)
print("90 leading=False failed", failed(rep))
print([r["check"] for r in rep if "exponent" in r["check"] or "non-singular" in r["check"] or "leading" in r["check"]])

# T-stress as crack with leading True: both non-singular n/a AND crack exponent FAIL
# TOPIC §7 item 7: "for a leading term with p<=0, a report that the field is not singular"
# and "leading=False skips these for higher-order or non-singular terms"
# Inconsistency: if you leave leading=True (default) on T-stress, you get FAIL on crack exponent
# even though they also emit the non-singular report. Documented that leading=False is required.

# Mode II crack opening not checked (n/a). Closing mode II?

# Check test_singularity_thresholds: 179 vs 181
print("I 179", williams_exponent(np.radians(179)/2, "I"))
print("I 181", williams_exponent(np.radians(181)/2, "I"))

# 257.45 labeled as 0 (double root). Computed pII = -2.6e-5 at 257.45 deg.
# table says 0. Honest rounding. At exact crit, finder crashes.

# g(lam) at double root: sin(2 lam beta) - lam sin(2 beta) = 0
# at lam=1: sin(2beta)-sin(2beta)=0 always? Wait sin(2*1*beta)-1*sin(2beta)=0 identically.
# They EXCLUDE lambda=1 for mode II. At the double root, the only remaining nearby root is
# also lambda=1, so after exclusion, no roots -> crash.
# Physics: at the critical angle the double root IS lambda=1, p=0. The code excludes it
# because "lambda=1 is rigid rotation except at isolated angles". The exception is not implemented.
# So table says p=0 (double root) [derived] but williams_exponent(crit, II) raises.
# test_singularity_thresholds uses crit-0.005 and crit+0.005, never the exact point.

from scipy.optimize import brentq
crit = brentq(lambda x: np.tan(x) - x, 4.0, 4.6)
print("pII crit-1e-4", williams_exponent((crit-1e-4)/2, "II"))
print("pII crit+1e-4", williams_exponent((crit+1e-4)/2, "II"))

# Does check_field at exact crit with a p=0 mode II field crash if law=linear notch?
# Construct a p=0 mode II field at that angle? Double root -> need log terms. Skip.

# field-scale: Cartesian uses max(max(mags), scale) but mags are derivatives of stress ~ p/r * stress
# for p=0, mag of derivatives is 0, so denom = scale. Good.

# Saint-Venant uses max(max(mags), estrain) NOT the stress scale. For T-stress estrain is O(1).

# Opening check uses ut_face < 0, not a tolerance. Roundoff can flip.

# hydrostatic on a crack: opening FAIL. Faces don't open (Ut=0). That's arguably correct
# (no opening) but hydrostatic also isn't a crack-face field (not traction-free).

print("done")
