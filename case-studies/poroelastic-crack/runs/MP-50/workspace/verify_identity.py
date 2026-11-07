"""Verify the reciprocity identity used to get the crack volume.

Claim: for the unit-pressure crack (face opening traction tau=1) in the
drained elastic solid,  Int (sigma_kk)^2 dA = 2*pi*L^2  (per unit thickness).

sigma_kk is the trace of the crack stress field:
    sigma_kk = 2*(Re[ z/sqrt(z^2-L^2) ] - 1)
Elliptic coords: z = L cosh(mu + i theta), mu in (0,inf), theta in (0,2pi).
    Re[ z/sqrt(z^2-L^2) ] = cosh(mu)sinh(mu) / D,
    D = sinh(mu)^2 cos(theta)^2 + cosh(mu)^2 sin(theta)^2
    dA = L^2 (sinh(mu)^2 + sin(theta)^2) dmu dtheta
"""
import numpy as np

L = 1.0

def sigma_kk(mu, th):
    D = np.sinh(mu)**2 * np.cos(th)**2 + np.cosh(mu)**2 * np.sin(th)**2
    return 2.0 * (np.cosh(mu) * np.sinh(mu) / D - 1.0)

def integrate(nmu, nth, mumax):
    # Gauss-Legendre in mu (transform to kill the 1/sqrt-ish tip behaviour:
    # use mu = t^2 to cluster points near mu=0)
    t, wt = np.polynomial.legendre.leggauss(nmu)
    t = 0.5 * (t + 1.0)          # t in (0,1)
    mu = mumax * t**2
    dmu = mumax * 2.0 * t * (0.5)
    wmu = wt * dmu
    q, wq = np.polynomial.legendre.leggauss(nth)
    th = 0.5 * (q + 1.0) * 2.0 * np.pi
    wth = wq * np.pi
    MU, TH = np.meshgrid(mu, th, indexing="ij")
    W = np.outer(wmu, wth)
    S = sigma_kk(MU, TH)
    jac = L**2 * (np.sinh(MU)**2 + np.sin(TH)**2)
    return np.sum(S**2 * jac * W)

for nmu, nth, mumax in [(400,400,20), (800,800,30), (1600,1600,40), (3000,3000,50)]:
    val = integrate(nmu, nth, mumax)
    print(f"nmu={nmu:5d} nth={nth:5d} mumax={mumax:3d}  Int sigma_kk^2 dA = {val:.8f}   2pi = {2*np.pi:.8f}   ratio={val/(2*np.pi):.8f}")