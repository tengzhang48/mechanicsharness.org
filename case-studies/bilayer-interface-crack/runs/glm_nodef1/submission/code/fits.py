"""Extraction of interface-crack quantities from FEM fields.

Conventions used (standard interface fracture):
    T(r) = sigma22 + i sigma12 = K r^{i eps} / sqrt(2 pi r)  ahead of tip,
with lengths made dimensionless by the reference length h (so at r = h,
arg T = arg K).  eps is measured from the log-slope of arg T.
"""
import numpy as np


def fit_series(r, y, powers):
    """Complex least squares y ~ sum c_p r^p.  Returns complex coeffs."""
    A = np.stack([r**p for p in powers], axis=1)
    c, *_ = np.linalg.lstsq(A, y.astype(complex), rcond=None)
    return c


def residual(r, y, powers):
    c = fit_series(r, y, powers)
    A = np.stack([r**p for p in powers], axis=1)
    model = A @ c
    return np.sum(np.abs(y-model)**2)/np.sum(np.abs(y)**2), c


def measure_eps(r, T, span=(-0.15, 0.15), npts=1201):
    """Scan eps: minimize residual of T sqrt(2 pi r) r^{-i eps} ~ K + c r^1/2 ..."""
    y0 = T*np.sqrt(2*np.pi*r)
    lnr = np.log(r)
    best = None
    for eps in np.linspace(*span, npts):
        y = y0*np.exp(-1j*eps*lnr)
        res, c = residual(r, y, (0.0, 0.5, 1.0, 1.5))
        if best is None or res < best[0]:
            best = (res, eps, c)
    return best[1], best[0], best[2]


def K_from_traction(r, T, eps, powers=(0.0, 0.5, 1.0)):
    y = T*np.sqrt(2*np.pi*r)*np.exp(-1j*eps*np.log(r))
    c = fit_series(r, y, powers)
    return c[0]


def amplitude_from_delta(r, d, eps, powers=(0.0, 0.5, 1.0)):
    """d(r) = A r^{1/2 + i eps}(1+...): return A."""
    y = d/np.sqrt(r)*np.exp(-1j*eps*np.log(r))
    c = fit_series(r, y, powers)
    return c[0]


def measure_eps_delta(r, d, span=(-0.15, 0.15), npts=1201):
    y0 = d/np.sqrt(r)
    lnr = np.log(r)
    best = None
    for eps in np.linspace(*span, npts):
        y = y0*np.exp(-1j*eps*lnr)
        res, c = residual(r, y, (0.0, 0.5, 1.0, 1.5))
        if best is None or res < best[0]:
            best = (res, eps, c)
    return best[1], best[0], best[2]


def measure_eps_unwrap(r, y_raw, sub=1.0):
    """Direct log-slope of unwrapped arg(y_raw * r^{1/2 sub correction}).

    y_raw = T*sqrt(2 pi r) (sub=1) or delta/sqrt(r) (sub=-1 handled by caller
    passing y_raw accordingly).  Fits arg ~ a + eps*ln r + c*sqrt(r) with
    arg unwrapped.
    """
    lnr = np.log(r)
    phi = np.unwrap(np.angle(y_raw))
    A = np.stack([np.ones_like(r), lnr, np.sqrt(r), r], axis=1)
    c, *_ = np.linalg.lstsq(A, phi, rcond=None)
    return c[1], c
