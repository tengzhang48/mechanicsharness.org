#!/usr/bin/env python3
"""Extract complex interface SIFs from saved crack-face displacement jumps."""
import json
import math
import sys

import numpy as np


def read_k(filename):
    q = json.load(open(filename))
    p = q['parameters']
    r = np.asarray(q['r'])
    eps = p['epsilon']
    D = p['D']
    c = math.cosh(math.pi * eps)
    output = []
    for case in 'AB':
        delta = np.asarray(q['delta' + case + '_real']) + 1j * \
                np.asarray(q['delta' + case + '_imag'])
        K = delta * 2 * (1 + 2j*eps) * c / D * np.sqrt(2*np.pi/r) * r**(-1j*eps)
        output.append(K)
    return q, r, output


def extrapolate(r, K, eps, rmax=0.02, model=1, rmin=0):
    take = (r >= rmin) & (r <= rmax)
    x, y = r[take], K[take]
    powers = [0]
    if model >= 1:
        powers.append(0.5 - 1j*eps)
    if model >= 2:
        powers.append(1.0)
    if model >= 3:
        powers.append(1.5 - 1j*eps)
    if model >= 4:
        powers.append(2.0)
    X = np.column_stack([np.ones(len(x)) if a == 0 else x**a for a in powers])
    coef, residual, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    error = np.max(np.abs(X @ coef-y))
    return coef[0], error, len(x)


if __name__ == '__main__':
    for fn in sys.argv[1:]:
        q, r, Ks = read_k(fn)
        print('\n', fn, q['parameters'])
        inner = max(10*q['parameters']['rmin'], 1e-5)
        for rmax in [0.003, 0.01, 0.02, 0.05]:
            print('rmax', rmax)
            K0s = []
            for case, K in zip('AB', Ks):
                K0, error, number = extrapolate(r, K, q['parameters']['epsilon'],
                                                rmax, 2, inner)
                K0s.append(K0)
                print(case, 'abs %.9f phase %.7f fitmax %.2g n %d' %
                      (abs(K0), np.angle(K0, deg=True), error, number))
            print('both phase %.7f' % np.angle(K0s[0]+K0s[1], deg=True))
