"""Independent Fourier collocation check of the initial stress intensity."""
import json
import numpy as np
from scipy.linalg import solve, toeplitz

nu, nuu, alpha = .2, .4, .89
eta = alpha*(1-2*nu)/(2*(1-nu))
b = (nuu-nu)/(alpha*(1-2*nu))
Du = 1/(1-nuu)
for n, s in ((2048, 1600.), (2048, 6400.), (4096, 6400.)):
    N = 8*n
    dx = 2/n
    k = np.abs(2*np.pi*np.fft.fftfreq(N, dx))
    m = np.sqrt(k*k+s)
    h = 2*eta*k/(m+k)
    a = 1-4*eta*b*k*k/(m*(m+k))
    symbols = [-Du*k*a, -h/m,
               2*b*Du*k*(m-k)/m, 1/m]
    blocks = [toeplitz(np.fft.ifft(z).real[:n]) for z in symbols]
    mat = np.block([[blocks[0], blocks[1]], [blocks[2], blocks[3]]])
    rhs = np.r_[-np.ones(n), np.zeros(n)]
    sol = solve(mat, rhs)
    x = -1+(np.arange(n)+.5)*dx
    v = sol[:n]
    elastic = solve(toeplitz(np.fft.ifft(-Du*k).real[:n]), -np.ones(n))
    # Drained crack opening: v=(1-nu)*K_normalized*sqrt(1-x*x).
    r = 1-x
    ratio = v/elastic*(1-nuu)/(1-nu)
    samples = {str(j): float(ratio[-j]) for j in (2, 4, 8, 16)}
    print(json.dumps({'n': n, 'ell': float(1/np.sqrt(s)), 'K_from_opening': samples,
                      'volume': float(2*dx*v.sum()),
                      'relative_residual': float(np.linalg.norm(mat@sol-rhs)/np.linalg.norm(rhs))}))
