"""Fixed-Talbot numerical Laplace inversion (Abate & Valko 2004)."""
import numpy as np

def talbot_nodes(t, M):
    r = 2.0*M/(5.0*t)
    th = np.arange(1, M)*np.pi/M
    cot = 1.0/np.tan(th)
    S = r*th*(cot + 1j)
    sig = th + (th*cot - 1.0)*cot
    return r, S, sig

def talbot_invert(Fvals_r, Fvals_S, t, M):
    r, S, sig = talbot_nodes(t, M)
    tot = 0.5*np.exp(r*t)*Fvals_r + np.sum(np.real(np.exp(t*S)*Fvals_S*(1 + 1j*sig)))
    return (r/M)*tot

if __name__ == '__main__':
    # tests: F = 1/s^{3/2} -> 2 sqrt(t/pi);  F = 1/(s+1) -> e^-t ; F = exp(-sqrt s)/s -> erfc(1/(2 sqrt t))
    from scipy.special import erfc
    for t in [1e-3, 0.1, 1.0, 10.0]:
        for M in [12, 20]:
            r, S, sig = talbot_nodes(t, M)
            f1 = talbot_invert(r**-1.5, S**-1.5, t, M); e1 = 2*np.sqrt(t/np.pi)
            f2 = talbot_invert(1/(r+1), 1/(S+1), t, M); e2 = np.exp(-t)
            f3 = talbot_invert(np.exp(-np.sqrt(r))/r, np.exp(-np.sqrt(S))/S, t, M); e3 = erfc(1/(2*np.sqrt(t)))
            print(t, M, abs(f1/e1-1), abs(f2/e2-1), abs(f3-e3))
