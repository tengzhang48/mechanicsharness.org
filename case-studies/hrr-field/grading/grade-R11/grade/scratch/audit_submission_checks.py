"""Diagnostic only: audit submitted reconstruction using reference angular fields.
This replaces an absent dependency and omits the submission's shooting check.
It is not a successful rerun of the unmodified submission.
"""
import sys
import types
from pathlib import Path
import numpy as np
import reference

hrr = types.ModuleType('hrr')
cases = {}
def get_case(n):
    if n not in cases:
        bvp = reference.collocation(n)
        pars, ivp = reference.shooting(n, [bvp.p[0], bvp.y[2, 0]])
        cases[n] = (pars, ivp, reference.properties(n, pars[0], ivp.sol))
    return cases[n]

def shoot(plane_strain, n):
    assert plane_strain
    pars, ivp, props = get_case(n)
    return pars[1], ivp

def solve_angular(n, plane_strain=True, n_theta=16001):
    assert plane_strain
    pars, ivp, props = get_case(n)
    theta = np.linspace(0, np.pi, n_theta)
    f, fp, fpp, fppp = ivp.sol(theta)
    s = 2-1/(n+1)
    sr, st, srt = s*f+fpp, s*(s-1)*f, (1-s)*fp
    se = np.sqrt(.75*(sr-st)**2+3*srt**2)
    c = 1/se.max()
    return dict(theta=theta, sr=c*sr, st=c*st, se=c*se, f=c*f,
                fp=c*fp, fpp=c*fpp, fppp=c*fppp, c=c, p=pars[1],
                _I=props['J_f0_one']*c**(n+1), _En=np.pi*props['C'])

hrr.singular_exponent = lambda n: 2-1/(n+1)
hrr._shoot = shoot
hrr.solve_angular = solve_angular
hrr.integral_I = lambda n, plane_strain=True, sol=None: sol['_I']
hrr.energy_coefficient_En = lambda n, plane_strain=True, sol=None, I=None: sol['_En']
sys.modules['hrr'] = hrr

source = Path('verify3.py').read_text()
start = source.index('    # --- two-parameter shooting:')
end = source.index('    # --- angular solution', start)
source = source[:start] + '    # Eigenvalue shooting omitted in this diagnostic.\n' + source[end:]
print('DIAGNOSTIC: supplied reference angular solver replaces absent hrr; original eigenvalue check omitted.', flush=True)
exec(compile(source, 'verify3_diagnostic.py', 'exec'))
