#!/usr/bin/env python3
"""Closed-book viscoelastic strip solution; NumPy/SciPy only.

Run: python solve.py case.json results.json
All bulk-loss values are direct spatial integrals, not energy differences.
"""
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad, quad_vec
from scipy.optimize import brentq
from scipy.special import gammainc
from numpy.polynomial.legendre import leggauss

PI = np.pi


def loading(epsilon):
    return 2 / PI * np.arccos(np.exp(-epsilon / 2))


def bottom_slope(r, epsilon):
    """-U_X at X=-epsilon+r, 0<r<epsilon, log-stable at rear."""
    den = -np.expm1(-epsilon)
    z = -np.expm1(-r) / den
    # Do not form c as sqrt(1-z): in long zones it would lose the
    # exponentially small slope by cancellation for r about 35 or larger.
    c = np.exp(-r/2)*np.sqrt(-np.expm1(-(epsilon-r))/den)
    return 2 / PI * (np.log1p(c) - 0.5 * np.log(z))


def opening(epsilon, V, lam):
    a = loading(epsilon)
    # This representation retains accuracy as V tends to zero.
    end = min(epsilon, 45 * V)
    memory = quad(lambda r: np.exp(-r / V) * bottom_slope(r, epsilon),
                  0, end, epsabs=2e-12, epsrel=2e-12, limit=300)[0]
    return PI * a * a / 2 - (1 - lam) * memory


def limiting_opening(V, lam):
    # Integrating the infinite-zone limiting slope gives the same creep
    # opening. The change r=V*t resolves arbitrarily small V.
    memory = quad(lambda t: np.exp(-t) * (2 / PI) *
                  (np.log1p(np.exp(-V*t/2)) -
                   0.5 * np.log(-np.expm1(-V*t))),
                  0, 45, epsabs=2e-12, epsrel=2e-12, limit=300)[0] * V
    return lam * PI/2 + (1-lam) * (PI/2 - memory)


def critical_speed(lam, delta):
    if delta >= 1:
        return 0.0
    if delta <= lam:
        return None
    target = PI*delta/2
    hi = 1.0
    while limiting_opening(hi, lam) > target:
        hi *= 2
    return brentq(lambda v: limiting_opening(v, lam)-target,
                  1e-8, hi, xtol=1e-11)


def stress_complex(X, Y, epsilon, derivative=False):
    """S=T2+i*T1, analytic in 0<Y<pi. Also supports the top edge.

    Branches are evaluated with q=e^z on the left and e^-z on the right,
    avoiding overflow and retaining the upper-half-strip continuation.
    """
    X, Y = np.broadcast_arrays(np.asarray(X, float), np.asarray(Y, float))
    z = X + 1j*Y
    out = np.empty(z.shape, complex)
    L = np.sqrt(np.expm1(epsilon))
    left = X <= 0
    q = np.exp(z[left])
    B = np.exp(z[left]/2) / np.sqrt(q - 1 + 0j)
    if derivative:
        out[left] = -L/PI * B / (np.exp(epsilon)*q - 1)
    else:
        # Small positive real part is the continuation from Y>0.
        arg = L*B
        arg = np.maximum(arg.real, 0) + 1j*arg.imag
        out[left] = 2/PI * np.arctan(arg)
    p = np.exp(-z[~left])
    if derivative:
        out[~left] = -L/PI * p / (np.sqrt(1-p)*(np.exp(epsilon)-p))
    else:
        out[~left] = 2/PI * np.arctan(L / np.sqrt(1-p))
    # Exact top condition (removes floating-point sin(pi)).
    out = np.where(Y == PI, out.real + 0j, out)
    return out


def stress(X, Y, epsilon):
    if Y == 0:
        if X == 0:
            return np.array([0., 1.])
        if X == -epsilon:
            raise ValueError('Logarithmic rear-edge stress is not finite.')
        if X > 0:
            return np.array([0., 2/PI*np.arctan(
                np.sqrt(np.expm1(epsilon))/np.sqrt(-np.expm1(-X)))])
        if X > -epsilon:
            return np.array([-bottom_slope(X+epsilon, epsilon), 1.])
        c = np.sqrt(np.expm1(epsilon))*np.exp(X/2)/np.sqrt(-np.expm1(X))
        return np.array([-2/PI*np.arctanh(c), 0.])
    S = stress_complex(X, Y, epsilon).item()
    return np.array([S.imag, S.real])


def creep_average(X, Y, epsilon, V):
    """Independent adaptive integration of the future exponential kernel."""
    a = loading(epsilon)
    far = np.array([0., a])
    end = max(35., X+1.)
    breaks = [X]
    # Kernel scales and both interface ends resolve narrow/long histories.
    breaks += [p for p in [-epsilon, 0., X+V, X+8*V, X+35*V]
               if X < p < end]
    breaks += [end]
    breaks = sorted(set(breaks))
    value = np.zeros(2)
    for lo, hi in zip(breaks[:-1], breaks[1:]):
        if (lo-X)/V > 45:
            break
        hi = min(hi, X+45*V)
        value += quad_vec(lambda s: np.exp(-(s-X)/V)/V *
                          (stress(s, Y, epsilon)-far), lo, hi,
                          epsabs=2e-11, epsrel=2e-11)[0]
    return far+value


def strain(X, Y, epsilon, V, lam):
    return lam*stress(X, Y, epsilon)+(1-lam)*creep_average(X,Y,epsilon,V)


def alternate_strain(X,Y,epsilon,V,lam):
    """Direct kernel in material-history variable t, no far-state subtraction."""
    points = sorted(set([0.,45.] + [(p-X)/V for p in [-epsilon,0.]
                                   if 0 < (p-X)/V < 45]))
    q = np.zeros(2)
    for lo,hi in zip(points[:-1],points[1:]):
        q += quad_vec(lambda t: np.exp(-t)*stress(X+V*t,Y,epsilon),
                      lo,hi,epsabs=3e-12,epsrel=3e-12)[0]
    return lam*stress(X,Y,epsilon)+(1-lam)*q


def spatial_grid(epsilon, max_step, log_count, back, ahead):
    """Core grid, with geometric refinement at both interface ends."""
    lo, hi = -epsilon-back, ahead
    x = [np.linspace(lo, hi, int(np.ceil((hi-lo)/max_step))+1)]
    for c in [-epsilon, 0.]:
        r = np.geomspace(1e-8, min(0.45*epsilon, 0.5), log_count)
        x.extend([c-r, np.array([c]), c+r])
    return np.unique(np.concatenate(x))


def spatial_loss(epsilon, V, lam, delta, max_step, ny, log_count,
                 back=32., ahead=32., check_points=None):
    """Integrate (T-E)^2 over X,Y, adding an explicit exponential tail.

    Y=pi*t^2 with Gauss-Legendre quadrature resolves the boundary corner.
    The stress and its analytic X derivative define intervalwise cubic
    Hermite polynomials. Their exponential history integrals are exact.
    This supplies E on the spatial grid without nested quadrature.
    """
    x = spatial_grid(epsilon, max_step, log_count, back, ahead)
    t, wt = leggauss(ny)
    t = (t+1)/2
    wy = wt/2 * 2*PI*t
    y = PI*t*t
    S = stress_complex(x[:, None], y[None, :], epsilon)
    D = stress_complex(x[:, None], y[None, :], epsilon, derivative=True)
    h = np.diff(x)
    u = h/V
    # W_n = integral_0^1 (h/V) exp(-h*s/V) s^n ds.
    weights = np.array([math.factorial(n)*gammainc(n+1,u)/u**n
                        for n in range(4)])
    decay = np.exp(-u)
    Q = np.full(ny, loading(epsilon), complex)
    compare = {}
    if check_points:
        for xp in [-4., -2., 0., 2., 4.]:
            j = int(np.argmin(abs(x-xp)))
            compare[j] = [int(np.argmin(abs(y-PI*f))) for f in [.1,.5,1.]]
    filter_checks = []
    density = np.empty(len(x))
    density[-1] = np.dot(wy, np.abs((1-lam)*(S[-1]-Q))**2)
    # Keep only y-integrated density, avoiding a second large field array.
    for j in range(len(x)-2, -1, -1):
        ds = S[j+1]-S[j]
        b1 = h[j]*D[j]
        b2 = 3*ds-h[j]*(2*D[j]+D[j+1])
        b3 = -2*ds+h[j]*(D[j]+D[j+1])
        Q = (decay[j]*Q + weights[0,j]*S[j] + weights[1,j]*b1 +
             weights[2,j]*b2 + weights[3,j]*b3)
        density[j] = np.dot(wy, np.abs((1-lam)*(S[j]-Q))**2)
        if j in compare:
            for k in compare[j]:
                ev = lam*S[j,k]+(1-lam)*Q[k]
                independent = strain(x[j],y[k],epsilon,V,lam)
                error = np.max(np.abs([ev.imag,ev.real]-independent))
                filter_checks.append({'X':float(x[j]),'Y':float(y[k]),
                                      'maximum_component_error':float(error)})
    core = np.trapezoid(density, x)
    # Left of x_min the stress is O(exp((X+epsilon)/2)); with T=0,
    # Q(X)=Q(x_min)*exp((X-x_min)/V). Integrate that tail analytically.
    tail = (1-lam)**2 * V/2 * np.dot(wy, np.abs(Q)**2)
    factor = 2/(PI*V*(1-lam)*delta)
    result = {'ratio': float(factor*(core+tail)),
            'core_ratio': float(factor*core), 'left_tail_ratio': float(factor*tail),
            'bounds_X': [float(x[0]), float(x[-1])], 'bounds_Y': [0., float(PI)],
            'maximum_X_step': float(np.max(h)), 'minimum_X_step': float(np.min(h)),
            'X_nodes': len(x), 'Y_nodes': ny, 'geometric_nodes_per_side': log_count,
            'left_boundary_max_abs_stress': float(np.max(np.abs(S[0]))),
            'right_boundary_max_abs_stress_minus_a':
                float(np.max(np.abs(S[-1]-loading(epsilon)))),
            'raw_integral_core': float(core), 'raw_integral_tail': float(tail)}
    if check_points:
        result['filter_vs_independent_adaptive_quadrature'] = filter_checks
        result['filter_maximum_component_error'] = max(q['maximum_component_error'] for q in filter_checks)
    return result


def vertical_U(X, epsilon):
    a = loading(epsilon)
    return PI*a-quad(lambda y: stress(X,y,epsilon)[1], 0, PI,
                     epsabs=2e-11, epsrel=2e-11, limit=250)[0]


def verification(epsilon, V, lam, delta, rows):
    points = [(-2., .3*PI), (0., .4*PI), (2., .5*PI),
              (-epsilon-.2, .15*PI)]
    eq, compat, material, local_energy = [], [], [], []
    for d in [2e-3, 1e-3, 5e-4]:
        re, rc, rm, rp = [], [], [], []
        for x,y in points:
            tx = (stress(x+d,y,epsilon)-stress(x-d,y,epsilon))/(2*d)
            ty = (stress(x,y+d,epsilon)-stress(x,y-d,epsilon))/(2*d)
            ex = (strain(x+d,y,epsilon,V,lam)-strain(x-d,y,epsilon,V,lam))/(2*d)
            ey = (strain(x,y+d,epsilon,V,lam)-strain(x,y-d,epsilon,V,lam))/(2*d)
            re.append(abs(tx[0]+ty[1]))
            rc.append(abs(ex[1]-ey[0]))
            residual = (strain(x,y,epsilon,V,lam)-V*ex -
                        stress(x,y,epsilon)+lam*V*tx)
            rm.append(np.max(np.abs(residual)))
            tt = stress(x,y,epsilon)
            ee = strain(x,y,epsilon,V,lam)
            qq = tt-ee
            # phi_D * G_inf*t_c/tau_c^2, calculated two different ways.
            power_minus_free_energy_rate = -V*(np.dot(tt,ex)-np.dot(ee,ex)-
                                        lam/(1-lam)*np.dot(qq,tx-ex))
            square_loss = np.dot(qq,qq)/(1-lam)
            rp.append(abs(power_minus_free_energy_rate-square_loss))
        eq.append({'step':d,'max_absolute_residual':float(max(re))})
        compat.append({'step':d,'max_absolute_residual':float(max(rc))})
        material.append({'step':d,'max_absolute_residual':float(max(rm))})
        local_energy.append({'step':d,'max_absolute_residual':float(max(rp))})
    # A second opening calculation: integrate the stress vertically to
    # reconstruct U, then apply creep to that displacement. This does not
    # use the bottom-slope formula or the integration-by-parts root formula.
    Urear = vertical_U(-epsilon,epsilon)
    end = min(epsilon,45*V)
    history_U = quad(lambda r: np.exp(-r/V)/V * vertical_U(-epsilon+r,epsilon),
                     0,end,epsabs=3e-10,epsrel=3e-10,limit=200)[0]
    Hrear = lam*Urear+(1-lam)*history_U
    front = []
    for r in [1e-2,1e-4,1e-6,1e-8]:
        front.append({'distance':r,'T2_ahead':float(stress(r,0,epsilon)[1]),
                      'T1_cohesive':float(stress(-r,0,epsilon)[0])})
    ahead_openings = [vertical_U(x,epsilon) for x in [.2,2.,5.]]
    interface_traction = [abs(stress(-epsilon*s,0,epsilon)[1]-1)
                          for s in [.1,.3,.7,.9]]
    crack_traction = [abs(stress(-epsilon-r,0,epsilon)[1]) for r in [.1,1.,4.]]
    top_e1 = [abs(strain(x,PI,epsilon,V,lam)[0]) for x in [-4.,0.,4.]]
    relaxed = stress(32., .4*PI,epsilon)-[0,loading(epsilon)]
    free_stress = stress(-epsilon-40., .4*PI,epsilon)
    rear_far_X = -epsilon-max(40.,40*V)
    free_strain = strain(rear_far_X,.4*PI,epsilon,V,lam)
    relaxed_strain = strain(32.,.4*PI,epsilon,V,lam)-[0,loading(epsilon)]
    numeric_top_H = PI*loading(epsilon)*(lam+(1-lam)*quad(lambda t:np.exp(-t),0,45)[0])
    sample_checks = []
    for row in rows:
        e = alternate_strain(row['X'],row['y_over_h']*PI,epsilon,V,lam)
        error = max(abs(e[0]-row['E1']),abs(e[1]-row['E2']))
        sample_checks.append({'X':row['X'],'y_over_h':row['y_over_h'],
                              'maximum_component_error':float(error)})
    return {'normalizations': {
                'equilibrium':'derivatives of T with respect to X,Y; units tau_c*pi/h',
                'compatibility':'derivatives of E with respect to X,Y; units tau_c*pi/(G_inf*h)',
                'material_law':'E-V*E_X-T+lambda*V*T_X; units tau_c/G_inf',
                'displacement':'U=pi*G_inf*w/(tau_c*h); H is its creep average'},
            'equilibrium_finite_difference':eq,
            'strain_compatibility_finite_difference':compat,
            'material_law_finite_difference':material,
            'local_power_minus_free_energy_rate_vs_square_loss':local_energy,
            'cohesive_traction_max_error':float(max(interface_traction)),
            'separated_traction_max_error':float(max(crack_traction)),
            'intact_U_max_abs':float(max(abs(v) for v in ahead_openings)),
            'top_E1_max_abs':float(max(top_e1)),
            'top_H_minus_pi_a':float(numeric_top_H-PI*loading(epsilon)),
            'rear_U_vertical':float(Urear),
            'rear_U_minus_pi_a_squared_over_2':float(Urear-PI*loading(epsilon)**2/2),
            'critical_delta_from_vertical_stress_and_creep':float(2*Hrear/PI),
            'critical_delta_absolute_error':float(abs(2*Hrear/PI-delta)),
            'critical_opening_root_absolute_residual':float(abs(opening(epsilon,V,lam)-PI*delta/2)),
            'bounded_front_sequence':front,
            'far_ahead_T_error':float(np.max(np.abs(relaxed))),
            'far_ahead_E_error':float(np.max(np.abs(relaxed_strain))),
            'far_behind_T_max_abs':float(np.max(np.abs(free_stress))),
            'far_behind_E_check_X':float(rear_far_X),
            'far_behind_E_max_abs':float(np.max(np.abs(free_strain))),
            'sample_field_alternative_kernel_checks':sample_checks,
            'sample_field_max_alternative_kernel_error':max(q['maximum_component_error'] for q in sample_checks)}


def main():
    started = time.time()
    source = Path(sys.argv[1] if len(sys.argv)>1 else 'case.json')
    dest = Path(sys.argv[2] if len(sys.argv)>2 else 'results.json')
    inp = json.loads(source.read_text())
    lam, delta = float(inp['lambda']), float(inp['delta_bar'])
    vc = critical_speed(lam,delta)
    out = {'critical_speed':vc,'cases':[],
           'parameters':{'lambda':lam,'delta_bar':delta},
           'method':'Analytic stress, adaptive creep; direct spatial bulk-loss quadrature',
           'environment':{'python':platform.python_version(),
                          'numpy':np.__version__,'scipy':scipy.__version__}}
    eps0 = -2*np.log(np.cos(PI*np.sqrt(delta)/2))
    out['slow_speed_limit'] = {'epsilon':float(eps0),'a':float(np.sqrt(delta)),
                              'R':1.,'bulk_dissipation_ratio':0.}
    if vc:
        dv = 1e-4
        derivative = (limiting_opening(vc+dv,lam)-limiting_opening(vc-dv,lam))/(2*dv)
        out['endpoint'] = {'a_limit':1.,'R_limit':1/delta,'epsilon_limit':None,
                           'opening_limit_derivative':float(derivative),
                           'bulk_dissipation_ratio_limit_energy_prediction':float(1/delta-1),
                           'length_divergence':'epsilon = -2 log(Vc-V) + 2 log(2/abs(Hinf_prime)) + o(1)'}
        series_checks = []
        for n in [1024,4096,16384]:
            odd = 2*np.arange(n)+1
            hs = lam*PI/2+(1-lam)*4/PI*np.sum(1/(odd**2*(1+odd*vc/2)))
            series_checks.append({'odd_terms':n,'limiting_opening':float(hs),
                                  'absolute_residual':float(abs(hs-PI*delta/2))})
        out['endpoint']['independent_series_checks'] = series_checks
        out['endpoint']['finite_length_sequence_at_Vc'] = [
            {'epsilon':e,'opening':float(opening(e,vc,lam)),
             'opening_shortfall':float(PI*delta/2-opening(e,vc,lam))}
            for e in [8.,16.,24.,32.,40.]]
        out['endpoint']['speed_side_checks'] = [
            {'V':float(vc+d),'limiting_opening_minus_target':
             float(limiting_opening(vc+d,lam)-PI*delta/2)} for d in [-1e-3,1e-3]]
    for v in inp['speeds']:
        V = float(v)
        admissible = V>0 and (vc is None or V<vc)
        case = {'V':V,'admissible':admissible}
        if not admissible:
            case.update({k:None for k in ['epsilon','a','R','bulk_dissipation_ratio']})
            case['fields'] = []
            out['cases'].append(case)
            continue
        hi = max(1.,eps0*2)
        while opening(hi,V,lam) < PI*delta/2:
            hi *= 2
        epsilon = brentq(lambda e:opening(e,V,lam)-PI*delta/2,
                          max(1e-10,eps0/2),hi,xtol=2e-12)
        a = loading(epsilon)
        case.update({'epsilon':float(epsilon),'a':float(a),'R':float(a*a/delta)})
        rows = []
        for x in inp['sample_X']:
            for yh in inp['sample_y_over_h']:
                t = stress(x,PI*yh,epsilon)
                e = strain(x,PI*yh,epsilon,V,lam)
                rows.append({'X':float(x),'y_over_h':float(yh),
                             'T1':float(t[0]),'T2':float(t[1]),
                             'E1':float(e[0]),'E2':float(e[1])})
        case['fields'] = rows
        # Save a reproducible partial result before the more costly evidence.
        case['bulk_dissipation_ratio'] = None
        out['cases'].append(case)
        dest.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
        refinements = []
        for step,ny,ng in [(.02,48,100),(.01,80,160),(.005,128,240),
                            (.0025,192,360)]:
            ev = spatial_loss(epsilon,V,lam,delta,step,ny,ng,
                              check_points=(step==.0025))
            refinements.append(ev)
            print(f'V={V:g}, dx={step:g}, ny={ny}: D/G0={ev["ratio"]:.12g}',flush=True)
        tails = []
        for back,ahead in [(24.,24.),(40.,40.)]:
            tails.append(spatial_loss(epsilon,V,lam,delta,.0025,192,360,back,ahead))
        separate_checks = []
        for label,step,ny,ng in [('X_only',.0025,128,240),
                                  ('Y_only',.005,192,240),
                                  ('corner_only',.005,128,360)]:
            ev = spatial_loss(epsilon,V,lam,delta,step,ny,ng)
            ev['changed_discretization'] = label
            separate_checks.append(ev)
        loss = refinements[-1]['ratio']
        case['bulk_dissipation_ratio'] = loss
        case['bulk_integration'] = {'refinements':refinements,'tail_checks':tails,
                                    'separate_refinement_checks':separate_checks,
                                    'separate_refinement_baseline':{'max_step':.005,'Y_nodes':128,'geometric_nodes_per_side':240},
                                    'density':'2/[pi V (1-lambda) delta_bar] * |T-E|^2 in upper half',
                                    'tail_method':'Analytic exponential left strain tail; negligible stress verified by moving bounds',
                                    'last_refinement_absolute_change':abs(loss-refinements[-2]['ratio']),
                                    'maximum_tail_check_absolute_change':max(abs(q['ratio']-loss) for q in tails)}
        case['energy_balance'] = {'signed_discrepancy':float(loss+1-case['R']),
                                  'relative_discrepancy':float(abs(loss+1-case['R'])/max(1,case['R']))}
        case['verification'] = verification(epsilon,V,lam,delta,rows)
        dest.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    out['runtime_seconds'] = time.time()-started
    out['resources'] = ['PROBLEM.md','CONTRACT.md','case.json',
                        'Python standard library','NumPy','SciPy']
    out['unresolved'] = []
    dest.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(f'Saved {dest}; runtime {out["runtime_seconds"]:.1f} s',flush=True)


if __name__ == '__main__':
    main()
