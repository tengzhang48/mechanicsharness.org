"""Reference answers of the trapped-front benchmark (closed form), written to reference.json.

Source of the theory: H. Gao and J. R. Rice (1986), "Shear stress intensity factors for a planar crack with slightly
curved front", J. Appl. Mech. 53, 774-778: Eq. (15) gives K1, K2, K3 along a slightly curved front of a half-plane
crack to first order in the front deviation; Eq. (2) is Irwin's relation. With stress intensity factors K_a^0 that are
uniform along a straight front and independent of its position (the task's loading):

  K1(z) = K1^0 + (1/2pi) K1^0 L[a](z)
  K2(z) = K2^0 - 2/(2-nu) K3^0 a'(z) + (2-3nu)/(2-nu) (1/2pi) K2^0 L[a](z)
  K3(z) = K3^0 + 2(1-nu)/(2-nu) K2^0 a'(z) + (2+nu)/(2-nu) (1/2pi) K3^0 L[a](z)
  L[a](z) = PV int (a(z') - a(z))/(z' - z)^2 dz'          (Fourier symbol: -pi |k|)

The slope terms cancel in G (Irwin), so G(z) = G0 + (F/pi) L[a](z) with
  G0 = (1-nu^2)/E [K1^2 + K2^2 + K3^2/(1-nu)],
  F  = (1-nu^2)/E [K1^2 + (2-3nu)/(2-nu) K2^2 + (2+nu)/((2-nu)(1-nu)) K3^2]   (Gao and Rice, Eqs. 19-21).

Part A (warm-up, prescribed front a = abar + A cos(2 pi z/lambda)):
  dK_a = (A/lambda) K1^0 [C_a cos + S_a sin],  dG = (A/lambda) G0 [C_G cos + S_G sin].
Part B (equilibrium G = Gamma with Gamma = Gamma0 [1 + eps g(z/lambda)], g = 2 sum q^n cos(2 pi n s)):
  G0 = Gamma0; each harmonic: -(F/pi) pi (2 pi n/lambda) A_n = Gamma0 eps 2 q^n  ->  A_n = -(eps lambda rho/pi) q^n/n,
  rho = Gamma0/F; summed: a - abar = eps lambda f(s), f(s) = (rho/2pi) ln(1 - 2q cos 2pi s + q^2) (zero mean);
  a'(z) = eps rho S(s), S = 2q sin(2pi s)/(1 - 2q cos 2pi s + q^2); (1/2pi) L[a] = eps (rho/2) g(s).
Normalized outputs: kappa_a(s) = (K_a(z) - K_a^0)/(eps K1^0).
Run: python3 -I reference.py [out.json]   (standard library only)
"""
import json, math, sys

NU = 0.3
R2, R3 = 0.6, 0.8          # K2^0/K1^0, K3^0/K1^0
Q = 0.5                    # toughness profile parameter
S_TABLE = [k / 16 for k in range(16)]


def weights(nu):
    m = (1.0, 1.0, 1.0 / (1.0 - nu))                                      # Irwin
    w = (1.0, (2 - 3 * nu) / (2 - nu), (2 + nu) / ((2 - nu) * (1 - nu)))  # Gao and Rice, Eq. 21
    return m, w


def g_profile(s, q=Q):
    return (1 - q * q) / (1 - 2 * q * math.cos(2 * math.pi * s) + q * q) - 1.0


def slope_profile(s, q=Q):
    return 2 * q * math.sin(2 * math.pi * s) / (1 - 2 * q * math.cos(2 * math.pi * s) + q * q)


def solve(nu=NU, r2=R2, r3=R3, q=Q):
    """All answers (the key). The wrong-answer controls have their own solver in scorer/check_controls.py: the closed
    forms below hold only for the right coefficients, for which the slope terms cancel in G."""
    m, w = weights(nu)
    c = {'slope2': 2 / (2 - nu), 'slope3': 2 * (1 - nu) / (2 - nu), 'pv2': w[1], 'pv3': (2 + nu) / (2 - nu),
         'pv1': 1.0, 'm3': m[2]}
    r = (1.0, r2, r3)
    sum_m = 1 + r2 ** 2 + c['m3'] * r3 ** 2
    sum_w = c['pv1'] + c['pv2'] * r2 ** 2 + c['pv3'] * r3 ** 2 / (1 - nu)
    rho = sum_m / sum_w
    out = {'nu': nu, 'K20_over_K10': r2, 'K30_over_K10': r3, 'q': q}
    # Part A: prescribed cosine front, coefficients of (A/lambda) K1^0 cos / sin
    out['partA'] = {'C1': -math.pi * c['pv1'], 'S1': 0.0,
                    'C2': -math.pi * c['pv2'] * r2, 'S2': 2 * math.pi * c['slope2'] * r3,
                    'C3': -math.pi * c['pv3'] * r3, 'S3': -2 * math.pi * c['slope3'] * r2,
                    'CG': -2 * math.pi / rho, 'SG': 0.0}
    # Part B
    f = [(rho / (2 * math.pi)) * math.log(1 - 2 * q * math.cos(2 * math.pi * s) + q * q) for s in S_TABLE]
    k1 = [c['pv1'] * (rho / 2) * g_profile(s, q) for s in S_TABLE]
    k2 = [c['pv2'] * r2 * (rho / 2) * g_profile(s, q) - c['slope2'] * r3 * rho * slope_profile(s, q) for s in S_TABLE]
    k3 = [c['pv3'] * r3 * (rho / 2) * g_profile(s, q) + c['slope3'] * r2 * rho * slope_profile(s, q) for s in S_TABLE]

    def kappa3(s):
        return c['pv3'] * r3 * (rho / 2) * g_profile(s, q) + c['slope3'] * r2 * rho * slope_profile(s, q)
    # maximum of kappa3 on one period: fine scan, then golden-section refinement
    n = 20000
    j = max(range(n), key=lambda i: kappa3(i / n - 0.5))
    lo, hi = (j - 1) / n - 0.5, (j + 1) / n - 0.5
    for _ in range(200):
        a1, a2 = lo + (hi - lo) * 0.381966, lo + (hi - lo) * 0.618034
        if kappa3(a1) < kappa3(a2):
            lo = a1
        else:
            hi = a2
    s_max = 0.5 * (lo + hi)
    out['partB'] = {'K10_scaled': 1 / math.sqrt(sum_m),                 # K1^0 sqrt((1-nu^2)/(E Gamma0))
                    'front_peak_to_peak': (rho / math.pi) * math.log((1 + q) / (1 - q)),
                    's_table': S_TABLE, 'f': f, 'kappa1': k1, 'kappa2': k2, 'kappa3': k3,
                    'kappa3_max': kappa3(s_max), 's_kappa3_max': s_max}
    out['derived'] = {'rho_Gamma0_over_F': rho, 'sum_irwin': sum_m, 'sum_F': sum_w}
    return out


if __name__ == '__main__':
    res = solve()
    path = sys.argv[1] if len(sys.argv) > 1 else 'reference.json'
    json.dump(res, open(path, 'w'), indent=1)
    A, B = res['partA'], res['partB']
    print('rho = Gamma0/F = %.9f' % res['derived']['rho_Gamma0_over_F'])
    print('Part A:', ', '.join('%s %.6f' % kv for kv in A.items()))
    print('Part B: K10_scaled %.6f  peak-to-peak %.6f  kappa3_max %.6f at s = %.6f' %
          (B['K10_scaled'], B['front_peak_to_peak'], B['kappa3_max'], B['s_kappa3_max']))
    for i, s in enumerate(S_TABLE):
        print('  s %.4f  f % .6f  k1 % .6f  k2 % .6f  k3 % .6f' % (s, B['f'][i], B['kappa1'][i], B['kappa2'][i], B['kappa3'][i]))
