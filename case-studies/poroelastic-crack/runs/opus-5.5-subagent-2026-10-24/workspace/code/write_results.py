"""Assemble results.json from closed-form expressions and the numerical extrapolations."""
import json, numpy as np
nu, nuu, al = 0.2, 0.4, 0.89
eta = al*(1 - 2*nu)/(2*(1 - nu)); beta = (nuu - nu)/(al*(1 - 2*nu)); gam = (nuu - nu)/(1 - nu)
fx = json.load(open('<opus-workspace>/code/final_extrapolation.json'))
a = fx['a']; A_tip = fx['A_tip']
A_face = 16*gam*(1 - nuu)/np.sqrt(np.pi)
res = {
  "K1_0": round(a, 5),
  "K1_inf": 1.0,
  "V1_0": round(np.pi*(1 - nuu), 6),
  "V1_inf": round(np.pi*(1 - nu), 6),
  "V1_short_coeff": round(A_face + A_tip, 3),
  "p_faces_0": round(2*(nuu - nu)/(al*(1 - 2*nu)), 6),
  "K2_0": 0.0,
  "K2_inf": round(-eta, 6),
  "V2_inf": round(-np.pi*al*(1 - 2*nu)/2, 6),
  "formulas": {
    "definitions": "eta = alpha(1-2nu)/(2(1-nu)); beta = (nu_u-nu)/(alpha(1-2nu)); gamma = (nu_u-nu)/(1-nu) = 2 beta eta; tau = c t/L^2; loading magnitude sigma (loading 1) or p (loading 2)",
    "K1_0": "K1(0+) = a(gamma) sigma sqrt(pi L), gamma=(nu_u-nu)/(1-nu). a = drained-tip/undrained-far-field factor of the self-similar crack-tip problem; Laplace-domain J-integral: (1-nu) a^2 = (1-nu_u) - (fluid tip-flux term) < (1-nu_u); weak coupling a = 1 - 3 gamma/4 + O(gamma^2) [a^2 = (1-gamma)(1-gamma/2) + O(gamma^2)]; numerical fit a = 1 - 0.75 gamma - 0.0634 gamma^2 - 0.0582 gamma^3; a(0.25) = 0.80759. (Instantaneous undrained value at t = 0 exactly: sigma sqrt(pi L).)",
    "K1_inf": "K1(inf) = sigma sqrt(pi L)",
    "V1_0": "V1(0+) = pi (1-nu_u) sigma L^2/G",
    "V1_inf": "V1(inf) = pi (1-nu) sigma L^2/G",
    "V1_short_coeff": "V1 = pi(1-nu_u) sigma L^2/G + A (sigma L^2/G) sqrt(c t/L^2) + O(t^(3/4)); A = (16/sqrt(pi)) gamma (1-nu_u) [drained face layers, exact] + A_tip [crack-tip regions]; A_tip = (2/sqrt(pi))(1-nu_u) chat(gamma), chat = 2 pi gamma + 0.434 gamma^2 + 0.37 gamma^3 (numerical; first order exact); first-order closed form A = 4(4+pi)(nu_u-nu)(1-nu_u)/(sqrt(pi)(1-nu))",
    "p_faces_0": "p(0+) next to faces = 2(nu_u-nu) sigma/(alpha(1-2nu)) = 2 B (1+nu_u) sigma/3",
    "K2_0": "K2(0+) = 0; K2 ~ -C (c t/L^2)^(1/4) p sqrt(pi L) at small t (C = 0.2599 for the given data)",
    "K2_inf": "K2(inf) = -eta p sqrt(pi L) = -alpha(1-2nu) p sqrt(pi L)/(2(1-nu)); approach K2 ~ -eta[1 - 1/(ln(16 c t/L^2) - gamma_E)]",
    "V2_inf": "V2(inf) = -pi eta (1-nu) p L^2/G = -pi alpha (1-2nu) p L^2/(2G)",
  },
  "supplementary": {
    "K1_0_instantaneous_undrained_t_exactly_0": 1.0,
    "K1_small_time": "K1 ~ 0.8076 + 0.1947 tau^(1/4)",
    "V1_short_coeff_face_layers_exact": round(A_face, 6),
    "V1_short_coeff_tip_regions_numerical": round(A_tip, 4),
    "V1_short_coeff_first_order_closed_form": round(4*(4 + np.pi)*gam*(1 - nuu)/np.sqrt(np.pi), 6),
    "V1_short_coeff_naive_uniaxial_layer_only_(incorrect)": round(8*gam/np.sqrt(np.pi), 6),
    "V2_small_time": "V2 ~ -(16/sqrt(pi)) eta (1-nu_u) tau^(1/2) = -1.8077 tau^(1/2)",
    "a_upper_bound_sqrt((1-nu_u)/(1-nu))": round(np.sqrt((1 - nuu)/(1 - nu)), 6),
    "uncertainty": "a = 0.807593 +/- 5e-6; A = 2.4401 +/- 0.0003 (numerical extrapolation s -> infinity, data to s = 1e7)"
  }
}
json.dump(res, open('<opus-workspace>/results.json', 'w'), indent=2, ensure_ascii=False)
print(json.dumps({k: v for k, v in res.items() if k not in ('formulas',)}, indent=1))
