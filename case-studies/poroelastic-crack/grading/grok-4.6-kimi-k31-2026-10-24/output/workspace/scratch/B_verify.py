import json
import math


nu = 0.2
nu_u = 0.4
alpha = 0.89
G = 1.0
kappa = 1.0

beta = (nu_u - nu) / (alpha * (1 - 2 * nu))
p0 = 2 * beta

results = {
    "K1_0": 1.0,
    "K1_inf": 1.0,
    "V1_0": math.pi * (1 - nu_u),
    "V1_inf": math.pi * (1 - nu),
    "V1_short_coeff": 8 * (nu_u - nu) / ((1 - nu) * math.sqrt(math.pi)),
    "p_faces_0": p0,
    "K2_0": 0.0,
    "K2_inf": -alpha,
    "V2_inf": -alpha * math.pi * (1 - nu),
    "formulas": {
        "K1_0": "1",
        "K1_inf": "1",
        "V1_0": "pi*(1-nu_u)",
        "V1_inf": "pi*(1-nu)",
        "V1_short_coeff": "8*(nu_u-nu)/((1-nu)*sqrt(pi))",
        "p_faces_0": "2*(nu_u-nu)/(alpha*(1-2*nu))",
        "K2_0": "0",
        "K2_inf": "-alpha",
        "V2_inf": "-alpha*pi*(1-nu)",
    },
}


def close(actual, expected, tolerance=1e-12):
    assert math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance), (
        actual,
        expected,
    )


# Undrained elimination: sigma_xx = sigma_yy = -1 at the crack faces.
sigma_trace = -2.0
close(-beta * sigma_trace, p0)
two_G_zeta = (
    alpha * (1 - 2 * nu) * sigma_trace
    + alpha**2 * (1 - 2 * nu) ** 2 * p0 / (nu_u - nu)
)
close(two_G_zeta, 0.0)
close(nu + alpha * (1 - 2 * nu) * beta, nu_u)

# In the one-dimensional diffusion layer, plane strain and equilibrium give
# Delta_sigma_xx = -alpha*(1-2*nu)/(1-nu)*Delta_p.  Check the resulting
# normal-strain coefficient and the pressure-diffusion coefficient.
dp = -0.123
dsxx = -alpha * (1 - 2 * nu) / (1 - nu) * dp
deyy = (-nu * dsxx + alpha * (1 - 2 * nu) * dp) / (2 * G)
D = alpha * (1 - 2 * nu) / (2 * G * (1 - nu))
close(deyy, D * dp)
pressure_storage = (
    alpha**2
    * (1 - 2 * nu) ** 2
    * (1 - nu_u)
    / ((nu_u - nu) * (1 - nu))
)
c_from_layer = 2 * G * kappa / pressure_storage
c_given = (
    2
    * kappa
    * G
    * (1 - nu)
    * (nu_u - nu)
    / (alpha**2 * (1 - 2 * nu) ** 2 * (1 - nu_u))
)
close(c_from_layer, c_given)

# Integral_0^infinity erfc(eta) d_eta = 1/sqrt(pi), evaluated independently.
n_diffusion = 200_000
eta_max = 10.0
h = eta_max / n_diffusion
erfc_integral = h * sum(math.erfc((i + 0.5) * h) for i in range(n_diffusion))
close(erfc_integral, 1 / math.sqrt(math.pi), 1e-9)

# The total opening is 2*(1-nu_eff)/G * sqrt(L^2-x^2).  With x/L = xi,
# its area normalization is integral_-1^1 2*sqrt(1-xi^2) d_xi = pi.
n_crack = 200_000
h = 2.0 / n_crack
ellipse_area = h * sum(
    2 * math.sqrt(max(0.0, 1 - (-1 + (i + 0.5) * h) ** 2))
    for i in range(n_crack)
)
close(ellipse_area, math.pi, 2e-6)

# Check the boundary-layer displacement calculation against the compact A.
# Per unit crack length, both faces contribute 4*D*p0*sqrt(c*t/pi).
D_p0 = D * p0
close(D_p0, (nu_u - nu) / (G * (1 - nu)))
close(8 * D_p0 / math.sqrt(math.pi), results["V1_short_coeff"] / G)
close(results["V1_inf"] - results["V1_0"], math.pi * (nu_u - nu))
close(results["K2_inf"], -alpha)
close(results["V2_inf"], -alpha * results["V1_inf"])

with open("results.json", "w", encoding="utf-8") as output:
    json.dump(results, output, indent=2)
    output.write("\n")

print(json.dumps(results, indent=2))
print("All constitutive, diffusion, quadrature, and normalization checks passed.")
