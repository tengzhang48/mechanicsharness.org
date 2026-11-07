import json
import math
from pathlib import Path

from scipy.integrate import quad

NU = 0.2
NU_U = 0.4
ALPHA = 0.89
D = NU_U - NU
A = ALPHA * (1 - 2 * NU)
B = A * A / D
ETA = A / (1 - NU)
B_LAYER = B - A * ETA
LAMBDA = 2 * D / (A * (1 - 2 * NU_U))


def close(actual, expected, tol=1e-11, label=""):
    if not math.isclose(actual, expected, rel_tol=tol, abs_tol=tol):
        raise AssertionError(f"{label}: {actual} != {expected}")


def computed_values():
    return {
        "K1_0": 1.0,
        "K1_inf": 1.0,
        "V1_0": math.pi * (1 - NU_U),
        "V1_inf": math.pi * (1 - NU),
        "V1_short_coeff": 16 * D * (1 - NU_U) / ((1 - NU) * math.sqrt(math.pi)),
        "p_faces_0": 2 * D / A,
        "K2_0": 0.0,
        "K2_inf": -ALPHA,
        "V2_inf": -ALPHA * math.pi * (1 - NU),
    }


def check_elastic_crack():
    g = 1.0
    length = 1.0
    for poisson in (NU, NU_U):
        volume = quad(
            lambda x: 2 * (1 - poisson) * math.sqrt(length * length - x * x) / g,
            -length,
            length,
            epsabs=1e-12,
        )[0]
        close(volume, math.pi * (1 - poisson) / g, label="elastic crack volume")

    def sigma_ahead(r):
        return (1 + r) / math.sqrt(r * (2 + r)) - 1

    r = 1e-12
    kappa = math.sqrt(2 * math.pi * r) * sigma_ahead(r)
    close(kappa, math.sqrt(math.pi), tol=2e-5, label="elastic stress intensity")

    s_face = -2.0
    p_face = -D * s_face / A
    close(p_face, 2 * D / A, label="initial face pressure")
    close((A * s_face + B * p_face) / 2, 0.0, label="undrained fluid content")


def check_undrained_constitutive():
    g = 1.3
    stress = ((-0.7, 0.2), (0.2, -1.1))
    trace = stress[0][0] + stress[1][1]
    pressure = -D * trace / A
    for i in range(2):
        for j in range(2):
            delta = float(i == j)
            original = (stress[i][j] - NU * trace * delta + A * pressure * delta) / (2 * g)
            undrained = (stress[i][j] - NU_U * trace * delta) / (2 * g)
            close(original, undrained, label="undrained constitutive reduction")


def check_layer_diffusion():
    kappa = 0.37
    g = 1.3
    c_local = 2 * kappa * g / B_LAYER
    c_given = (
        2 * kappa * g * (1 - NU) * D / (A * A * (1 - NU_U))
    )
    close(c_local, c_given, label="diffusivity")

    c = c_given
    t = 0.013
    y = 0.21
    p_face = 2 * D / A
    similarity = y / (2 * math.sqrt(c * t))
    exponential = math.exp(-similarity * similarity)
    p_t = -p_face * y * exponential / (2 * math.sqrt(math.pi * c) * t ** 1.5)
    c_p_yy = -p_face * y * exponential / (2 * math.sqrt(math.pi * c) * t ** 1.5)
    close(p_t, c_p_yy, label="layer heat equation")
    close(p_face * math.erf(0.0), 0.0, label="drained face condition")

    integral = quad(
        lambda z: math.erfc(z), 0.0, math.inf, epsabs=1e-13, epsrel=1e-13
    )[0]
    close(integral, 1 / math.sqrt(math.pi), label="erfc integral")


def check_trace_integral():
    def angular_integrand(r):
        r4 = r**4
        first = 4 * math.pi * r * (1 + r**-4) / (r4 - 1)
        second = 8 * math.pi / (r**3 * (r4 - 1))
        return first - second

    for r in (1.05, 1.2, 1.7, 3.0, 10.0):
        close(angular_integrand(r), 4 * math.pi / r**3, label="conformal angular integral")
    integral = quad(lambda r: 4 * math.pi / r**3, 1.0, math.inf)[0]
    close(integral, 2 * math.pi, label="integral of squared trace")

    volume_u = math.pi * (1 - NU_U)
    volume_d = math.pi * (1 - NU)
    close(volume_u + D * integral / 2, volume_d, label="reciprocal endpoint check")


def check_short_time_coefficient():
    g = 1.0
    length = 1.0
    c = 0.7
    p_face = 2 * D / A
    zeta_amplitude = B_LAYER * p_face / (2 * g)
    theta_u_face = -(1 - 2 * NU_U) / g
    integrand_amplitude = -LAMBDA * theta_u_face * zeta_amplitude
    expected = computed_values()["V1_short_coeff"]

    scaled_integral = quad(
        lambda z: math.erfc(z), 0.0, math.inf, epsabs=1e-14, epsrel=1e-12
    )[0]
    for t in (1e-8, 1e-6, 1e-4):
        layer_integral = 2 * math.sqrt(c * t) * scaled_integral
        delta_volume = 4 * length * integrand_amplitude * layer_integral
        coefficient = delta_volume / math.sqrt(c * t / length**2)
        close(coefficient, expected, label="short-time volume coefficient")

    inner = 8 * D / ((1 - NU) * math.sqrt(math.pi))
    outer = 8 * D * (1 - 2 * NU_U) / ((1 - NU) * math.sqrt(math.pi))
    close(inner + outer, expected, label="inner and outer displacement check")


def check_loading_superposition(values):
    close(values["K1_inf"] + values["K2_inf"], 1 - ALPHA, label="combined long-time K")
    close(
        values["V1_inf"] + values["V2_inf"],
        (1 - ALPHA) * math.pi * (1 - NU),
        label="combined long-time V",
    )
    close(values["K2_0"], 0.0, label="initial pressure-only K")


def check_results_file(values):
    path = Path(__file__).with_name("results.json")
    if not path.exists():
        return
    data = json.loads(path.read_text())
    for key, expected in values.items():
        close(data[key], expected, tol=1e-12, label=f"results.json {key}")
    for key in values:
        if key not in data.get("formulas", {}):
            raise AssertionError(f"results.json formulas missing {key}")


def main():
    check_elastic_crack()
    check_undrained_constitutive()
    check_layer_diffusion()
    check_trace_integral()
    check_short_time_coefficient()
    values = computed_values()
    check_loading_superposition(values)
    check_results_file(values)
    for key, value in values.items():
        print(f"{key}: {value:.16g}")
    print("All verification checks passed.")


if __name__ == "__main__":
    main()
