import numpy as np

# ---- Build the nonlocal operator from the PROBLEM STATEMENT kernel ----
#  s(x) = 1 - (1/pi) PV Int_0^1 S(x,y) w(y) dy
#  S(x,y) = 1/(x-y) - 1/(x+y) - 2y(x-y)/(x+y)^3      (PV at y=x)
# For hat basis phi_j(y) we need  Int phi_j(y) S(x,y) dy.
# Antiderivatives (w.r.t. y), with z = x fixed:
#   Int 1/(z-y) dy          = -log|z-y|
#   Int -1/(z+y) dy         = -log|z+y|
#   Int -2y(z-y)/(z+y)^3 dy : simplify -2y(z-y)/(z+y)^3.
#     Write -2y(z-y)/(z+y)^3 = 2z(y-z)/(z+y)^3... keep as is and integrate by parts.
#   Note y-z = (z+y) - 2z, so -2y(z-y)/(z+y)^3 = 2y(y-z)/(z+y)^3 = 2y[(z+y)-2z]/(z+y)^3
#        = 2y/(z+y)^2 - 4z y/(z+y)^3.
#   Int 2y/(z+y)^2 dy = 2[ log|z+y| + z/(z+y) ]         (std)
#   Int -4z y/(z+y)^3 dy : y/(z+y)^3 = 1/(z+y)^2 - z/(z+y)^3
#        => -4z[ Int (z+y)^-2 dy - z Int (z+y)^-3 dy ]
#        = -4z[ -1/(z+y) + z/(2(z+y)^2) ] = 4z/(z+y) - 2z^2/(z+y)^2
# So total 3rd-term antiderivative G(y) = 2 log|z+y| + 2z/(z+y) + 4z/(z+y) - 2z^2/(z+y)^2
#                                          = 2 log|z+y| + 6z/(z+y) - 2z^2/(z+y)^2
# Full antiderivative: F(y) = -log|z-y| - log|z+y| + 6z/(z+y) - 2z^2/(z+y)^2   (+const)

def F_state(z, y):
    with np.errstate(divide='ignore', invalid='ignore'):
        return (-np.log(np.abs(z - y) + 1e-300) - np.log(np.abs(z + y) + 1e-300)
                + 6.0 * z / (z + y) - 2.0 * z ** 2 / (z + y) ** 2)


def _S_int(z, y):
    """The statement kernel integrand at a point."""
    return 1.0 / (z - y) - 1.0 / (z + y) - 2.0 * y * (z - y) / (z + y) ** 3


def _num_check_F(z, y, h=1e-7):
    return (F_state(z, y + h) - F_state(z, y - h)) / (2 * h)


if __name__ == "__main__":
    for (z, y) in [(1.2, 0.3), (0.6, 0.15), (0.9, 1.4)]:
        print("z=%.2f y=%.2f  dF/dy=%.8f  S_stmt=%.8f" % (z, y, _num_check_F(z, y), _S_int(z, y)))