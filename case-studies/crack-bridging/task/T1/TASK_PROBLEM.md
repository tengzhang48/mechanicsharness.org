# A finite flaw under opening tension

An infinite isotropic plate in **plane strain** contains a straight, traction-free
slit on $y=0$, from $x=-a$ to $x=a$, with $a=2.40\ \mathrm{mm}$.
A uniform remote tensile stress $\Sigma$ acts normal to the slit. The solid has

\[
E=130\ \mathrm{GPa},\qquad \nu=0.28.
\]

As the two crack fronts advance along the slit line, intact connectors remain
between the separated faces on $a<|x|<a+L$. Their closing normal traction under
monotonic opening is

\[
t(\delta)=160\ \mathrm{MPa}\quad(0<\delta<9.00\ \mu\mathrm m),
\]

and becomes zero irreversibly when the opening first reaches
$\delta_c=9.00\ \mu\mathrm m$. Here $\delta$ is the **full face-to-face gap**,
and $L$ is the added crack length **at each end**. Advancing either outer crack
front also consumes an intrinsic work $\Gamma_0=22.0\ \mathrm{J/m^2}$, separate
from connector work. Area means area of the crack plane, counting a pair of new
faces once. Treat the intrinsic front process as spatially negligible. The bulk
outside the separated region is linear elastic; assume symmetric, quasistatic
growth and no other damage.

Following the monotonic branch from crack initiation, find:

1. The remote stress $\Sigma$, in MPa, when a connector adjacent to the original
   slit first fails.
2. The corresponding $L$, in mm per end.

Give a brief method and one physical or numerical check. Each numerical answer
must be within **±0.3%** of the reference. Budget: 45 minutes; Python with NumPy
and SciPy is available.
