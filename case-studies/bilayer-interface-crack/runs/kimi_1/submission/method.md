# Method

## Units and material constants

I used \(E_2=h=1\), hence \(H=2\) and \(E_1=9\). All reported quantities are dimensionless, so this fixes no physical scale. Also,
\[
 \nu_1=\nu_2=\frac13,\qquad \kappa_1=\kappa_2=\frac53,\qquad
 \mu_1=\frac{27}{8},\qquad \mu_2=\frac38 .
\]
Consequently,
\[
 \epsilon=\frac{1}{2\pi}\log\frac{2}{3}
 =-0.0645317762067041\ldots .
\]
The reference length in the phase definition is \(h=1\), so \(h^{i\epsilon}=1\) in these calculations.

## Energy release rates and magnitudes of K

Far behind the tip, the exact plane-strain axial strain is
\(\varepsilon_{11}=\sigma_{11}/\bar E_k\), where
\(\bar E_k=E_k/(1-\nu_k^2)\). Here \(\bar E_1=81/8\) and
\(\bar E_2=9/8\). For steady quasi-static crack advance, the crack-tip energy flux equals the elastic energy per unit length in the detached beams far behind the tip:
\[
 G=\sum_{k=1}^2\int\frac{\sigma_{11}^2}{2\bar E_k}\,dx_2.
\]
The bonded material far ahead is stress-free and the other boundaries are traction-free, so they add no far-field contribution.

Directly integrating the two stated linear stress profiles gives
\[
 \frac{G_AE_2h}{P^2}
 =\frac{4}{81}+\frac{31}{18}
 =\frac{287}{162}
 =1.771604938271605\ldots,
\]
and
\[
 \frac{G_BE_2h^3}{M^2}
 =\frac{16}{27}+\frac23
 =\frac{34}{27}
 =1.2592592592592593\ldots .
\]
The factors include the complete bending energy in layer 2; in case A its moment is \(P(h+H)/2\).

For the stress-intensity convention in the task,
\[
 G=C|K|^2,\qquad
 C=\frac{S}{16\cosh^2(\pi\epsilon)},\qquad
 S=\frac{\kappa_1+1}{\mu_1}+\frac{\kappa_2+1}{\mu_2}
 =\frac{640}{81}.
\]
For this particular epsilon, \(\cosh^2(\pi\epsilon)=25/24\), hence \(C=64/135\). Therefore the reported magnitudes follow directly from the energies:
\[
 \frac{|K_A|\sqrt h}{P}=\sqrt{\frac{G_A}{C}}
 =\sqrt{\frac{1435}{384}}
 =1.9331267849436742\ldots,
\]
\[
 \frac{|K_B|h^{3/2}}{M}=\sqrt{\frac{G_B}{C}}
 =\sqrt{\frac{85}{32}}
 =1.6298006013006623\ldots .
\]
`verify_results.py` also obtains the two energies by independent numerical quadrature and checks these identities.

## Finite-element calculation of the phases

The phases cannot be obtained from the scalar energies. I solved the two load cases separately with a small-strain, plane-strain Q4 finite-element model (`bilayer_fe.py`). The model has the following features.

- Rectangular layers occupy \(-L\le x_1\le L\), \(0\le x_2\le1\), and \(-2\le x_2\le0\), with the crack tip at the origin.
- Interface nodes for \(x_1<0\) are duplicated to make traction-free crack faces. Interface nodes for \(x_1\ge0\) are shared by the layers.
- The outer surfaces and the right end are traction-free. Exact consistent nodal tractions for the stated \(\sigma_{11}\) profiles are applied at the left ends of the detached layers.
- Three displacement restraints at the far right remove only rigid-body motion.
- Two-by-two Gauss integration and the plane-strain constitutive matrix are used. Both load vectors are solved with the same sparse LU factorization.
- The grid is geometrically graded toward the tip in both coordinate directions. The finest run had \(511\times135\times119\) grid lines and 259,073 active displacement unknowns; its smallest tip cell size was \(10^{-5}h\).

The upper-minus-lower crack-face displacement jump has the leading asymptotic form
\[
 \Delta u_2+i\Delta u_1
 =\frac{S}{2\sqrt{2\pi}\cosh(\pi\epsilon)}
 \frac{K}{1+2i\epsilon}r^{1/2+i\epsilon},\qquad r=-x_1>0.
\]
I normalized the finite-element jumps by this expression and extrapolated the resulting \(K(r)\) to \(r=0\). In addition to the constant term, the fits included the next Williams-series terms. For example, the three-column fit was
\[
 K(r)=K+A r^{1/2-i\epsilon}+B r,
\]
and the longer fit added \(r^{3/2-i\epsilon}\) and \(r^2\). Fits excluded \(r\le10^{-3}h\), where the shared-tip Q4 displacement approximation is least reliable, and were repeated over several outer radii.

As an independent extraction, `check_extraction.py` forms
\[
 K(r)=\sqrt{2\pi r}\,(\sigma_{22}+i\sigma_{12})r^{-i\epsilon}
\]
from averaged Q4 stresses immediately above and below the bonded interface, and applies the same series fits. This does not use the crack-opening relation.

Selected converged estimates are:

| extraction | outer fit radius | case A: fitted magnitude, phase (deg) | case B: fitted magnitude, phase (deg) |
|---|---:|---:|---:|
| finest openings, 3 terms | 0.16 | 1.929994, 52.2353 | 1.627141, 4.2648 |
| finest openings, 3 terms | 0.32 | 1.931137, 52.2093 | 1.628979, 4.2562 |
| finest openings, 3 terms | 0.64 | 1.933222, 52.1343 | 1.633559, 4.2264 |
| finest openings, 5 terms | 0.32 | 1.927751, 52.2603 | 1.624932, 4.2675 |
| medium-grid tractions, 3 terms | 0.08 | 1.924535, 52.1023 | 1.626499, 4.2400 |
| medium-grid tractions, 3 terms | 0.16 | 1.927885, 52.0902 | 1.629877, 4.2404 |

I selected the central phase estimates
\[
 \boxed{\psi_A=52.15^\circ},\qquad \boxed{\psi_B=4.24^\circ}.
\]
They are supported by both kinematic quantities (displacements) and static quantities (tractions). The fitted magnitudes agree with the analytical energy-based values to better than 0.45% in all shown checks; the final submitted magnitudes are the analytical values, not individual fit results.

## Combined loading and checks

For \(M=Ph\), the non-dimensional unit choices give \(P=M=1\), and linearity gives \(K=K_A+K_B\). Using the analytical magnitudes and the phases above,
\[
 \psi_{\rm both}=\arg\!\left[
 |K_A|e^{i\psi_A}+|K_B|e^{i\psi_B}\right]
 =30.3611229031^\circ .
\]

The following checks are saved in this folder and can be repeated with `python3 verify_results.py`:

1. Independent quadrature of both far-field energies and verification of the \(G\)-\(K\) identity.
2. Verification of the combined-load phase by direct complex superposition.
3. The opening and independent traction fits shown above.
4. Refinement from 45,783 to 131,015 to 259,073 unknowns.
5. Truncation-length sensitivity: changing \(L/h\) from 12 to 6 changed the representative complex-\(K\) fit by only \(7.2\times10^{-6}\).

The principal commands used were:

```sh
python3 bilayer_fe.py --growth 1.08 --largest-x 0.08 --largest-y 0.06 \
  --first 0.00001 --out run_fine.npz
python3 fit_openings.py run_fine.npz
python3 check_extraction.py
python3 make_results.py
python3 verify_results.py
```

## Remaining numerical uncertainty

The energy release rates and magnitudes are analytical within the stated ideal model. The phase angles are numerical extrapolations, not closed-form results. Based on the shown fit and mesh spread, I estimate the residual numerical uncertainty at roughly \(0.2^\circ\) for case A and \(0.1^\circ\) for case B; the combined phase has comparable or smaller sensitivity. No comparison with an external reference solution was available.
