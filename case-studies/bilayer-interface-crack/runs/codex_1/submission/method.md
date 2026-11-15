I used exact beam energies for the energy release rates and stress-intensity
magnitudes, and a two-dimensional plane-strain finite-element solution for the
phases. All calculations were local, with NumPy and SciPy; no web resources or
external data were used.

I set h = E2 = 1, so H = 2 and E1 = 9. The two right-hand sides are (P,M) =
(1,0) and (0,1). In these units the requested phase is simply arg(K), since
h^(i eps) = 1. This choice also makes the computed G and |K| equal to their
requested dimensionless normalizations. The crack faces remain traction-free,
including where this allows interpenetration, as specified in the question.

For the given elastic constants, kappa1 = kappa2 = 5/3, mu1 = 27/8, and
mu2 = 3/8. Substitution in the prescribed definition gives

\[
\epsilon=\frac{\ln(2/3)}{2\pi}=-0.0645317762067041.
\]

Write Ebar_k = E_k/(1-nu_k^2). In a remote separated layer the stress state
has only sigma11, and its strain energy density is sigma11^2/(2 Ebar_k).
Because the bonded state far ahead is unstressed, the steady energy release
rate is the total remote strain energy per unit beam length:

\[
G=\frac{1}{2\overline E_1}
  \left(\frac{P^2}{h}+\frac{12M^2}{h^3}\right)
 +\frac{1}{2\overline E_2}
  \left(\frac{P^2}{H}
       +\frac{12[M+P(h+H)/2]^2}{H^3}\right).
\]

This also follows by evaluating the J integral on a contour far behind and
ahead of the tip. The rear cross section contributes the beam strain energy;
the forward cross section and free surfaces contribute zero. Consequently,

\[
G_A=\frac{287}{162},\qquad G_B=\frac{34}{27},\qquad
G_{P=M=1}=\frac{815}{162}=G_A+G_B+2.
\]

For the traction definition of K in the question, the plane-strain
interface-crack energy relation is

\[
G=C|K|^2,\qquad
C=\frac{(\kappa_1+1)/\mu_1+(\kappa_2+1)/\mu_2}
        {16\cosh^2(\pi\epsilon)}=\frac{64}{135}.
\]

The reported magnitudes therefore come from the exact expressions
|K_A| = sqrt(1435/384) = 1.9331267849436742 and
|K_B| = sqrt(85/32) = 1.629800601300662.

The numerical model in `solve.py` uses nine-node isoparametric quadrilateral
elements for the displacement formulation. The plane-strain material matrix,
with engineering shear strain, is

\[
D_k=E_k\begin{pmatrix}
1.5&0.75&0\\0.75&1.5&0\\0&0&0.375
\end{pmatrix}.
\]

The domain is [-L,L] x [-2,1]. Nodes on y=0 are shared for x>=0 and duplicated
for x<0. Each left end receives traction (-sigma11,0) from the prescribed
remote stress distribution. The top, bottom, crack faces, and right end have
zero prescribed traction. Three scalar displacement constraints at the right
end remove rigid translations and rotation. The applied forces and moments
balance, so these constraints introduce no physical loading; numerical
reactions are included in the reported residual check.

Inside [-1,1] x [-1,1], the mesh consists of geometrically scaled square rings.
There are `nside` elements along each side of a ring. The first ring connects
to the common crack-tip node through elements with collapsed inner edges.
Their Gauss integration points have positive Jacobians. Thus the model has
an actual crack tip, rather than a small hole. The first positive radius is
`rmin`; subsequent radii have an approximately constant ratio until radius 1.
A rectangular mesh covers the remaining domain. The finest calculation used
nside=28, radius ratio about 1.10, rmin=1e-8, L=24, and 5x5 Gauss integration.
It had 25,508 elements and 102,685 displacement nodes.

I extracted K primarily from the displacement jump behind the tip. Define

\[
\delta(r)=[u_2](-r)+i[u_1](-r),\qquad
S=\frac1{\overline E_1}+\frac1{\overline E_2}.
\]

The brackets mean upper-layer displacement minus lower-layer displacement.
The leading crack-opening field, with the traction convention in the task, is

\[
\delta(r)=\frac{4S}{\cosh(\pi\epsilon)}
 \frac{K}{1+2i\epsilon}\sqrt{\frac{r}{2\pi}}r^{i\epsilon}.
\]

I inverted this expression at crack-face nodes, then fitted the resulting
complex K estimates to K0 + a r + b r^2 over 1e-4 <= r <= 1e-2. The phase of
K0 is the reported phase for each individual load. The positive real opening
prefactor affects magnitudes but not phases. As a separate check of the
complex convention and phase sign, I extracted
(sigma22+i sigma12) sqrt(2 pi r) r^(-i eps) on the bonded interface, evaluating
tractions independently from each material.

For M=Ph, linearity gives K_both = K_A + K_B in the chosen units. I used the
exact magnitudes and the fine-mesh individual phases in this sum. The result
is psi_both = 30.4332583 degrees, reported as 30.4333 degrees. Adding or
averaging the individual phase angles would not give this result.

The following checks passed; detailed numbers are saved in `checks.json`.

| Mesh | Nodes | psi A (degrees) | psi B (degrees) | psi of fitted K_A+K_B (degrees) |
| --- | ---: | ---: | ---: | ---: |
| nside=12, ratio=1.25, L=16 | 18,257 | 52.250795 | 4.277864 | 30.434086 |
| nside=20, ratio=1.15, L=24 | 52,749 | 52.249667 | 4.278676 | 30.433515 |
| nside=28, ratio=1.10, L=24 | 102,685 | 52.249396 | 4.278694 | 30.433307 |

- The medium-to-fine phase changes were below 0.0003 degrees. Reducing L
  from 24 to 12 on the medium mesh changed all three phases by less than
  0.000001 degrees.
- Moving the fine-mesh extrapolation interval to [1e-5,1e-3] changed phases
  by at most 0.00055 degrees. The alternative [1e-3,0.03] interval is also
  recorded in the FEM summary files.
- The fine-mesh opening-derived magnitudes differed from the exact energy
  values by less than 0.0006 percent.
- The exact energy cross term requires
  cos(psi_A-psi_B) = 1/sqrt(G_A G_B), so the magnitude of the phase difference
  must be 47.97051474 degrees. The fine-mesh phases differ by
  47.97070200 degrees, an error of 0.00019 degrees.
- Direct traction extrapolation from the upper material gave phases
  52.24162 and 4.33705 degrees; from the lower material it gave 52.25725 and
  4.31626 degrees. These converge toward the opening-derived phases and
  confirm the prescribed normal-plus-i-shear convention. The largest
  fine-mesh traction-derived magnitude error was 0.082 percent.
- J integrals on square contours with radii between 0.001 and 0.5 agreed
  with the exact energies within 0.087 percent on the fine mesh.
- Resultant applied forces and moments balanced to approximately 1e-15.
  The largest fine-mesh linear-system residual, including constrained
  degrees of freedom, was 1.1e-7 for loads normalized to unity.

To reproduce all saved calculations and regenerate `results.json` and
`checks.json`, run the following in this folder:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 solve.py --nside 12 --ratio 1.25 --rmin 1e-7 --length 16 --output fem_coarse
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 solve.py --nside 20 --ratio 1.15 --rmin 1e-8 --length 24 --output fem_medium
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 solve.py --nside 28 --ratio 1.10 --rmin 1e-8 --length 24 --quadrature 5 --output fem_fine
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 solve.py --nside 20 --ratio 1.15 --rmin 1e-8 --length 12 --output fem_short
OPENBLAS_NUM_THREADS=1 python3 checks.py
```

The `.npz` files retain crack-opening estimates, interface-traction estimates,
and contour J integrals, so the checks can be repeated without refactoring
the stiffness matrices. The individual `fem_*.json` files record parameters,
fits, force balance, and residuals. The calculation used Python 3.12.3,
NumPy 1.26.4, and SciPy 1.11.4.

The energies, epsilon, and the energy-derived magnitudes have analytic
expressions above. The absolute phases have numerical convergence evidence,
independent traction and energy checks, and domain-size checks; there is no
rigorous numerical error bound or comparison with an external reference
phase solution. Effects excluded by the stated linear, traction-free model,
including crack-face contact, were not evaluated.
