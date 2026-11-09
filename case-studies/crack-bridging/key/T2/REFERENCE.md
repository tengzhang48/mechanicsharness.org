# Reference solution

At the first complete connector failure:

| Quantity | Reference | Allowed relative error |
|---|---:|---:|
| Total couple $M$ on each strip | **0.2126257874 N·m** | ±0.5% |
| Retained-zone length $L$ | **10.80010237 mm** | ±0.5% |
| Full opening $\delta(L/2)$ | **38.21889793 µm** | ±0.5% |

Use the unrounded values in [reference_results.json](reference_results.json).
The standalone program [reference.py](reference.py) requires NumPy and SciPy:

    python T2/reference.py

## Adapted mechanical description

For bending stiffness per unit width, use

\[
D_i=\frac{E_i h_i^3}{12},\qquad
C=\frac1{D_1}+\frac1{D_2}.
\]

Numerically, $D_1=1.980000000\ \mathrm{N\,m}$,
$D_2=11.54250000\ \mathrm{N\,m}$ and
$C=0.5916868490\ (\mathrm{N\,m})^{-1}$. Both arms contribute to the relative
opening. Averaging the moduli or thicknesses would not give the correct stiffness.
Define $m=M/b$, the couple per unit width on one arm.

Equilibrium of the two arms gives the opening equation

\[
\delta''''=-Ct(\delta).
\]

At the advancing front the faces meet with the same slope:
$\delta(0)=\delta'(0)=0$. Intrinsic front work sets the relative curvature:

\[
\delta''(0)=\sqrt{2C\Gamma_0}.
\]

This follows because the two equal and opposite tip couples have combined
energy release $\tfrac12 C m_{\rm tip}^2=\Gamma_0$, and
$\delta''(0)=C m_{\rm tip}$. Setting this curvature to zero would incorrectly
discard the supplied front work.

At the original cut root, first complete failure gives $\delta(L)=\delta_c$.
No end force means zero shear there, $\delta'''(L)=0$; the unbridged starter cut
transmits the pure couple, giving $\delta''(L)=Cm$.
Thus a fourth-order boundary-value problem with unknown $L$ determines the
profile; its rear curvature determines $M$.

An energy first integral provides a second way to obtain the couple:

\[
\frac{Cm^2}{2}
 =\Gamma_0+\int_0^{\delta_c}t(\delta)\,d\delta
 =\Gamma_0+\frac{\sigma_b\delta_c}{2}
 =21.40\ \mathrm{J/m^2}.
\]

Consequently $m=8.505031495\ \mathrm N$ and $M=bm=0.2126257874\ \mathrm{N\,m}$.
This energy relation alone cannot determine $L$ or the midpoint opening.
The 20 mm starter-cut length does not enter the retained-zone solution under
pure couples and zero shear.

For numerical solution let

\[
\ell=\left(\frac{\delta_c}{C\sigma_b}\right)^{1/4}
 =6.411756501\ \mathrm{mm},\quad
u=\frac{\delta}{\delta_c},\quad z=\frac x\ell,\quad
\gamma=\frac{\Gamma_0}{\sigma_b\delta_c}=0.035,\quad
\kappa=\sqrt{2\gamma}.
\]

Then $u''''=u-1$, with $u(0)=u'(0)=0$, $u''(0)=\kappa$,
$u(\lambda)=1$, and $u'''(\lambda)=0$, where $\lambda=L/\ell$.
The primary calculation uses SciPy's boundary-value solver on $x/L\in[0,1]$
with $\lambda$ as an unknown parameter.

## Independent calculation and achieved agreement

A separate exact solution is

\[
u(z)=1-\frac{\cosh z+\cos z}{2}
       +\frac{\kappa}{2}(\cosh z-\cos z)
       +\frac A2(\sinh z-\sin z),
\]

where imposing the two rear conditions gives the characteristic equation

\[
1+\cosh\lambda\cos\lambda-\kappa\sinh\lambda\sin\lambda=0
\]

and

\[
A=\frac{\sinh\lambda+\sin\lambda-\kappa(\sinh\lambda-\sin\lambda)}
         {\cosh\lambda+\cos\lambda}.
\]

The smallest positive, physically admissible root is
$\lambda=1.68442179119909$, and $A=1.18543116372457$.
This route uses a scalar root solve and the exact profile, plus the energy
formula for $M$, rather than the boundary-value discretization.

| Output | Boundary-value route | Exact/energy route |
|---|---:|---:|
| $M$ (N·m per arm) | 0.212625787376721 | 0.212625787376721 |
| $L$ (mm) | 10.8001023708973 | 10.8001023708972 |
| $\delta(L/2)$ (µm) | 38.2188979260933 | 38.2188979260938 |

Relative differences are respectively $2.2\times10^{-16}$, $1.2\times10^{-14}$
and $1.4\times10^{-14}$. Maximum full-profile disagreement is
$5.9\times10^{-14}$ in $\delta/\delta_c$. The solver's maximum RMS ODE
residual is $1.0\times10^{-10}$.

Additional executed checks:

- The full-zone first integral
  $\tfrac12u''^2-u'u'''=\gamma+u-u^2/2$ has maximum absolute residual
  $6.5\times10^{-14}$; external energy balance agrees within
  $3.4\times10^{-16}$ relative error.
- At $\Gamma_0=0$, the softening problem recovers the known first cantilever
  characteristic root $1.875104068712$; relative length error is
  $5.8\times10^{-15}$.
- Replacing the law by constant closing traction and taking $\Gamma_0=0$
  gives $\delta_c=C\sigma_bL^4/8$, hence $L/\ell=8^{1/4}$; relative
  length error is $2.2\times10^{-16}$.
- The computed opening and slope are nonnegative and increasing. The maximum
  opening occurs at the rear, where the cutoff is first reached. Work increases
  with rear opening, so the couple increases up to that event even though the
  local traction softens.

The result has $L/h_2=12.0001$ and rear relative rotation $0.04338$ rad,
consistent with the prescribed strip idealization. All reference outputs were
independently checked.
