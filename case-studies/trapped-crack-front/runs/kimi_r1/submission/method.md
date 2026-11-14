# Method

## 1. Coplanar front perturbation

Let \(h(z)=a(z)-\bar a\), and define the finite-part operator

\[
\mathcal R h(z)=\frac{1}{2\pi}\,\mathrm{PV}\!\int_{-\infty}^{\infty}
       \frac{h(z')-h(z)}{(z'-z)^2}\,dz'.
\]

I used the first-order Rice--Gao--Rice result for a **coplanar** perturbation of a
straight half-plane crack front:

\[
\begin{aligned}
\delta K_1&=K_1^0\,\mathcal R h,\\
\delta K_2&=\frac{2}{2-\nu}K_2^0\,\mathcal R h
             -\frac32 K_3^0h_{,z},\\
\delta K_3&=\frac{2(1-\nu)}{2-\nu}K_3^0\,\mathcal R h
             +\frac{2-\nu}{2}K_2^0h_{,z}.
\end{aligned}
\]

These are in-plane front-waviness formulae, not out-of-plane kink formulae.
Mirror symmetry about \(y=0\) excludes mode-I/shear coupling at this order.
The shear signs were checked directly from the stated local frame:

\[
\mathbf n\simeq \mathbf e_x-h_{,z}\mathbf e_z,\qquad
\mathbf t=\mathbf n\times\mathbf e_y
          \simeq \mathbf e_z+h_{,z}\mathbf e_x.
\]

Thus simple frame rotation contributes \(-K_3^0h_{,z}\) to \(\delta K_2\) and
\(+K_2^0h_{,z}\) to \(\delta K_3\); the remaining induced-field contributions
are included in the coefficients above.

For \(h=A\cos(kz)\), \(k=2\pi/\lambda\),

\[
\mathcal R h=-\frac{k}{2}h=-\pi\frac{A}{\lambda}\cos(kz),
\qquad h_{,z}=-2\pi\frac{A}{\lambda}\sin(kz).
\]

Substitution gives the Part-A SIF coefficients. To evaluate the energy
coefficients, put \(r_2=K_2^0/K_1^0=0.6\), \(r_3=K_3^0/K_1^0=0.8\),
\(c_I=1-\nu^2\), \(c_{III}=1+\nu\), and \(c_2=2/(2-\nu)\). Define

\[
D=c_I(1+r_2^2)+c_{III}r_3^2,\qquad
W=c_I\bigl[1+c_2(r_2^2+r_3^2)\bigr],
\]
\[
P=\frac{2\pi W}{D},\qquad
B=\frac{r_2r_3[3c_I-(1+\nu)(2-\nu)]}{D}.
\]

Since \(G^0=DK_1^{02}/E\), substitution into Irwin's relation gives
\(C_G=-P\) and \(S_G=2\pi B\). The shear nonlocal terms combine using
\(c_{III}[2(1-\nu)/(2-\nu)]=c_Ic_2\).

## 2. Equilibrium front

The condition \(G^0=\Gamma_0\) fixes the loading:

\[
\boxed{\frac{K_1^0\sqrt{1-\nu^2}}{\sqrt{E\Gamma_0}}
      =\sqrt{\frac{c_I}{D}}}.
\]

For a general periodic \(f(s)\), define \(|D_s|\) on each harmonic by
\(|D_s|\cos(2\pi ms)=m\cos(2\pi ms)\) and likewise for sine. Because
\(h=\varepsilon\lambda f\) and \(h_{,z}=\varepsilon f_s\),

\[
\begin{aligned}
\kappa_1&=-\pi|D_s|f,\\
\kappa_2&=-\pi c_2r_2|D_s|f-\frac32r_3f_s,\\
\kappa_3&=-\pi\frac{2(1-\nu)}{2-\nu}r_3|D_s|f
            +\frac{2-\nu}{2}r_2f_s.
\end{aligned}
\]

The first-order equilibrium equation is

\[
\mathcal Lf\equiv-P|D_s|f-Bf_s=g(s),
\]

where the zero-mean conditions on \(f\) and \(g\) require no additional
constant loading correction at this order. The toughness has the Fourier
series

\[
g(s)=2\sum_{m=1}^{\infty}q^m\cos(2\pi ms),\qquad q=0.5.
\]

Writing \(f=\sum_{m\ge1}(F_m^c\cos2\pi ms+F_m^s\sin2\pi ms)\), equilibrium
reduces to

\[
P F_m^c+2\pi B F_m^s=-\frac{2q^m}{m},
\qquad P F_m^s-2\pi B F_m^c=0.
\]

Consequently, with \(\Delta=P^2+(2\pi B)^2\),

\[
F_m^c=-\frac{2Pq^m}{m\Delta},\qquad
F_m^s=-\frac{4\pi Bq^m}{m\Delta}.
\]

The skew (sine) part is a consequence of the mixed-mode shear term in
Irwin's relation. I summed 256 Fourier modes; the remaining coefficients are
smaller than \(O(q^{257})\), far below double-precision accuracy here.
Stationary points were found with bracketed root solves of \(f_s=0\) and
\(\kappa_{3,s}=0\) over the full period. The front extremes are

\[
f_{\min}=-0.228242379768\quad(s=0.00999563097133),
\]
\[
f_{\max}= 0.134050828278\quad(s=-0.470091573956),
\]

and the largest \(\kappa_3\) is at \(s=0.0657666104810\).

## 3. Checks

`solve.py` performs all calculations and writes `results.json`. Its checks include:

- direct reconstruction of \(C_G,S_G\) from all six SIF coefficients;
- a 4096-midpoint periodic principal-value quadrature checking
  \(\mathcal R[\cos(2\pi s)]=-\pi\cos(2\pi s)\);
- evaluation of the equilibrium residual using Irwin's relation directly from
  the computed \(\kappa_1,\kappa_2,\kappa_3\), independently of the reduced
  operator \(\mathcal L\), at 4096 points (maximum residual below \(4.5\times10^{-14}\));
- agreement between the direct-energy and reduced-operator expressions;
- zero mean and periodicity of \(f\);
- comparison of all bracketed extrema with a dense period-wide grid.

These checks test the implementation and internal consistency. The elastic
input is the standard first-order coplanar weight-function result; no separate
finite-amplitude or higher-order elastic simulation was attempted.
