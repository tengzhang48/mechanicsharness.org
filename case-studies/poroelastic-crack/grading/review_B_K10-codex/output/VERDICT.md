# Verdict: the initial tip intensity and the small-time volume

For the stated permeable, stationary Biot crack, **MP-37 and MP-38 are correct on both contested quantities**. MP-39 is consistent with the correct intensity but is wrong on the leading volume correction.

Let \(\tau=ct/L^2\), and use the full opening \(w=u_y(x,0^+)-u_y(x,0^-)\). The quantities settled here are

\[
k_0=\frac{1}{\sqrt{\pi L}}\lim_{t\downarrow0}
       \lim_{r\downarrow0}\sqrt{2\pi r}\,\sigma_{yy}(L+r,0,t),
\qquad
\frac{GV_1(t)}{L^2}=\pi(1-\nu_u)+A\sqrt{\tau}+o(\sqrt{\tau}).
\]

| Quantity | My value for the specified material | Estimated numerical uncertainty | Side supported |
|---|---:|---:|---|
| \(k_0\) | **0.8075930** | \(\pm 0.000002\) | GPT: MP-37, MP-38; consistent with MP-39's provisional estimate |
| \(A\) | **2.4402196** | \(\pm 0.00002\) | GPT: MP-37 and MP-38; against the paper and MP-39 |

These are numerical error estimates, **not certified interval bounds**. The factor calculation agrees much more closely across methods and refinements than these estimates. A separate finite-crack calculation validates the model and the disputed size of the effects, at lower precision.

The decisive points are:

* The drained inner tip opening has an amplitude **1.07679067 times** its matched undrained outer amplitude. The paper's Eq. (29) effectively sets this matching factor to one.
* The volume coefficient is the sum of **1.3540550 from the flat faces** and **1.0861646 from the two tips**. The latter survives at order \(\sqrt{\tau}\), even though the tip region shrinks.

The reversed stress limit, taking time to zero at a fixed positive distance first, is exactly **1**. The paper's **0.75 is neither of these two stress limits**.

All evidence used is local: the supplied problem, submissions, and paper. The original files were preserved; submitted programs were run on copies in [scratch/](scratch/). A machine-readable numerical summary is [scratch/summary.json](scratch/summary.json).

## 1. Independent reduction of the governing equations

The following derivation starts from PROBLEM.md, rather than assuming the paper's asymptotic formulas or a moving-crack result. To avoid confusing the paper's poroelastic coefficient with the coupling parameter used in MP-37, define

\[
d=\nu_u-\nu,\quad a=\alpha(1-2\nu),\quad b=\frac{a}{1-\nu},
\quad D_f=\frac{a^2(1-\nu_u)}{d(1-\nu)},\quad
\beta=\frac{a}{D_f},\quad
\gamma=b\beta=\frac{d}{1-\nu_u},\quad
C=\frac{G}{2(1-\nu)}.
\]

Write \(T=\sigma_{xx}+\sigma_{yy}\), and \(h=T+bp\). Compatibility of the supplied strain law, together with equilibrium, gives

\[
\nabla^2h=0,\qquad
2G\zeta=D_f(p+\beta h).
\]

For a Laplace parameter \(s>0\), let \(q^2=s/c\). Mass conservation then gives

\[
\nabla^2p=q^2(p+\beta h),\qquad
(\nabla^2-q^2)(p+\beta h)=0.
\]

These equations include the initial zero fluid-content condition through the Laplace transform.

Use the Fourier convention \(\widehat f(\xi)=\int f(x)e^{i\xi x}\,dx\), and set \(k=|\xi|\), \(m=\sqrt{k^2+q^2}\). Decaying upper-half-space modes can be written

\[
h=H e^{-ky},\quad p=P e^{-my}-\beta H e^{-ky},
\]
\[
\Phi=(A_0+B_0y)e^{-ky}+C_0e^{-my},
\quad B_0=-\frac{(1+\gamma)H}{2k},\quad
C_0=-\frac{bP}{q^2},\quad A_0=\frac{B_0-mC_0}{k}.
\]

Here \(\Phi\) is an Airy stress function. The last relation enforces zero shear traction. Integrating the \(yy\) strain from infinity gives

\[
H=-\frac{Gk}{1-\nu}\widehat w.
\]

Consequently the exact surface relations are

\[
\widehat S=-Ck\left[1+\gamma\frac{m-k}{m+k}\right]\widehat w
              -\frac{bk}{m+k}\widehat p_0,
\]
\[
\widehat Q=-m\widehat p_0+2\beta Ck(m-k)\widehat w,
\qquad S=\sigma_{yy}(x,0^+),\quad Q=p_{,y}(x,0^+).
\tag{1}
\]

On the crack, \(S=-1,\ p_0=0\) for the unit resolvent problem. On its intact continuation, symmetry requires **\(w=Q=0\)**. Pressure on that continuation is an output, not an imposed zero boundary value. Here \(Q\) is a pressure derivative; Darcy flux is \(-\kappa Q\).

Eliminate \(p_0\), put \(\lambda=2\beta Cq^2/b\), and use the unknown \(F=-Q/\sqrt{\lambda}\). The symmetric physical boundary operator is

\[
\binom{-S}{p_0/\sqrt{\lambda}}=
\begin{pmatrix}
CkD & \sqrt{2\gamma C}\,q\,k/[m(m+k)]\\
\sqrt{2\gamma C}\,q\,k/[m(m+k)] & 1/m
\end{pmatrix}
\binom{w}{F},
\quad
D=1+\gamma-\frac{2\gamma k^2}{m(m+k)}.
\tag{2}
\]

This independently derived operator is the basis of my finite-crack code. Its Airy-mode algebra was checked symbolically: stress trace, shear traction, diffusion, displacement, normal stress, and pressure-gradient residuals all vanish. See [analytic_checks.py](scratch/analytic_checks.py) and [its output](scratch/analytic_checks.jsonl).

## 2. Matching the initial stress intensity

Away from a shrinking diffusion region, the initial response is undrained elasticity:

\[
w_{\rm out}(x)=\frac{2(1-\nu_u)}{G}\sqrt{L^2-x^2},
\quad K_{\rm out}=\sqrt{\pi L}.
\]

At a positive time, arbitrarily near a permeable tip the pressure tends to zero, and the singular mechanical field is drained. Its opening relation is

\[
w(L-r,t)\sim \frac{4(1-\nu)}{G\sqrt{2\pi}}K(t)\sqrt r.
\tag{3}
\]

The missing operation in the paper is matching these two fields through the region of size \(\delta_s=\sqrt{c/s}\). A pointwise limit of the outer opening does not determine its inner square-root coefficient.

### The matched tip problem

Scale \(X=(x-L)/\delta_s\), normalize the opening so that \(W(X)\sim\sqrt{-X}\) as \(X\to-\infty\), and retain the leading singular outer field. The finite face traction is lower order in this scaling. Opening and flux are supported on \(X<0\); the unforced leading stress and face pressure are supported on \(X>0\).

After rescaling (2), its mechanical and hydraulic diagonal factors are \(k\) and \(1/m\), now with \(m=\sqrt{1+\xi^2}\). Introduce

\[
k_\pm=\sqrt{\mp i\xi},\quad m_\pm=\sqrt{1\mp i\xi},\quad
t_\xi=\frac{|\xi|}{\sqrt{1+\xi^2}+|\xi|}.
\]

Square roots have positive real part. Factoring these known scalar multipliers leaves the matrix

\[
N(\xi)=
\begin{pmatrix}
1+\gamma-2\gamma |\xi|t_\xi/\sqrt{1+\xi^2}
&
\sqrt{2\gamma}\,t_\xi/(m_+k_+)\\
\sqrt{2\gamma}\,t_\xi/(m_-k_-)&1
\end{pmatrix}.
\tag{4}
\]

In particular,

\[
N(0)=\operatorname{diag}(1+\gamma,1),\quad N(\infty)=I,\quad
\det N=1+\gamma\frac{m-k}{m+k}>0.
\]

It is Hermitian positive definite on the real line and has a canonical factorization

\[
N=N_+N_-,\qquad N_\pm(\infty)=I,
\]

with each factor and inverse analytic in its corresponding half-plane. Define

\[
\binom{n_{11}}{n_{21}}=N_-(0)\binom10.
\]

The transform of the outer opening is
\(\Gamma(3/2)(i\xi)^{-3/2}\). Multiplication by \(k_-\) produces its simple pole. Analytic continuation of the Wiener–Hopf equation therefore has pole residue

\[
\Gamma(3/2)N_-(0)(1,0)^T.
\]

At large frequency \(N_-\to I\). Dividing the mechanical component by \(k_-\) shows that the inner opening is

\[
W(X)\sim n_{11}\sqrt{-X}\quad(X\uparrow0).
\]

The hydraulic component similarly gives the tip pressure

\[
P_{\rm tip}\sim B_p\sqrt R\cos(\theta/2),
\qquad B_p=-\sqrt{\frac2\gamma}\,n_{21}.
\tag{5}
\]

The normalized pressure used here has dimensional multiplier
\((d/a)\sqrt{2L}\,\delta_s^{-1/2}\); \(-\pi<\theta<\pi\) is measured around the right tip from the intact line.

Matching (3) to the outer opening now gives the independently derived answer

\[
\boxed{\frac{K(0^+)}{\sqrt{\pi L}}
        =\frac{n_{11}}{1+\gamma}.}
\tag{6}
\]

The parameter-free tip problem is independent of \(s\). Its constant leading resolvent intensity is therefore the initial positive-time intensity of the step response. Extraction at the drained tip precedes the large-\(s\), or small-time, limit.

### Numerical evaluation and an energy check

I computed (4) with an independently written **Wilson right spectral factorization**, \(N=F^*F\), rather than MP-37's projected inverse iteration. The lower analytic factor is normalized at physical infinity. On mapped frequency grids it gives

\[
n_{11}=1.07679067070,\quad n_{21}=0.41695945226,\quad
\frac{n_{11}}{4/3}=0.80759300302.
\]

See [wilson_factor.py](scratch/wilson_factor.py), [refinement output](scratch/wilson_refined.jsonl), and the finite-crack checks below.

There is also a conserved configurational integral for the **Laplace-domain** problem. Its density is the elastic-pressure enthalpy minus \(\kappa|\nabla p|^2/(2s)\). For example, using the equivalent strain-pressure variables, it is

\[
\mathcal F=
G e_{ij}e_{ij}+\frac{\lambda_d}{2}e_{kk}^2-\alpha pe_{kk}
-\frac{S_0}{2}p^2-\frac{\kappa}{2s}|\nabla p|^2,
\]

where \(\lambda_d=2G\nu/(1-2\nu)\) and
\(S_0=\alpha^2(1-2\nu)(1-2\nu_u)/(2Gd)\).
Its Euler equations are equilibrium and \(\nabla^2p=s\zeta/\kappa\). For \(p\sim A_p\sqrt r\cos(\theta/2)\), the inner contour contributes the hydraulic term \(\pi\kappa A_p^2/(4s)\). Thus

\[
\frac{1-\nu_u}{2G}K_{\rm out}^2
=\frac{1-\nu}{2G}K_{\rm tip}^2+\frac{\pi\kappa}{4s}A_p^2,
\]

or, with the above scales,

\[
\boxed{n_{11}^2+n_{21}^2=1+\gamma.}
\tag{7}
\]

About **13.0391%** of this outer configurational flux is in the hydraulic term. Ignoring that term would give a different modulus-based answer. This is a transformed-field identity; it is not an assertion that ordinary elastic energy release alone is conserved during the transient.

## 3. The leading volume correction, including the tips

Eliminating pressure from the constitutive equations yields

\[
2Ge_{ij}=\sigma_{ij}-\nu_uT\delta_{ij}
                 +\frac{2Gd}{a}\zeta\delta_{ij}.
\]

Apply elastic reciprocity with the unit-loaded undrained slit solution, whose trace is \(T^a\). The exact identity is

\[
\boxed{V_1-V_{1,0}=\frac da\int_\Omega\zeta T^a\,dA},
\qquad V_{1,0}=\frac{\pi(1-\nu_u)L^2}{G}.
\tag{8}
\]

One explicit auxiliary trace is
\(T^a=2\operatorname{Re}[z/\sqrt{z^2-L^2}]-2\), with \(z=x+iy\) and the branch chosen to give zero stress at infinity.

### Smooth face contribution

Away from the tips, \(T^a=-2\) at the faces, and the initial undrained pressure is \(2d/a\). Integrating the one-dimensional erfc drainage layer over both faces, using
\(\int_0^\infty\operatorname{erfc}(y/(2\sqrt{ct}))\,dy=2\sqrt{ct/\pi}\),
gives

\[
A_{\rm faces}=\frac{16d}{(1+\gamma)\sqrt\pi}
               =1.35405500051.
\tag{9}
\]

This is the paper's Eq. (18), and MP-39's answer, **as a face contribution**.

### Why a shrinking tip region still contributes at this order

Near either tip,

\[
T^a\sim\sqrt{\frac{2L}{r}}\cos(\theta/2).
\]

The outer undrained pressure has the same inverse-square-root singularity with the opposite sign. In a region of size \(\delta_s\), both the fluid-content response and the auxiliary trace have magnitude proportional to \(\delta_s^{-1/2}\). Their product integrated over area \(\delta_s^2\) is proportional to **\(\delta_s\)**, the same order as the face term in (8). Estimating only the area under a local opening correction misses this weighted coupling to the rest of the crack.

In the normalized tip cell let \(\mathcal H\) be its harmonic stress-pressure field, \(P\) its pressure, and \(Z=P+\mathcal H\). The equations are

\[
\nabla^2\mathcal H=0,\qquad(\nabla^2-1)Z=0,\qquad\nabla^2P=Z.
\]

Fluid content is proportional to \(Z\). Put
\(f=R^{-1/2}\cos(\theta/2)\). Green's identity gives

\[
I=\int_{\text{slit plane}}Zf\,dA=-\pi B_p
       =\pi\sqrt{\frac2\gamma}\,n_{21}.
\tag{10}
\]

Here both \(P\) and \(f\) vanish on the slit faces, the outer-circle term tends to zero by undrained matching, and the excluded inner circle contributes \(-\pi B_p\). The angular integral is over the **whole slit plane**, not just its upper half.

Restoring scales in (8), one tip contributes
\(dL\delta_s I/[G(1+\gamma)]\) to the unit resolvent volume. There are two tips. Inverting the step-response term uses
\(\mathcal L^{-1}(s^{-3/2})=2\sqrt{t/\pi}\). Therefore

\[
\boxed{A=
\frac{4d}{(1+\gamma)\sqrt\pi}
\left[4+\pi\sqrt{\frac2\gamma}\,n_{21}\right].}
\tag{11}
\]

Numerically,

\[
A_{\rm tips}=1.0861645789,\qquad
A=1.3540550005+1.0861645789=2.4402195794.
\]

The paper's poroelastic reciprocity relation between **loading-2 volume and loading-1 leak-off** is not disproved. The invalid step for loading-1 volume is replacing its entire spatially varying initial pressure response by a uniform face pressure-drop multiple of loading 2. Its singular tip pressure requires the extra calculation above.

### An explicit analytical check independent of the numerical factorization

At weak coupling, the leading harmonic tip field is \(f\). The pressure problem has the exact solution

\[
\mathcal H=f,\quad Z=e^{-R}f,\quad P=(e^{-R}-1)f.
\]

It satisfies the face condition, bounded inner pressure, and undrained outer matching. Hence \(B_p=-1\), and directly

\[
I=\int_0^\infty e^{-R}\,dR
       \int_{-\pi}^{\pi}\cos^2(\theta/2)\,d\theta=\pi.
\]

Equations (5) and (7) imply

\[
n_{21}=\sqrt{\gamma/2}+O(\gamma^{3/2}),\quad
n_{11}=1+\gamma/4+O(\gamma^2),
\]
\[
k_0=1-\frac34\gamma+O(\gamma^2),\qquad
A=\frac{16+4\pi}{\sqrt\pi}d+O(d^2)
\quad\text{at fixed }\nu.
\tag{12}
\]

The paper instead predicts slopes \(k_0=1-\gamma+O(\gamma^2)\) and
\(A/d\to16/\sqrt\pi\). The correct weak limit is **16.11684874**, rather than **9.02703334**. This is an explicit countercheck to the face-only argument, not a fit to the contested finite-coupling answers.

There is a direct physical-time check as well. The decoupled tip pressure is proportional to

\[
-r^{-1/2}\operatorname{erf}\!\left(\frac{r}{2\sqrt{ct}}\right)
 \cos(\theta/2).
\]

It solves the heat equation, has zero pressure at the tip and faces for positive time, and tends to the singular undrained initial pressure at fixed \(r\) as time vanishes. Its fluid-content correction is proportional to \(r^{-1/2}\operatorname{erfc}(r/(2\sqrt{ct}))\cos(\theta/2)\). Weighting by the auxiliary singularity gives the exact radial integral \(2\sqrt{ct/\pi}\) and angular integral \(\pi\). The additional \(\sqrt t\) contribution therefore also follows without numerical transform inversion. The physical-time heat-equation residual was checked symbolically in analytic_checks.py.

My independent factorization at \(\gamma=10^{-4}\) gives
\((1-k_0)/\gamma=0.74993191\) and \(A/d=16.11529080\), tending to (12). Zero coupling gives \(k_0=1\). See [weak_wilson.jsonl](scratch/weak_wilson.jsonl).

## 4. Step-by-step audit of the submitted derivations and code

### MP-37

1. **Bulk reduction:** its \(h,D_f,\beta\) equations and undrained elastic outer solution agree with the independent derivation above. Its positive full opening and volume normalization are correct.
2. **Surface symbol:** both equations agree with (1), including the pressure-gradient sign and the no-flux intact continuation. Its Airy representation is valid.
3. **Tip factor:** its matrix is (4). The ordering \(N=N_+N_-\) matters; factoring only the scalar determinant would not determine the needed first column. Its pole residue and opening amplitude agree with the matching calculation above.
4. **factor_check.py:** the circle map, exclusion of the constant and Nyquist modes from the lower projection, and the projected inverse equation are appropriate. If \(B\) is the analytic inverse found by the iteration, the physical-infinity normalization \(B(\infty)B(0)^{-1}\) is the required \(N_-(0)\). For these parameters the iteration converges in 42 iterations. The code does not reject exhaustion of its 200-iteration cap, so it is not a generally certified solver for arbitrary inputs; this did not occur here.
5. **solve.py:** the weights \((1/3,-2,8/3)\) on successive fourfold grids correctly cancel the first two powers of the zero-frequency \(N_{\rm grid}^{-1/2}\) cusp error. No contested amplitude is hard-coded. Regeneration of results.json is identical to the supplied JSON.
6. **Energy and volume derivations:** (7), (8), and the pressure Green identity reproduce its formulas. The face and two-tip factors, and the inverse-Laplace factor, are correct. The pressure coefficient is negative, consistent with the initial tensile-tip suction.
7. **finite_check.py:** eliminating pressure and symmetrizing with a scaled flux gives (2), up to the harmless opposite sign of the flux unknown. Restriction to the crack and its CG preconditioner are consistent. Subtraction of a same-grid undrained reference is necessary. I ran it on copied files and independently reproduced its finite-\(s\) volumes with my own operator. Its fitted limit still has periodic-box, mesh, and extrapolation errors.
8. **check_results.py:** reproduction, schema, energy-volume consistency, and weak-coupling checks pass. The energy relation is useful but not by itself an independent validation of the boundary operator; the finite-crack and symbolic checks provide that.

Fresh outputs: [mp37_reproduction.txt](scratch/mp37_reproduction.txt), [mp37_check_results.txt](scratch/mp37_check_results.txt), and [mp37_finite_check.txt](scratch/mp37_finite_check.txt).

### MP-38

1. **Constitutive reduction:** its drained Lamé modulus, storage coefficient, undrained modulus, and diffusivity follow from the stated two-index equations.
2. **Boundary conditions and Galerkin restriction:** it correctly uses crack-supported opening and flux, with pressure prescribed on the crack. Taking the hydraulic inverse after restriction is essential.
3. **Basis transforms:** the transforms of \(\sqrt{1-x^2}U_{2j}\) and \(T_{2j}/\sqrt{1-x^2}\) are its Bessel expressions, including the alternating signs absorbed into coefficients. The diagonal elastic matrix is \(\pi(2j+1)/2\).
4. **A, B, C and Schur complement:** converting (2) from full opening to upper-face displacement and rescaling flux gives its three integral matrices and right-hand side. Their signs and normalization are correct.
5. **Quadrature:** Simpson weights include the \(1/\pi\) Fourier factor. The special value of the first opening transform at zero is correct. Using direct Bessel evaluation at low argument and forward recurrence only at high argument avoids unstable low-argument recurrence.
6. **Fourier tails:** oscillatory averages give exactly its mechanical tail \(-d_{\rm code}\,\mathrm{outer}((2j+1)(-1)^j)/(\sqrt{R^2+q^2}+R)\) and hydraulic tail \(d_{\rm code}q\,\mathrm{asinh}(q/R)\,\mathrm{outer}((-1)^j)\). The leading mixed average vanishes, but finite-cutoff oscillatory errors remain, especially in tip extraction.
7. **Extraction and extrapolation:** volume is \(\pi u_0\); the intensity sum is divided by the drained \(1-\nu\). Its \(2/\sqrt\pi\) volume conversion is correct. The infinite-basis limit must precede \(q\to\infty\). Its exclusion of the last, insufficiently resolved intensity point is justified by its reported convergence tests. The polynomial fits are empirical, not error proofs.
8. **Fresh reproduction:** the copied production code regenerates exactly \(k_0=0.808,\ A=2.44\), at its reported precision. The independent factorization validates those rounded answers.

See [mp38_reproduction.txt](scratch/mp38_reproduction.txt). Its use of a different symbol named eta is a notation difference, not a discrepancy in the physics.

### MP-39

1. **Bulk limits and surface equations:** its two equations are algebraically equivalent to (1). Its mixed tip operator has the correct intact symmetry conditions.
2. **tip_check.py:** the periodic Toeplitz collocation system correctly discretizes these surface multipliers. Its same-grid elastic opening calibration uses the right drained/undrained conversion. I reproduced its three finest runs, including their volumes and residuals, in [mp39_tip_check.jsonl](scratch/mp39_tip_check.jsonl).
3. **Initial intensity:** solve.py fits recorded mesh-extrapolated values; it does not recompute those fields. The \(\ell=0.025\) correction was estimated rather than established by a second mesh. Its approximately 0.807 result is consequently provisional, as its answer acknowledges, but supports the correct side.
4. **Volume derivation:** the expansion giving a spatially constant interior displacement correction captures (9). The statement that its nonuniform endpoint expansion nevertheless gives the full leading integral is unjustified and contradicted by (8)–(12). Tip drainage has the singular auxiliary weight and creates an additional leading contribution.
5. **A check that failed:** its own reported successive finite-crack volume divided differences, which remove an unknown constant periodic-box initial volume, are **1.26970, 1.48907, 1.69308** per unit Laplace diffusion length. These do not establish its claimed limiting slope **1.2**; they move toward the full slope **2.1625883**. Finite box and mesh errors prevent using those few values alone for a precise limit.
6. **Reproduction:** its saved JSON is reproduced exactly, including the incorrect face-only volume coefficient. Reproducing arithmetic on stored data does not validate that coefficient.

See [mp39_reproduction.txt](scratch/mp39_reproduction.txt).

## 5. Independent numerical checks and checks that could fail

### Factorization refinement and independent algorithm

My Wilson code imports no submitted solver. It solves \(N=F^*F\) with \(F\) analytic in the lower half-plane, using the Wilson half-constant Fourier projection. The map is \(\xi=\rho\tan(\theta/2)\). Three grids, 16384, 65536, and 262144, and the two-term cusp extrapolation give:

| Frequency map scale \(\rho\) | \(k_0\) | \(A\) |
|---:|---:|---:|
| 0.01 | 0.807593003023 | 2.440219578445 |
| 0.1 | 0.807593003694 | 2.440219572522 |
| 1 | 0.807593018068 | 2.440219450435 |
| MP-37 independent inverse iteration, finer grids | 0.807593002749 | 2.440219579443 |

For the finest \(\rho=0.01\) grid, the full factor residual is \(1.17\times10^{-9}\) and unwanted positive Fourier coefficients have maximum \(2.0\times10^{-12}\). These are discretization floors: 8, 30, and 60 Wilson iterations give the same factor entries to roundoff. Extra iterations alone would not remove the cusp/grid error.

**Could fail:** wrong half-plane projection, reversed matrix-factor ordering, wrong normalization at infinity, or inadequate resolution could change the amplitudes. A small norm residual alone cannot distinguish an incorrectly oriented factor. This is why the independent finite-crack calculation and weak-coupling pressure solution are also used.

### Finite crack, no factorization input

[independent_fft.py](scratch/independent_fft.py) directly solves (2) on a finite slit, with FFT convolution on a periodic box and crack restriction. It uses no \(n_{11},n_{21}\), no candidate \(K\), and no candidate \(A\) in the solve. Both prescribed crack conditions and both intact symmetry conditions are enforced.

I checked box lengths \(32L,64L,128L\), and refined the \(64L\) box to 2097152 and 4194304 points. The following are first-order mesh extrapolations from those last two grids, before extrapolating the box:

| \(q=\sqrt{sL^2/c}\) | \(2q(\overline V-V_0)/\sqrt\pi\) | \(K/\sqrt{\pi L}\), from opening | From ahead-tip stress |
|---:|---:|---:|---:|
| 20 | 1.83244089 | 0.84801065 | 0.84801066 |
| 80 | 2.11607072 | 0.82759881 | 0.82759893 |
| 320 | 2.27232433 | 0.81752838 | 0.81753038 |
| 1280 | 2.35488913 | 0.81253754 | 0.81257025 |

Here \(\overline V=s\widehat V\) is the unit resolvent volume. Stress was obtained by applying the mechanical surface operator on the intact line. The same-grid undrained elastic solution calibrates the endpoint sampling bias; the known constant stress term is removed before forming the stress ratio. The near-tip sampling distance is refined at fixed \(q\) **before** extrapolating \(q\). At the last point, the finest \(q\,dx=0.01953\).

Quadratic/cubic extrapolations in \(q^{-1/2}\) give intensities around **0.80755–0.80767**, and volume limits around **2.438–2.441**, including fit and box sensitivity. Their precision is much lower than the factor calculation, but they independently reject 0.75 and 1.354. Box extrapolation of the finite-\(q\) volume agrees with MP-38's Bessel Galerkin solve to a few \(10^{-6}\) for the first three \(q\)'s and a few \(10^{-5}\) at the last one.

See [fine-grid output](scratch/independent_fft_stress_fine.jsonl), [box-128 output](scratch/independent_fft_stress_128.jsonl), and [summary.json](scratch/summary.json). The saved boundary arrays permit inspection without rerunning the solve.

**Could fail, and actually did at coarser resolution:** at large \(q\,dx\), the first opening samples underestimate the intensity and approach the undrained-opening/drained-modulus result. Stress samples outside the sufficiently small tip region give still different apparent intensities. Also, subtracting the infinite-domain \(V_0\) from a fixed periodic-grid volume gives a spurious growing short-time estimator; the log's raw estimators show this. I used a same-grid reference and varied the box. Tiny linear-system residuals alone do not establish accuracy.

## 6. Can the paper's numerical method resolve the smallest Fig. 9 times?

**No: its spatial approximation does not resolve this tip layer at the smallest plotted times.** This is not a criticism of adaptive kernel quadrature as such.

The PDF shows Fig. 9 beginning at \(\tau=10^{-4}\). Page 1655 specifies six equal quadratic elements over \(2L\), and Appendix B specifies stress collocation at the quarter and three-quarter positions of each element. Therefore:

* element length is \(L/3\);
* the nearest stress collocation point is \(L/12=0.08333L\) from a tip;
* the next density/pressure node after the endpoint is \(L/6=0.16667L\).

At this time, the physical diffusion length is \(0.01L\). The six Stehfest evaluations are \(s_nL^2/c=n\log2/\tau\), with layer lengths

\[
\delta_{s_n}/L=
0.01201,\ 0.008493,\ 0.006935,\ 0.006006,\ 0.005372,\ 0.004904.
\]

The nearest stress point is therefore **6.94–16.99 layer lengths** from the tip, and one element spans **28–68 layer lengths**. Factoring the elastic \(r^{-1/2}\) opening-slope singularity does not remove the variation of its regular amplitude across this layer.

Moreover, (5) implies a nonzero \(p_{,y}\sim r^{-1/2}\) face flux. Appendix B factors the opening slope but approximates the fluid-source density by ordinary piecewise quadratics; it does not give that density the appropriate endpoint singular factor. Refinement or a specifically enriched tip representation is needed.

Appendix B explicitly recognizes narrow Bessel-kernel peaks and improves their integration with adaptive Romberg quadrature. That fixes **kernel integration**, but it does not add degrees of freedom to the density approximation or supply its missing singular flux behavior.

A reproducible illustration is [coarse_basis_check.py](scratch/coarse_basis_check.py): retaining four smooth opening basis functions in the otherwise correct MP-38 method gives

| \(q\) | 80 | 320 | 1280 | 5120 | 20480 |
|---:|---:|---:|---:|---:|---:|
| Apparent normalized \(K\) | 0.78950 | 0.76257 | 0.75346 | 0.75090 | 0.75023 |

This is not a reproduction of the paper's IMSL implementation. It demonstrates the same failure of limit ordering: at fixed spatial basis, the opening converges to its undrained ellipse, and evaluating its endpoint coefficient with a drained modulus gives 0.75. Taking the resolved spatial limit first gives 0.80759.

The paper's **six-term Stehfest inversion itself** reproduces a constant step exactly. For \(s^{-3/2}\), its factor is 1.12911044 instead of \(2/\sqrt\pi=1.12837917\), a relative error of only **+0.064807%**. That cannot explain the contested coefficient or correct unresolved spatial fields. The Fig. 6 agreement at small time is not a volume-coefficient convergence test: the disputed increment is small compared with the initial volume 1.8849556.

## 7. What remains unproved, and reproduction

I have not reproduced the exact six-element IMSL program, which was not supplied. I cannot apportion the error in every plotted point among density interpolation, the unrepresented source singularity, quadrature tolerances, and transform inversion. The resolution assessment and independent continuum limits do not require that apportionment.

I have not established a certified numerical error bound, an elementary closed form for the finite-coupling factor entries, or a full direct time-stepping solution. The canonical factor problem (4), matching argument, and integrals (8)–(11) specify the continuum answer; the independent algorithms and limiting checks establish its numerical value to the stated estimated uncertainty. No steadily moving-crack formula was assumed.

Representative verification commands, run from this task folder:

~~~bash
OPENBLAS_NUM_THREADS=1 python3 scratch/analytic_checks.py
OPENBLAS_NUM_THREADS=1 python3 scratch/wilson_factor.py \
  --counts 16384 65536 262144 --scales 0.01 0.1 1
OPENBLAS_NUM_THREADS=1 python3 scratch/weak_wilson.py
OPENBLAS_NUM_THREADS=1 python3 scratch/independent_fft.py \
  --sizes 524288 1048576 --box 64 --qs 20 80 320 1280
~~~

The finer finite-crack command replaces the two sizes by 2097152 and 4194304. Submission reproductions and their copied code are retained under scratch/MP-37, scratch/MP-38, and scratch/MP-39. [PROCESS.md](PROCESS.md) records the UTC work log.
