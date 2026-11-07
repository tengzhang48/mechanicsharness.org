# Independent Finite-Crack Check

Only `TASK.md` and `factor_check.py` were read. No factorization routine or its
numerical results were used as solver inputs. The expected 2.4397 and face-only
1.354 are comparison labels only. Only numpy and scipy are used.

## Method

Set C=.625, b=.6675, beta=.4993757802746567, eta=1/3. Let
h=b*k/(m+k), D=2*beta*C*k*(m-k), lambda=2*beta*C*s/b,
so D=lambda*h. Eliminate p0 from Q=-m*p0+D*w. Write the unknown
crack-supported Q=sqrt(lambda)*z. The two equations on the crack are

    [ C*k*(1+eta*(m-k)/(m+k)) + lambda*h*h/m, -sqrt(lambda)*h/m ] [w] = [1]
    [ -sqrt(lambda)*h/m,                              1/m   ] [z]   [0]

The second equation enforces p0=0. Off-crack w=Q=0 are imposed exactly on
the collocation grid. This is the full mixed problem, not a face-only
approximation. Its Fourier symbol is symmetric positive definite except
for the elastic rigid mode at k=0; crack support removes that nullspace.
CG uses the restricted inverse full-space symbol as a preconditioner, with
the zero-mode regularization confined to the preconditioner. Unit traction
is used, not 1/s. Half-grid tip alignment avoids sampling the singular flux
at a tip.

Opening is integrated by the midpoint rule. For each grid, also solve the
undrained elastic problem C*(1+eta)*|k|*w=1. Define

    A_grid(s) = 2*sqrt(s/pi)*(V_grid(s)-V0_grid).

The infinite-domain elastic reference is V0=pi/[2*C*(1+eta)]=1.88495559215.
Same-grid subtraction cancels the dominant elastic sqrt-tip grid bias and
the s-independent periodic-box bias. The residual first-order grid error
is checked by doubling N and extrapolated as 2*A_fine-A_coarse. Neither
correction is derived from factor_check.

## Completed Reproduction

Command:

```bash
OPENBLAS_NUM_THREADS=1 python3 finite_check.py --box 64 --sizes 1048576 2097152 --s 100000 1000000 10000000
```

Exact pasted output:

```text
Independent full mixed-boundary FFT check; expected A~2.4397; face-only A~1.354
GRID box=64 N=1048576 dx=6.10351562e-05 crack_points=32768 V0_grid=1.88575682365 V0_infinite=1.88495559215
RESULT s=100000 V=1.8921127108 A_grid=2.267938310 A_raw=2.553837596 sqrt_s_dx=0.019301 iterations=10 relres=6.79e-11 info=0 seconds=2.74
RESULT s=1e+06 V=1.88782363585 A_grid=2.332147838 A_raw=3.236240762 sqrt_s_dx=0.0610352 iterations=9 relres=1.81e-10 info=0 seconds=2.52
RESULT s=1e+07 V=1.88641464428 A_grid=2.347267300 A_raw=5.206260155 sqrt_s_dx=0.19301 iterations=9 relres=1.53e-10 info=0 seconds=2.62
GRID box=64 N=2097152 dx=3.05175781e-05 crack_points=65536 V0_grid=1.88573494592 V0_infinite=1.88495559215
RESULT s=100000 V=1.89209563014 A_grid=2.269650024 A_raw=2.547742792 sqrt_s_dx=0.00965051 iterations=11 relres=6.95e-12 info=0 seconds=8.61
RESULT s=1e+06 V=1.88780694387 A_grid=2.337999326 A_raw=3.217405874 sqrt_s_dx=0.0305176 iterations=9 relres=1.29e-10 info=0 seconds=7.30
RESULT s=1e+07 V=1.88639815291 A_grid=2.366487192 A_raw=5.147414873 sqrt_s_dx=0.0965051 iterations=10 relres=7.26e-12 info=0 seconds=8.43
First-order dx extrapolation from last two grids (assumed, not guaranteed):
EXTRAP s=100000 A_dx0=2.271361738
EXTRAP s=1e+06 A_dx0=2.343850813
EXTRAP s=1e+07 A_dx0=2.385707083
FIT empirical powers s^(-1/4) degree=1 A_limit=2.438239292 max_fit_residual=0.000454
FIT empirical powers s^(-1/4) degree=2 A_limit=2.440642227 max_fit_residual=4.88e-15
Failure modes: finite periodic box; unresolved diffusion/tip layer when sqrt(s)*dx is large; midpoint sqrt-tip bias and singular flux; finite-s asymptotics; CG residual is not a discretization error bound. A_grid subtracts a same-grid undrained reference, not a factorization-derived correction.
```

## Box Convergence

The same spacings and s values, in a box half as long:

```bash
OPENBLAS_NUM_THREADS=1 python3 finite_check.py --box 32 --sizes 524288 1048576 --s 100000 1000000 10000000
```

Exact pasted output:

```text
Independent full mixed-boundary FFT check; expected A~2.4397; face-only A~1.354
GRID box=32 N=524288 dx=6.10351562e-05 crack_points=32768 V0_grid=1.8880352156 V0_infinite=1.88495559215
RESULT s=100000 V=1.89439785902 A_grid=2.270349114 A_raw=3.369235205 sqrt_s_dx=0.019301 iterations=10 relres=6.44e-11 info=0 seconds=1.45
RESULT s=1e+06 V=1.89010422635 A_grid=2.334628629 A_raw=5.809611566 sqrt_s_dx=0.0610352 iterations=10 relres=1.13e-11 info=0 seconds=1.72
RESULT s=1e+07 V=1.88869372753 A_grid=2.349734021 A_raw=13.338594934 sqrt_s_dx=0.19301 iterations=8 relres=1.95e-10 info=0 seconds=1.26
GRID box=32 N=1048576 dx=3.05175781e-05 crack_points=65536 V0_grid=1.88801328494 V0_infinite=1.88495559215
RESULT s=100000 V=1.89438073807 A_grid=2.272065336 A_raw=3.363126024 sqrt_s_dx=0.00965051 iterations=11 relres=7.25e-12 info=0 seconds=3.00
RESULT s=1e+06 V=1.89008749458 A_grid=2.340494947 A_raw=5.790731784 sqrt_s_dx=0.0305176 iterations=9 relres=1.18e-10 info=0 seconds=3.14
RESULT s=1e+07 V=1.8886771966 A_grid=2.369001602 A_raw=13.279608475 sqrt_s_dx=0.0965051 iterations=10 relres=7.25e-12 info=0 seconds=3.39
First-order dx extrapolation from last two grids (assumed, not guaranteed):
EXTRAP s=100000 A_dx0=2.273781559
EXTRAP s=1e+06 A_dx0=2.346361265
EXTRAP s=1e+07 A_dx0=2.388269183
FIT empirical powers s^(-1/4) degree=1 A_limit=2.440866973 max_fit_residual=0.000455
FIT empirical powers s^(-1/4) degree=2 A_limit=2.443271410 max_fit_residual=5.33e-15
Failure modes: finite periodic box; unresolved diffusion/tip layer when sqrt(s)*dx is large; midpoint sqrt-tip bias and singular flux; finite-s asymptotics; CG residual is not a discretization error bound. A_grid subtracts a same-grid undrained reference, not a factorization-derived correction.
```

An assumed O(box^-2) extrapolation (4*A_box64-A_box32)/3 gives
A_limit=2.437363398 (linear fit) and 2.439765833 (quadratic fit).
The empirical subleading variable s^-1/4 is not a proven expansion here;
the quadratic fit has only three points and therefore its tiny fit residual
is NOT an accuracy estimate.

## Additional Refinement And Failure

Command `OPENBLAS_NUM_THREADS=1 python3 finite_check.py --sizes 1048576 2097152 --s 100000 1000000 10000000 100000000`
was run before the extrapolation-printing code was added. Selected exact
output at s=1e8:

```text
RESULT s=1e+08 V=1.88822155942 A_grid=2.350125882 A_raw=36.852494257 sqrt_s_dx=0.305176 iterations=9 relres=1.05e-10 info=0 seconds=4.58
RESULT s=1e+08 V=1.88821333952 A_grid=2.381105081 A_raw=36.759742565 sqrt_s_dx=0.152588 iterations=10 relres=1.47e-11 info=0 seconds=7.48
```

The large raw coefficient is a failed use of the infinite-domain V0 in a
finite box; even small fixed volume biases grow like sqrt(s). The falling
high-s A_grid on coarse grids is another failure, caused by tip-layer
underresolution, not a physical asymptote.

Command `OPENBLAS_NUM_THREADS=1 python3 finite_check.py --box 16 --sizes 2097152 4194304 --s 1000000 10000000 100000000 1000000000`
exceeded the tool's 120-second timeout, after completing its first grid and
two points on the second grid. It does NOT constitute a completed check.
Selected exact output:

```text
RESULT s=1e+08 V=1.89741227222 A_grid=2.407073371 A_raw=140.558582829 sqrt_s_dx=0.0762939 iterations=10 relres=4.16e-11 info=0 seconds=13.65
RESULT s=1e+09 V=1.89726582344 A_grid=2.386178256 A_raw=439.259610334 sqrt_s_dx=0.241263 iterations=10 relres=3.15e-11 info=0 seconds=7.59
GRID box=16 N=4194304 dx=3.81469727e-06 crack_points=524288 V0_grid=1.8971961827 V0_infinite=1.88495559215
RESULT s=1e+06 V=1.89928396771 A_grid=2.355813108 A_raw=16.167840474 sqrt_s_dx=0.0038147 iterations=11 relres=7.36e-11 info=0 seconds=20.57
RESULT s=1e+07 V=1.89786774423 A_grid=2.396298250 A_raw=46.073763830 sqrt_s_dx=0.0120631 iterations=11 relres=6.37e-11 info=0 seconds=20.26
```

## Conclusion

Independent finite-crack numerics support A approximately 2.44, consistent
with 2.4397 and inconsistent with face-only 1.354. The difference between
empirical fit choices is about .0024; box and grid extrapolation introduce
additional uncertainty. This is a few-parts-in-a-thousand validation, not
a certified four-decimal result. No K extraction was attempted.

The parent-owned ANSWER.md, PROCESS.md and results.json were not edited.

## Complete Terminal Audit

This appendix contains the full available tool-returned output of EVERY
terminal check run, in submission order. It supersedes the selected-output
presentation above for audit completeness. There were eight terminal runs.
Runs 2/3, 4/5, and 6/7 were submitted in parallel pairs, so submission order
does not imply completion order. All ran in the task folder. No new numerical
checks were run to prepare this appendix.

No `date -u` commands (or other date commands) were actually run during this
subagent's work. There are therefore no date outputs or timestamped steps
to supply to the parent PROCESS. The `seconds` values below are solver
durations from `time.monotonic()`, not UTC timestamps.

All completed runs below have their full available output, with no known
truncation. Run 6 has its full available partial output and the tool timeout
metadata; subsequent numerical output is unavailable because the tool
terminated the command. No missing output or times have been reconstructed.

### Run 1: Initial Grid And Laplace Convergence

Compared N=8192, 16384, 32768, 65536 at box=32 and s=100 through 100000,
including the same-grid versus infinite-domain elastic subtraction.
Failure modes: finite-s corrections; unresolved diffusion/tip layer on
coarse grids; periodic-box error in V0; sqrt-tip midpoint bias. Small CG
residuals alone do not establish continuum accuracy.

Command:

```bash
OPENBLAS_NUM_THREADS=1 python3 finite_check.py
```

Full available pasted output:

```text
Independent full mixed-boundary FFT check; expected A~2.4397; face-only A~1.354
GRID box=32 N=8192 dx=0.00390625 crack_points=512 V0_grid=1.89079774309 V0_infinite=1.88495559215
RESULT s=100 V=2.03447002484 A_grid=1.621168096 A_raw=1.687089710 sqrt_s_dx=0.0390625 iterations=7 relres=1.6e-10 info=0 seconds=0.01
RESULT s=1000 V=1.94495917465 A_grid=1.932614324 A_raw=2.141076772 sqrt_s_dx=0.123526 iterations=7 relres=2.03e-11 info=0 seconds=0.01
RESULT s=10000 V=1.90932439079 A_grid=2.090508330 A_raw=2.749724471 sqrt_s_dx=0.390625 iterations=7 relres=1.01e-10 info=0 seconds=0.01
RESULT s=100000 V=1.89658462723 A_grid=2.064903909 A_raw=4.149528386 sqrt_s_dx=1.23526 iterations=7 relres=1.03e-11 info=0 seconds=0.01
GRID box=32 N=16384 dx=0.001953125 crack_points=1024 V0_grid=1.88939473674 V0_infinite=1.88495559215
RESULT s=100 V=2.03314041178 A_grid=1.621996251 A_raw=1.672086633 sqrt_s_dx=0.0195312 iterations=8 relres=9.1e-11 info=0 seconds=0.02
RESULT s=1000 V=1.94375028264 A_grid=1.939540806 A_raw=2.097940504 sqrt_s_dx=0.0617632 iterations=8 relres=5.23e-11 info=0 seconds=0.03
RESULT s=10000 V=1.90818690455 A_grid=2.120469066 A_raw=2.621372893 sqrt_s_dx=0.195312 iterations=7 relres=1.64e-10 info=0 seconds=0.03
RESULT s=100000 V=1.89546280895 A_grid=2.165238795 A_raw=3.749235777 sqrt_s_dx=0.617632 iterations=7 relres=1.82e-10 info=0 seconds=0.01
GRID box=32 N=32768 dx=0.0009765625 crack_points=2048 V0_grid=1.88869309257 V0_infinite=1.88495559215
RESULT s=100 V=2.0324750414 A_grid=1.622405557 A_raw=1.664578733 sqrt_s_dx=0.00976562 iterations=8 relres=7.43e-11 info=0 seconds=0.04
RESULT s=1000 V=1.94314504962 A_grid=1.942980995 A_raw=2.076344287 sqrt_s_dx=0.0308816 iterations=8 relres=2.87e-11 info=0 seconds=0.09
RESULT s=10000 V=1.90761820958 A_grid=2.135470776 A_raw=2.557202537 sqrt_s_dx=0.0976562 iterations=7 relres=5.46e-11 info=0 seconds=0.05
RESULT s=100000 V=1.89491174694 A_grid=2.218970246 A_raw=3.552603172 sqrt_s_dx=0.308816 iterations=8 relres=6.58e-11 info=0 seconds=0.09
GRID box=32 N=65536 dx=0.00048828125 crack_points=4096 V0_grid=1.88834223517 V0_infinite=1.88495559215
RESULT s=100 V=2.03214221363 A_grid=1.622608999 A_raw=1.660823174 sqrt_s_dx=0.00488281 iterations=10 relres=6.77e-11 info=0 seconds=0.17
RESULT s=1000 V=1.94284221886 A_grid=1.944694705 A_raw=2.065538534 sqrt_s_dx=0.0154408 iterations=10 relres=1.76e-12 info=0 seconds=0.27
RESULT s=10000 V=1.9073337061 A_grid=2.142958016 A_raw=2.525099758 sqrt_s_dx=0.0488281 iterations=9 relres=5.27e-13 info=0 seconds=0.26
RESULT s=100000 V=1.89463754099 A_grid=2.246321387 A_raw=3.454759681 sqrt_s_dx=0.154408 iterations=7 relres=5.47e-11 info=0 seconds=0.20
Failure modes: finite periodic box; unresolved diffusion/tip layer when sqrt(s)*dx is large; midpoint sqrt-tip bias and singular flux; finite-s asymptotics; CG residual is not a discretization error bound. A_grid subtracts a same-grid undrained reference, not a factorization-derived correction.
```

### Run 2: Larger Grids And Higher s

Compared N=131072, 262144, 524288 at box=32 and s=10000 through 10000000
to Run 1 and to each other. Failure modes: especially underresolution at
s=1e7, in addition to periodic-box and finite-s bias. The turnover of A_grid
on coarse grids is not evidence of a physical large-s limit.

Command:

```bash
OPENBLAS_NUM_THREADS=1 python3 finite_check.py --sizes 131072 262144 524288 --s 10000 100000 1000000 10000000
```

Full available pasted output:

```text
Independent full mixed-boundary FFT check; expected A~2.4397; face-only A~1.354
GRID box=32 N=131072 dx=0.000244140625 crack_points=8192 V0_grid=1.88816679762 V0_infinite=1.88495559215
RESULT s=10000 V=1.90719139292 A_grid=2.146695699 A_raw=2.509041435 sqrt_s_dx=0.0244141 iterations=10 relres=4.39e-12 info=0 seconds=0.26
RESULT s=100000 V=1.89450057516 A_grid=2.260049050 A_raw=3.405886874 sqrt_s_dx=0.077204 iterations=7 relres=1.1e-10 info=0 seconds=0.21
RESULT s=1e+06 V=1.89020472495 A_grid=2.299554734 A_raw=5.923012088 sqrt_s_dx=0.244141 iterations=8 relres=1.1e-11 info=0 seconds=0.24
RESULT s=1e+07 V=1.8887939617 A_grid=2.237877096 A_raw=13.696255339 sqrt_s_dx=0.77204 iterations=7 relres=1.73e-10 info=0 seconds=0.22
GRID box=32 N=262144 dx=0.000122070312 crack_points=16384 V0_grid=1.88807907664 V0_infinite=1.88495559215
RESULT s=10000 V=1.90712021807 A_grid=2.148562731 A_raw=2.501010214 sqrt_s_dx=0.012207 iterations=10 relres=1.83e-11 info=0 seconds=0.49
RESULT s=100000 V=1.89443209902 A_grid=2.266916088 A_raw=3.381452889 sqrt_s_dx=0.038602 iterations=9 relres=9.7e-11 info=0 seconds=0.59
RESULT s=1e+06 V=1.89013770024 A_grid=2.322907980 A_raw=5.847382806 sqrt_s_dx=0.12207 iterations=10 relres=1.81e-12 info=0 seconds=0.67
RESULT s=1e+07 V=1.88872689886 A_grid=2.311590478 A_raw=13.456958482 sqrt_s_dx=0.38602 iterations=8 relres=3.95e-11 info=0 seconds=0.59
GRID box=32 N=524288 dx=6.10351562e-05 crack_points=32768 V0_grid=1.8880352156 V0_infinite=1.88495559215
RESULT s=10000 V=1.90708462572 A_grid=2.149495753 A_raw=2.496994047 sqrt_s_dx=0.00610352 iterations=10 relres=9.15e-11 info=0 seconds=1.17
RESULT s=100000 V=1.89439785902 A_grid=2.270349114 A_raw=3.369235205 sqrt_s_dx=0.019301 iterations=10 relres=3.22e-11 info=0 seconds=1.43
RESULT s=1e+06 V=1.89010422635 A_grid=2.334628629 A_raw=5.809611566 sqrt_s_dx=0.0610352 iterations=10 relres=1.12e-11 info=0 seconds=1.44
RESULT s=1e+07 V=1.88869372752 A_grid=2.349734012 A_raw=13.338594924 sqrt_s_dx=0.19301 iterations=7 relres=1.89e-10 info=0 seconds=1.03
Failure modes: finite periodic box; unresolved diffusion/tip layer when sqrt(s)*dx is large; midpoint sqrt-tip bias and singular flux; finite-s asymptotics; CG residual is not a discretization error bound. A_grid subtracts a same-grid undrained reference, not a factorization-derived correction.
```

### Run 3: First Box Doubling

Compared box=64, N=262144 to box=32, N=131072 in Run 2, at the same
dx and s=10000, 100000, 1000000. Failure modes: this single-spacing box
comparison cannot independently establish either the continuum grid limit
or the box-error power; finite-s and tip-layer errors remain.

Command:

```bash
OPENBLAS_NUM_THREADS=1 python3 finite_check.py --box 64 --sizes 262144 --s 10000 100000 1000000
```

Full available pasted output:

```text
Independent full mixed-boundary FFT check; expected A~2.4397; face-only A~1.354
GRID box=64 N=262144 dx=0.000244140625 crack_points=8192 V0_grid=1.88588808809 V0_infinite=1.88495559215
RESULT s=10000 V=1.90489271706 A_grid=2.144442741 A_raw=2.249663639 sqrt_s_dx=0.0244141 iterations=10 relres=5.73e-12 info=0 seconds=0.98
RESULT s=100000 V=1.89221518515 A_grid=2.257665291 A_raw=2.590402987 sqrt_s_dx=0.077204 iterations=7 relres=1.19e-10 info=0 seconds=0.93
RESULT s=1e+06 V=1.88792389536 A_grid=2.297162514 A_raw=3.349371498 sqrt_s_dx=0.244141 iterations=9 relres=1.05e-11 info=0 seconds=0.74
Failure modes: finite periodic box; unresolved diffusion/tip layer when sqrt(s)*dx is large; midpoint sqrt-tip bias and singular flux; finite-s asymptotics; CG residual is not a discretization error bound. A_grid subtracts a same-grid undrained reference, not a factorization-derived correction.
```

### Run 4: Million-Point Refinement

Compared N=1048576 and 2097152 at box=32 and s=100000 through 100000000,
also to Run 2. Failure modes: s=1e8 is still tip-layer underresolved enough
to show appreciable grid sensitivity. A_raw diverges from the desired
coefficient because fixed periodic-volume bias is multiplied by sqrt(s).

Command:

```bash
OPENBLAS_NUM_THREADS=1 python3 finite_check.py --sizes 1048576 2097152 --s 100000 1000000 10000000 100000000
```

Full available pasted output:

```text
Independent full mixed-boundary FFT check; expected A~2.4397; face-only A~1.354
GRID box=32 N=1048576 dx=3.05175781e-05 crack_points=65536 V0_grid=1.88801328494 V0_infinite=1.88495559215
RESULT s=100000 V=1.89438073807 A_grid=2.272065336 A_raw=3.363126024 sqrt_s_dx=0.00965051 iterations=11 relres=7.25e-12 info=0 seconds=3.56
RESULT s=1e+06 V=1.89008749458 A_grid=2.340494947 A_raw=5.790731784 sqrt_s_dx=0.0305176 iterations=9 relres=1.18e-10 info=0 seconds=3.65
RESULT s=1e+07 V=1.8886771966 A_grid=2.369001602 A_raw=13.279608475 sqrt_s_dx=0.0965051 iterations=10 relres=7.25e-12 info=0 seconds=4.49
RESULT s=1e+08 V=1.88822155942 A_grid=2.350125882 A_raw=36.852494257 sqrt_s_dx=0.305176 iterations=9 relres=1.05e-10 info=0 seconds=4.58
GRID box=32 N=2097152 dx=1.52587891e-05 crack_points=131072 V0_grid=1.88800231957 V0_infinite=1.88495559215
RESULT s=100000 V=1.89437217729 A_grid=2.272923354 A_raw=3.360071327 sqrt_s_dx=0.00482525 iterations=11 relres=1.44e-11 info=0 seconds=8.56
RESULT s=1e+06 V=1.89007912941 A_grid=2.343428949 A_raw=5.781292697 sqrt_s_dx=0.0152588 iterations=11 relres=1.39e-11 info=0 seconds=7.52
RESULT s=1e+07 V=1.88866893955 A_grid=2.378665546 A_raw=13.250145277 sqrt_s_dx=0.0482525 iterations=9 relres=9.21e-11 info=0 seconds=6.26
RESULT s=1e+08 V=1.88821333952 A_grid=2.381105081 A_raw=36.759742565 sqrt_s_dx=0.152588 iterations=10 relres=1.47e-11 info=0 seconds=7.48
Failure modes: finite periodic box; unresolved diffusion/tip layer when sqrt(s)*dx is large; midpoint sqrt-tip bias and singular flux; finite-s asymptotics; CG residual is not a discretization error bound. A_grid subtracts a same-grid undrained reference, not a factorization-derived correction.
```

### Run 5: Box Halving

Compared box=16, N=524288 to box=32, N=1048576 in Run 4 at the same
spacing and s=100000, 1000000, 10000000. Failure modes: the smaller box
has a larger systematic error; three box sizes at differing spacings do
not by themselves prove an exact box-error expansion.

Command:

```bash
OPENBLAS_NUM_THREADS=1 python3 finite_check.py --box 16 --sizes 524288 --s 100000 1000000 10000000
```

Full available pasted output:

```text
Independent full mixed-boundary FFT check; expected A~2.4397; face-only A~1.354
GRID box=16 N=524288 dx=3.05175781e-05 crack_points=65536 V0_grid=1.89721556008 V0_infinite=1.88495559215
RESULT s=100000 V=1.90361049522 A_grid=2.281871601 A_raw=6.656532487 sqrt_s_dx=0.00965051 iterations=11 relres=8.45e-12 info=0 seconds=1.22
RESULT s=1e+06 V=1.89929874785 A_grid=2.350625680 A_raw=16.184518073 sqrt_s_dx=0.0305176 iterations=9 relres=1.03e-10 info=0 seconds=1.06
RESULT s=1e+07 V=1.89788233203 A_grid=2.379207836 A_raw=46.125816703 sqrt_s_dx=0.0965051 iterations=10 relres=8.21e-12 info=0 seconds=1.33
Failure modes: finite periodic box; unresolved diffusion/tip layer when sqrt(s)*dx is large; midpoint sqrt-tip bias and singular flux; finite-s asymptotics; CG residual is not a discretization error bound. A_grid subtracts a same-grid undrained reference, not a factorization-derived correction.
```

Between Runs 5 and 6, code was manually patched to print first-order dx
extrapolations and empirical large-s fits. Earlier runs did not print those
lines; their absence above is not missing output.

### Run 6: Timed-Out Finest Refinement

Attempted to compare N=2097152 and 4194304 at box=16 and s=1e6 through
1e9, also to Run 5. The first grid completed, but the second grid returned
only two s values. Failure modes: actual tool timeout made the intended
comparison incomplete; the first-grid s=1e9 value is underresolved;
small-box bias remains; no full-run extrapolation or fits were returned.

Command (tool timeout 120000 ms):

```bash
OPENBLAS_NUM_THREADS=1 python3 finite_check.py --box 16 --sizes 2097152 4194304 --s 1000000 10000000 100000000 1000000000
```

Full available pasted partial output, including tool metadata:

```text
Independent full mixed-boundary FFT check; expected A~2.4397; face-only A~1.354
GRID box=16 N=2097152 dx=7.62939453e-06 crack_points=262144 V0_grid=1.8971989509 V0_infinite=1.88495559215
RESULT s=1e+06 V=1.89928607911 A_grid=2.355071991 A_raw=16.170222936 sqrt_s_dx=0.00762939 iterations=11 relres=3.09e-11 info=0 seconds=8.45
RESULT s=1e+07 V=1.89786982758 A_grid=2.393854516 A_raw=46.081197724 sqrt_s_dx=0.0241263 iterations=11 relres=3.04e-11 info=0 seconds=8.27
RESULT s=1e+08 V=1.89741227222 A_grid=2.407073371 A_raw=140.558582829 sqrt_s_dx=0.0762939 iterations=10 relres=4.16e-11 info=0 seconds=13.65
RESULT s=1e+09 V=1.89726582344 A_grid=2.386178256 A_raw=439.259610334 sqrt_s_dx=0.241263 iterations=10 relres=3.15e-11 info=0 seconds=7.59
GRID box=16 N=4194304 dx=3.81469727e-06 crack_points=524288 V0_grid=1.8971961827 V0_infinite=1.88495559215
RESULT s=1e+06 V=1.89928396771 A_grid=2.355813108 A_raw=16.167840474 sqrt_s_dx=0.0038147 iterations=11 relres=7.36e-11 info=0 seconds=20.57
RESULT s=1e+07 V=1.89786774423 A_grid=2.396298250 A_raw=46.073763830 sqrt_s_dx=0.0120631 iterations=11 relres=6.37e-11 info=0 seconds=20.26


<shell_metadata>
shell tool terminated command after exceeding timeout 120000 ms. If this command is expected to take longer and is not waiting for interactive input, retry with a larger timeout value in milliseconds.
</shell_metadata>
```

### Run 7: Completed Box-64 Extrapolation

Compared N=1048576 and 2097152 at box=64 and s=1e5, 1e6, 1e7; printed
dx extrapolation and two empirical fits. Compared the fitted limit to
expected 2.4397 and face-only 1.354. Failure modes: first-order dx behavior
is assumed, not guaranteed; empirical s^-1/4 fit may misrepresent the true
subleading expansion; quadratic fit interpolates only three points and its
near-zero residual is not a confidence bound; periodic-box bias remains.

Command:

```bash
OPENBLAS_NUM_THREADS=1 python3 finite_check.py --box 64 --sizes 1048576 2097152 --s 100000 1000000 10000000
```

Full available pasted output:

```text
Independent full mixed-boundary FFT check; expected A~2.4397; face-only A~1.354
GRID box=64 N=1048576 dx=6.10351562e-05 crack_points=32768 V0_grid=1.88575682365 V0_infinite=1.88495559215
RESULT s=100000 V=1.8921127108 A_grid=2.267938310 A_raw=2.553837596 sqrt_s_dx=0.019301 iterations=10 relres=6.79e-11 info=0 seconds=2.74
RESULT s=1e+06 V=1.88782363585 A_grid=2.332147838 A_raw=3.236240762 sqrt_s_dx=0.0610352 iterations=9 relres=1.81e-10 info=0 seconds=2.52
RESULT s=1e+07 V=1.88641464428 A_grid=2.347267300 A_raw=5.206260155 sqrt_s_dx=0.19301 iterations=9 relres=1.53e-10 info=0 seconds=2.62
GRID box=64 N=2097152 dx=3.05175781e-05 crack_points=65536 V0_grid=1.88573494592 V0_infinite=1.88495559215
RESULT s=100000 V=1.89209563014 A_grid=2.269650024 A_raw=2.547742792 sqrt_s_dx=0.00965051 iterations=11 relres=6.95e-12 info=0 seconds=8.61
RESULT s=1e+06 V=1.88780694387 A_grid=2.337999326 A_raw=3.217405874 sqrt_s_dx=0.0305176 iterations=9 relres=1.29e-10 info=0 seconds=7.30
RESULT s=1e+07 V=1.88639815291 A_grid=2.366487192 A_raw=5.147414873 sqrt_s_dx=0.0965051 iterations=10 relres=7.26e-12 info=0 seconds=8.43
First-order dx extrapolation from last two grids (assumed, not guaranteed):
EXTRAP s=100000 A_dx0=2.271361738
EXTRAP s=1e+06 A_dx0=2.343850813
EXTRAP s=1e+07 A_dx0=2.385707083
FIT empirical powers s^(-1/4) degree=1 A_limit=2.438239292 max_fit_residual=0.000454
FIT empirical powers s^(-1/4) degree=2 A_limit=2.440642227 max_fit_residual=4.88e-15
Failure modes: finite periodic box; unresolved diffusion/tip layer when sqrt(s)*dx is large; midpoint sqrt-tip bias and singular flux; finite-s asymptotics; CG residual is not a discretization error bound. A_grid subtracts a same-grid undrained reference, not a factorization-derived correction.
```

### Run 8: Matched-Spacing Box-32 Extrapolation

Compared box=32, N=524288/1048576 to Run 7 at matched spacings and s
values. Compared first-order dx extrapolations and fitted limits across
boxes. Failure modes: same fit and grid assumptions as Run 7; assuming
O(box^-2) with just these two fitted limits is not an independently
certified error model. Derived box-extrapolated numbers in the preceding
sections were arithmetic from Runs 7/8, not another terminal check; they
have no separate terminal output.

Command:

```bash
OPENBLAS_NUM_THREADS=1 python3 finite_check.py --box 32 --sizes 524288 1048576 --s 100000 1000000 10000000
```

Full available pasted output:

```text
Independent full mixed-boundary FFT check; expected A~2.4397; face-only A~1.354
GRID box=32 N=524288 dx=6.10351562e-05 crack_points=32768 V0_grid=1.8880352156 V0_infinite=1.88495559215
RESULT s=100000 V=1.89439785902 A_grid=2.270349114 A_raw=3.369235205 sqrt_s_dx=0.019301 iterations=10 relres=6.44e-11 info=0 seconds=1.45
RESULT s=1e+06 V=1.89010422635 A_grid=2.334628629 A_raw=5.809611566 sqrt_s_dx=0.0610352 iterations=10 relres=1.13e-11 info=0 seconds=1.72
RESULT s=1e+07 V=1.88869372753 A_grid=2.349734021 A_raw=13.338594934 sqrt_s_dx=0.19301 iterations=8 relres=1.95e-10 info=0 seconds=1.26
GRID box=32 N=1048576 dx=3.05175781e-05 crack_points=65536 V0_grid=1.88801328494 V0_infinite=1.88495559215
RESULT s=100000 V=1.89438073807 A_grid=2.272065336 A_raw=3.363126024 sqrt_s_dx=0.00965051 iterations=11 relres=7.25e-12 info=0 seconds=3.00
RESULT s=1e+06 V=1.89008749458 A_grid=2.340494947 A_raw=5.790731784 sqrt_s_dx=0.0305176 iterations=9 relres=1.18e-10 info=0 seconds=3.14
RESULT s=1e+07 V=1.8886771966 A_grid=2.369001602 A_raw=13.279608475 sqrt_s_dx=0.0965051 iterations=10 relres=7.25e-12 info=0 seconds=3.39
First-order dx extrapolation from last two grids (assumed, not guaranteed):
EXTRAP s=100000 A_dx0=2.273781559
EXTRAP s=1e+06 A_dx0=2.346361265
EXTRAP s=1e+07 A_dx0=2.388269183
FIT empirical powers s^(-1/4) degree=1 A_limit=2.440866973 max_fit_residual=0.000455
FIT empirical powers s^(-1/4) degree=2 A_limit=2.443271410 max_fit_residual=5.33e-15
Failure modes: finite periodic box; unresolved diffusion/tip layer when sqrt(s)*dx is large; midpoint sqrt-tip bias and singular flux; finite-s asymptotics; CG residual is not a discretization error bound. A_grid subtracts a same-grid undrained reference, not a factorization-derived correction.
```
