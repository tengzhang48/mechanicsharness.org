# Task

Solve the problem below. Write your complete answer to `ANSWER.md` in this folder, with these
parts: the requested numbers or decisions; the model you chose and why; any assumption you
changed and why; how you verified the result; and what is missing or uncertain. You may write
and run code in this folder (`python3` with numpy and scipy). Budget: 45 minutes; write
`ANSWER.md` early, update it, and make sure it is complete by 40 minutes. Keep `PROCESS.md`:
one line per step, each starting with a time pasted from `date -u +%H:%M:%S` (never type a time).

Rules: use only this folder and your own knowledge. Do not open files outside this folder, and do
not use web search or web fetch. The session log will be checked.
Never stop or kill a process you did not start. Kill only by the PID you recorded when you started it; never use pkill, killall, or pgrep-based loops (this machine is shared with other runs).

Provided material: `./knowledge/` contains a mechanics topic (crack-tip and notch-tip asymptotic
fields), checks you can run on your own solution, and a knowledge entry for one paper. See
`./knowledge/README.md`. Use it as you see fit; you are free to choose your own method.

---

# Slip and strain contours near a stationary slit

Use mm, N, and MPa (= N/mm²). Angles are measured in radians unless a degree symbol is shown. A homogeneous body contains a straight, stationary slit along the negative x-axis with its tip at the origin. Write x = r cos θ, y = r sin θ, with material occupying −π < θ < π. All fields are independent of z, and the only displacement is u_z = w(x,y). Loading gives w(x,−y) = −w(x,y), with the upper face moving in the positive z direction. The slit faces are traction free. Loading is quasistatic and increases monotonically with a fixed loading pattern; there are no body forces.

The material follows the following small-strain constitutive law throughout the range in question. Define the **engineering** shear vector γ = (2ε_xz, 2ε_yz) = (w_,x, w_,y), the shear stress vector τ = (σ_xz, σ_yz), and τ = |τ| when it occurs as a scalar. Then

\[
\boldsymbol\gamma
 = \frac{\gamma_0}{\tau_0}
   \left(\frac{\tau}{\tau_0}\right)^{n-1}\boldsymbol\tau,
\qquad \tau_0=80\ {\rm MPa},\quad \gamma_0=0.002,\quad n=3.5.
\]

Take γ = 0 at τ = 0. This is the specified shear law, with no additional elastic term. Strain and engineering shear are dimensionless. The stored-work density is W = ∫ τ · dγ along the proportional loading path, in N/mm².

The applied loading supplies the following energy flux per unit length of the straight crack front:

\[
J = \int_\Gamma \left(W N_x - t_z w_{,x}\right)\,ds
  = 0.120\ {\rm N/mm},\qquad
t_z=\sigma_{xz}N_x+\sigma_{yz}N_y.
\]

Here Γ traverses a circular arc from the lower face (θ = −π) to the upper face (θ = +π), N is the radial outward unit normal, and ds is positive arc length. Work on the faces is zero. Geometry outside r = 0.30 mm introduces negligible corrections to the leading local field at all the positions and contour intersections requested below. Use that leading field; its displacement has no additive rigid translation. Ignore inertia, damage, and changes of slit shape in this stationary calculation.

Determine:

1. The relative face slip Δw(r) = w(r,+π) − w(r,−π), at r = 0.0200 mm and 0.0800 mm, in mm.
2. The **signed** stresses τ_rz and τ_θz at r = 0.0200 mm, θ = 60°. The positive basis vectors are e_r = (cos θ, sin θ) and e_θ = (−sin θ, cos θ). Report MPa.
3. The radii where |γ| = 0.0100 intersects the rays θ = 0°, 90°, and the upper face θ = 180°, in mm. The face value means the limit from within the material.

Give a short justification of your governing model and how the loading fixes its magnitude. Include one quantitative verification, such as a boundary residual, an energy-flux calculation, or a known limiting solution. Finally, state whether these results alone determine the load at which the slit starts to advance, and explain your answer in one sentence.

All seven requested numbers have a **2% relative tolerance**. A correct method earns full credit whether it uses numerical integration, an analytical solution, or another formulation. This is intended to be solvable in 45 minutes using Python with numpy and scipy.
