# MP-41 process log

- 21:07:08  read TASK.md; identified problem: Biot poroelastic suddenly-loaded /
            pressurized mode-I crack (Rice-Cleary type). Confirmed python3 +
            numpy 2.4.3 / scipy 1.18.0 available.
- 21:07     derived governing transport law (T); established Beltrami harmonic
            combination T2 = sigma_kk + [alpha(1-2nu)/(1-nu)] p; established that K
            is controlled by the total face traction.
- 21:08:06  wrote solve.py; ran it. All quantities computed; numerical
            boundary-layer check vs analytic integral agrees to 7.8e-5.
            Earlier explicit FTCS check rejected (stability dt<=dy^2/2c);
            replaced by stiff BDF (method of lines) after a timed-out first run
            (reduced ny/tmax).
- 21:08     wrote results.json, ANSWER.md; final read-back verification.
