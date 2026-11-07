# REPRODUCE

From this folder, single command:

    python3 reproduce_all.py

It runs the Case A energy-difference FEM (fem_energy2.py / femA_H.py) and the
Case B half-domain FEM (femb2.py) at ell/H = 0.5, 1, 1.5, 2 and prints
gamma_A, gamma_B, KI_A/(sigma0 sqrt ell), KI_B/(sigma0 sqrt ell), matching
results.json. (Runtime is several minutes; each script prints incrementally.)

Note: compute.py contains an early, WRONG closed-form attempt and is NOT part
of the reported results.