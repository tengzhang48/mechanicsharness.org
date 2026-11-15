"""Timing of the first route in the agents' container: NumPy and SciPy only, 2 CPUs. Task parameters.

Organizer only. Run on 13 Nov 2026 (23:59 UTC) as
  docker run --rm --network none --cpus 2 -v <this folder>:/key:ro magnet-bench-agent:20261111 python3 /key/container_timing.py
(the script was then mounted from another folder; the output is container_timing.log).
"""
import math, resource, sys, time
sys.path.insert(0, '/key')
import fe_bilayer as fe
import numpy, scipy
print('python', sys.version.split()[0], 'numpy', numpy.__version__, 'scipy', scipy.__version__, flush=True)
for refine in (0.25, 0.5, 1.0):
    t0 = time.time()
    oA = fe.solve(9.0, 1/3, 1.0, 1/3, 2.0, -1.0, 0.0, refine=refine)
    oB = fe.solve(9.0, 1/3, 1.0, 1/3, 2.0, 0.0, 1.0, refine=refine)
    print('refine %.2f ndof %d  psiA %.4f  psiB %.4f  both solves %.0f s  peak memory %.2f GB' %
          (refine, oA['ndof'], oA['psi_deg'], oB['psi_deg'], time.time() - t0, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6), flush=True)
