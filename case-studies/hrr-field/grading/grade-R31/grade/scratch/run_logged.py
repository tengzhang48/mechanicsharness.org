"""Run an unmodified grading input in scratch and record its PID/output."""
from datetime import datetime, timezone
from pathlib import Path
import os
import subprocess
import sys

script = Path(sys.argv[1]).resolve()
log = script.with_suffix('.run.log')
env = os.environ.copy()
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    env[name] = '1'
with log.open('w') as stream:
    process = subprocess.Popen(
        [sys.executable, '-u', str(script)], cwd=script.parent,
        stdout=stream, stderr=subprocess.STDOUT, env=env,
    )
    message = f"{datetime.now(timezone.utc):%H:%M:%S} — Started {script.name}, PID {process.pid}; log {log.relative_to(Path.cwd())}"
    print(message, flush=True)
    with (Path.cwd() / 'scratch/pids.log').open('a') as registry:
        registry.write(message + '\n')
    status = process.wait()
print(f"{script.name}: exit status {status}", flush=True)
print(log.read_text(), end='', flush=True)
sys.exit(status)
