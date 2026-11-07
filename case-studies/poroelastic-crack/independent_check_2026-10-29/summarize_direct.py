"""Print the box, mesh and large-q sensitivity tables from saved direct solves.
python3 summarize_direct.py.  Recreate raw logs with finite_crack_check.py;
REPORT.md records the exact commands.  Does not read any factorization output.
"""
import json
from pathlib import Path
import numpy as np

face=16*.2/((4/3)*np.sqrt(np.pi))
records={}
for box,filename in [(16,'finite_crack_box16.jsonl'),(32,'finite_crack_results.jsonl'),(64,'finite_crack_box64.jsonl')]:
    rows=[json.loads(x) for x in Path(filename).read_text().splitlines()]
    records[box]=rows
    selected=[r for r in rows if r['kind']=='dx_extrapolation']
    assert len(selected)==5
    print(json.dumps(dict(kind='box_dx_table',box=box,rows=selected)))

# Fit infinite-domain values at every q BEFORE extrapolating q -> infinity.
qs=[r['q'] for r in records[32] if r['kind']=='dx_extrapolation']
aa=[];kk=[]
for j,q in enumerate(qs):
    vals=[[r for r in records[b] if r['kind']=='dx_extrapolation'][j] for b in (16,32,64)]
    a2=(4*vals[2]['A']-vals[1]['A'])/3
    k2=(4*vals[2]['K']-vals[1]['K'])/3
    # Three boxes remove both W^-2 and W^-4 terms, a sensitivity check.
    design=np.array([[1,b**-2,b**-4] for b in (16.,32.,64.)])
    a3=np.linalg.solve(design,np.array([v['A'] for v in vals]))[0]
    aa.append(a2);kk.append(k2)
    print(json.dumps(dict(kind='infinite_box',q=q,A=a2,K=k2,three_box_A=a3,box_model_difference=a3-a2)))
q=np.array(qs,dtype=float)
for degree in (1,2,3):
    for start in (0,1):
        design=np.array([q**(-j/2) for j in range(degree+1)]).T[start:]
        A=np.linalg.lstsq(design,np.array(aa)[start:],rcond=None)[0][0]
        K=np.linalg.lstsq(design,np.array(kk)[start:],rcond=None)[0][0]
        print(json.dumps(dict(kind='infinite_box_large_q',degree=degree,points=5-start,
                              A=A,tip=A-face,K=K)))
assert abs(K-.807593)<.0002,(K,'trusted intensity calibration failed')
