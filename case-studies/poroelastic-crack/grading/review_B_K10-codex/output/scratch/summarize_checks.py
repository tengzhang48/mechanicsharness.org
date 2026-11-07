import json
from pathlib import Path
import numpy as np

def lines(name):
    return [json.loads(s) for s in (Path('scratch')/name).read_text().splitlines()]

wilson=[r for r in lines('wilson_refined.jsonl') if 'quadratic_extrap_n11' in r]
finite=lines('independent_fft_stress_fine.jsonl')
rows=[r for r in finite if r.get('kind')=='result']
coarse=[r for r in rows if r['N']==2097152]
fine=[r for r in rows if r['N']==4194304]
res=[]
for old,new in zip(coarse,fine):
    q=new['q']
    out={'q':q,'A_dx0_box64':2*new['A_same_grid']-old['A_same_grid'],
         'K_opening_dx0_box64':2*new['K_samples']['1']-old['K_samples']['1'],
         'K_stress_dx0_box64':2*new['K_stress_samples']['1']-old['K_stress_samples']['1']}
    res.append(out)
x=np.array([r['q']**(-.5) for r in res])
fits=[]
for key in ['A_dx0_box64','K_opening_dx0_box64','K_stress_dx0_box64']:
    for degree in [1,2,3]:
        c=np.polynomial.polynomial.polyfit(x,[r[key] for r in res],degree)
        fits.append(dict(quantity=key,degree=degree,limit=c[0]))
    c=np.polynomial.polynomial.polyfit(x[:3],[r[key] for r in res[:3]],2)
    fits.append(dict(quantity=key,degree=2,points='first3',limit=c[0]))
boxes={}
for box in [32,64,128]:
    names={32:'independent_fft_stress_32.jsonl',64:'independent_fft_stress_64.jsonl',128:'independent_fft_stress_128.jsonl'}
    boxes[str(box)]=[r for r in lines(names[box]) if r.get('kind')=='dx_extrapolation']
box_ex=[]
for r64,r128 in zip(boxes['64'],boxes['128']):
    box_ex.append(dict(q=r64['q'],A_box_inf=(4*r128['A']-r64['A'])/3,
                       K_opening_box_inf=(4*r128['K_nearest']-r64['K_nearest'])/3))
repro=[]
for name in ['MP-37','MP-38','MP-39']:
    orig=json.loads((Path('gpt_answers')/name/'results.json').read_text())
    copy=json.loads((Path('scratch')/name/'results.json').read_text())
    repro.append(dict(run=name,identical=orig==copy,K=copy['K1_0'],A=copy['V1_short_coeff']))
summary=dict(wilson_extrapolations=wilson,finite_crack_fine_dx_extrapolation=res,
             finite_crack_empirical_time_fits=fits,finite_crack_box_extrapolations=box_ex,
             original_reproductions=repro)
Path('scratch/summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
