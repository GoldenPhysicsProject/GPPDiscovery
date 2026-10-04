"""Collect shard JSONs into results/euler_search_summary.json and a markdown table."""
import json, glob, os, numpy as np
rows={}
for f in glob.glob('artifacts/**/*.json',recursive=True):
    d=json.load(open(f)); key=(d.get('mode','complex'),d['lam']); rows.setdefault(key,dict(primes=d['primes'],runs=[]))['runs']+=d['runs']
os.makedirs('discovery/odd_weil/results',exist_ok=True)
summ={}
lines=["| mode | lambda | starts | best margin | positive-margin starts (> 1e-9) | best point: max abs(alpha_p - 1) |","| --- | --- | --- | --- | --- | --- |"]
for key in sorted(rows):
    mode,lam=key; rs=rows[key]['runs']; m=len(rows[key]['primes'])
    best=max(rs,key=lambda r:r['margin']); x=np.array(best['x'])
    alpha=x[:m]*np.exp(1j*x[m:]) if mode=='complex' else x[:m].astype(complex); dev=float(np.max(np.abs(alpha-1)))
    npos=sum(r['margin']>1e-9 for r in rs)
    summ[f'{mode}_{lam}']=dict(mode=mode,lam=lam,starts=len(rs),best_margin=best['margin'],n_positive=npos,best_alpha=[[a.real,a.imag] for a in alpha],max_dev=dev,primes=rows[key]['primes'])
    lines.append(f"| {mode} | {lam} | {len(rs)} | {best['margin']:.3e} | {npos} | {dev:.4f} |")
json.dump(summ,open('discovery/odd_weil/results/euler_search_summary.json','w'),indent=1)
open('discovery/odd_weil/results/euler_search_summary.md','w').write("\n".join(lines)+"\n")
print("\n".join(lines))
