"""Collect shard JSONs into results/euler_complex_summary.json and a markdown table."""
import json, glob, os, numpy as np
rows={}
for f in glob.glob('artifacts/**/*.json',recursive=True):
    d=json.load(open(f)); rows.setdefault(d['lam'],dict(primes=d['primes'],runs=[]))['runs']+=d['runs']
os.makedirs('discovery/odd_weil/results',exist_ok=True)
summ={}
lines=["| lambda | starts | best margin | positive-margin starts (> 1e-9) | best point: max abs(alpha_p - 1) |","| --- | --- | --- | --- | --- |"]
for lam in sorted(rows):
    rs=rows[lam]['runs']; m=len(rows[lam]['primes'])
    best=max(rs,key=lambda r:r['margin']); x=np.array(best['x'])
    alpha=x[:m]*np.exp(1j*x[m:]); dev=float(np.max(np.abs(alpha-1)))
    npos=sum(r['margin']>1e-9 for r in rs)
    summ[lam]=dict(starts=len(rs),best_margin=best['margin'],n_positive=npos,best_alpha=[[a.real,a.imag] for a in alpha],max_dev=dev,primes=rows[lam]['primes'])
    lines.append(f"| {lam} | {len(rs)} | {best['margin']:.3e} | {npos} | {dev:.4f} |")
json.dump(summ,open('discovery/odd_weil/results/euler_complex_summary.json','w'),indent=1)
open('discovery/odd_weil/results/euler_complex_summary.md','w').write("\n".join(lines)+"\n")
print("\n".join(lines))
