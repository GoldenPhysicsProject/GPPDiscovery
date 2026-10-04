"""One shard of the complex-local-root search (Actions matrix job).
Usage: python euler_complex_shard.py LAMBDA N SEED0 NSEEDS OUT.json [complex|real]"""
import sys, json, numpy as np
from scipy.optimize import minimize
from dh_control import build, Gsym_matrix
lam, N, s0, ns, out = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
mode = sys.argv[6] if len(sys.argv) > 6 else 'complex'
def primes_upto(M): return [p for p in range(2,M+1) if all(p%q for q in range(2,int(p**0.5)+1))]
L=np.log(lam); W=build('zeta',lam,N)['W']; W=(W+W.T)/2
ps=[p for p in primes_upto(int(lam**2)) if np.log(p)<2*L]; m=len(ps)
terms=[]
for i,p in enumerate(ps):
    k=1
    while np.log(p**k)<2*L:
        terms.append((i,k,2*np.log(p)/np.sqrt(p**k)*Gsym_matrix(N,L,np.log(p**k)))); k+=1
def margin(x):
    if mode=='real': r,th=x[:m],np.zeros(m)
    else: r,th=x[:m],x[m:]
    A=W.copy()
    for i,k,Mt in terms: A-=(r[i]**k*np.cos(k*th[i]))*Mt
    return np.linalg.eigvalsh((A+A.T)/2/L)[0]
res=[]
for seed in range(s0,s0+ns):
    rng=np.random.default_rng(seed)
    s=np.r_[np.abs(1+0.4*rng.standard_normal(m)),rng.uniform(-np.pi,np.pi,m)]
    if mode=='real': s=1+0.4*rng.standard_normal(m)
    rr=minimize(lambda x:-margin(x),s,method='Powell',options={'maxiter':8000,'xtol':1e-7,'ftol':1e-13})
    res.append(dict(seed=seed,margin=float(-rr.fun),x=rr.x.tolist()))
    print(seed, -rr.fun, flush=True)
json.dump(dict(lam=lam,N=N,mode=mode,primes=ps,runs=res),open(out,'w'))
