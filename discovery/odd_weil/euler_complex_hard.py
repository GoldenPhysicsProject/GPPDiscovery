"""Hard multi-start search with complex local roots at lambda = 5, 6 (2 worker processes)."""
import numpy as np, sys, json
from multiprocessing import Pool
from scipy.optimize import minimize
from dh_control import build, Gsym_matrix
def primes_upto(M): return [p for p in range(2,M+1) if all(p%q for q in range(2,int(p**0.5)+1))]
SETUP={}
def setup(lam,N):
    L=np.log(lam); W=build('zeta',lam,N)['W']; W=(W+W.T)/2
    ps=[p for p in primes_upto(int(lam**2)) if np.log(p)<2*L]
    terms=[]
    for i,p in enumerate(ps):
        k=1
        while np.log(p**k)<2*L:
            terms.append((i,k,2*np.log(p)/np.sqrt(p**k)*Gsym_matrix(N,L,np.log(p**k)))); k+=1
    return L,W,ps,terms
def run(args):
    lam,N,seed=args
    if lam not in SETUP: SETUP[lam]=setup(lam,N)
    L,W,ps,terms=SETUP[lam]; m=len(ps); rng=np.random.default_rng(seed)
    def margin(x):
        r,th=x[:m],x[m:]; A=W.copy()
        for i,k,Mt in terms: A-=(r[i]**k*np.cos(k*th[i]))*Mt
        return np.linalg.eigvalsh((A+A.T)/2/L)[0]
    s=np.r_[np.abs(1+0.4*rng.standard_normal(m)),rng.uniform(-np.pi,np.pi,m)]
    rr=minimize(lambda x:-margin(x),s,method='Powell',options={'maxiter':8000,'xtol':1e-7,'ftol':1e-13})
    return lam,-rr.fun,rr.x.tolist()
if __name__=='__main__':
    jobs=[(lam,N,s) for lam,N in [(5.0,22),(6.0,26)] for s in range(int(sys.argv[1]))]
    out=[]
    with Pool(2) as pool:
        for lam,mg,x in pool.imap_unordered(run,jobs):
            out.append(dict(lam=lam,margin=mg,x=x))
            json.dump(out,open('euler_complex_hard.json','w'))
            print(lam,f"{mg:.3e}",flush=True)
