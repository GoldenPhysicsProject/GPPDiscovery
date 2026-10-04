"""Complex local roots alpha_p = r_p e^{i theta_p}: b(p^k) = log p * r_p^k cos(k theta_p).
Global search for positive-margin points at the zeta archimedean place."""
import numpy as np
from scipy.optimize import minimize
from dh_control import build, Gsym_matrix
rng=np.random.default_rng(11)
def primes_upto(M): return [p for p in range(2,M+1) if all(p%q for q in range(2,int(p**0.5)+1))]
for lam,N in [(4.0,20),(5.0,22)]:
    L=np.log(lam); W=build('zeta',lam,N)['W']; W=(W+W.T)/2
    ps=[p for p in primes_upto(int(lam**2)) if np.log(p)<2*L]; m=len(ps)
    terms=[]
    for i,p in enumerate(ps):
        k=1
        while np.log(p**k)<2*L:
            terms.append((i,k,2*np.log(p)/np.sqrt(p**k)*Gsym_matrix(N,L,np.log(p**k)))); k+=1
    def margin(x):
        r,th=x[:m],x[m:]; A=W.copy()
        for i,k,Mt in terms: A-=(r[i]**k*np.cos(k*th[i]))*Mt
        return np.linalg.eigvalsh((A+A.T)/2/L)[0]
    best=(-np.inf,None)
    starts=[np.r_[np.ones(m),np.zeros(m)]]+[np.r_[1+0.2*rng.standard_normal(m),rng.uniform(-np.pi,np.pi,m)] for _ in range(40)]
    for s in starts:
        rr=minimize(lambda x:-margin(x),s,method='Powell',options={'maxiter':6000,'xtol':1e-6,'ftol':1e-12})
        if -rr.fun>best[0]: best=(-rr.fun,rr.x)
    r,th=best[1][:m],np.angle(np.exp(1j*best[1][m:]))
    print(f"lambda={lam}: best margin {best[0]: .3e};  r = {np.round(r,4).tolist()}  theta = {np.round(th,4).tolist()}",flush=True)
