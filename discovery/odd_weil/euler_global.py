"""Global search: maximize the min eigenvalue of the odd semilocal form over degree-one Euler data
b(p^k) = log p * alpha_p^k (alpha_p real), primes p < lambda^2. Many random starts + Powell."""
import numpy as np, sys
from scipy.optimize import minimize
from dh_control import build, Gsym_matrix
rng=np.random.default_rng(7)
def primes_upto(M):
    return [p for p in range(2,M+1) if all(p%q for q in range(2,int(p**0.5)+1))]
for lam,N in [(3.0,16),(4.0,20),(5.0,22),(6.0,26)]:
    L=np.log(lam); W=build('zeta',lam,N)['W']; W=(W+W.T)/2
    M=int(np.floor(lam**2)); ps=[p for p in primes_upto(M) if np.log(p)<2*L]
    terms=[]  # (prime index, k, matrix 2/sqrt(n) log p G_n)
    for i,p in enumerate(ps):
        k=1
        while np.log(p**k)<2*L:
            n=p**k; terms.append((i,k,2*np.log(p)/np.sqrt(n)*Gsym_matrix(N,L,np.log(n)))); k+=1
    def margin(a):
        A=W.copy()
        for i,k,Mt in terms: A-=(a[i]**k)*Mt
        return np.linalg.eigvalsh((A+A.T)/2/L)[0]
    m1=margin(np.ones(len(ps)))
    best=(-np.inf,None)
    starts=[np.ones(len(ps))]+[1+0.3*rng.standard_normal(len(ps)) for _ in range(60)]
    for s in starts:
        r=minimize(lambda a:-margin(a),s,method='Powell',options={'maxiter':4000,'xtol':1e-6,'ftol':1e-12})
        if -r.fun>best[0]: best=(-r.fun,r.x)
    print(f"lambda={lam}: primes {ps}")
    print(f"   margin at alpha=1: {m1: .3e}    best margin found: {best[0]: .3e}")
    print(f"   best alpha: {np.round(best[1],4).tolist()}   max|alpha-1| = {np.max(np.abs(best[1]-1)):.4f}",flush=True)
