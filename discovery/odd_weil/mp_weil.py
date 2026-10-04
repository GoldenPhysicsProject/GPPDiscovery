"""High-precision odd semilocal Weil matrix, all in u-space.
A = W_inf - 2 sum c_n G(log n), with
W_inf(f) = (-log pi - gammaE) G(0) + int_0^inf 2[e^{-2u}G(0) - e^{-u/2}G(u)]/(1-e^{-2u}) du
(derived from psi(s) = -gamma + int_0^inf (e^{-t}-e^{-st})/(1-e^{-t}) dt, t=2u)."""
import mpmath as mp
mp.mp.dps = 50

def C(j,k,u,L):
    aj=j*mp.pi/L; ak=k*mp.pi/L
    lo=-L; hi=L-u
    def I(c,d):
        if d==0: return mp.cos(c)*(hi-lo)
        return (mp.sin(c+d*hi)-mp.sin(c+d*lo))/d
    return (I(aj*u,aj-ak)-I(aj*u,aj+ak))/2

def Gsym(j,k,u,L):
    if u>=2*L: return mp.mpf(0)
    return (C(j,k,u,L)+C(k,j,u,L))/2

def vm(n):
    for p in range(2,n+1):
        if n%p==0:
            m=n
            while m%p==0: m//=p
            return mp.log(p) if m==1 else 0
    return 0

def build(lam,N):
    L=mp.log(lam)
    G0=lambda j,k: L if j==k else mp.mpf(0)
    W=mp.matrix(N,N); Pr=mp.matrix(N,N); Lap=mp.matrix(N,N)
    terms=[]
    n=2
    while mp.log(n)<2*L:
        c=vm(n)
        if c: terms.append((n,c/mp.sqrt(n)))
        n+=1
    S=sum(c for _,c in terms)
    for j in range(1,N+1):
        for k in range(j,N+1):
            g0=G0(j,k)
            f=lambda u: 2*(mp.exp(-2*u)*g0-mp.exp(-u/2)*Gsym(j,k,u,L))/(-mp.expm1(-2*u))
            val=mp.quad(f,mp.linspace(0,2*L,9))
            val+= -g0*mp.log(-mp.expm1(-4*L))  # tail u>2L
            val+= (-mp.log(mp.pi)-mp.euler)*g0
            W[j-1,k-1]=W[k-1,j-1]=val
            p=sum(2*c*Gsym(j,k,mp.log(n),L) for n,c in terms)
            Pr[j-1,k-1]=Pr[k-1,j-1]=p
            l=sum(c*(2*g0-2*Gsym(j,k,mp.log(n),L)) for n,c in terms)
            Lap[j-1,k-1]=Lap[k-1,j-1]=l
    A=W-Pr
    Mass=W-2*S*mp.eye(N)*L
    return dict(A=A,W=W,Pr=Pr,Lap=Lap,Mass=Mass,S=S,L=L,terms=terms)
