import numpy as np
from scipy.special import psi
from dh_control import build, S_basis, dirichlet_coeffs
# 1. zeta via u-space code reproduces 50-digit results
for lam,ref in [(2.0,1.3684151e-5),(3.0,None)]:
    d=build('zeta',lam,14); ev=np.linalg.eigvalsh(d['A']/d['L'])
    print('zeta',lam,'minA',ev[0],'ref',ref)
# 2. archimedean u-space vs r-space for Gamma((s+1)/2), q=5
def arch_r(N,L,kinf,q,R=4000.,dr=0.002):
    r=np.arange(0,R,dr)+dr/2
    w=np.real(psi(0.25+kinf/2+0.5j*r))+np.log(q)
    S=np.array([S_basis(k,r,L) for k in range(1,N+1)])
    return (S*w)@S.T*dr*2/(2*np.pi) - np.log(np.pi)*L*np.eye(N)
L=np.log(3.0)
d=build('dh',3.0,8)
print('arch u vs r (DH):',np.abs(d['W']-arch_r(8,L,1,5.0)).max())
# 3. b(n) table
_,bz=dirichlet_coeffs('zeta',40); _,bd=dirichlet_coeffs('dh',40); _,bc=dirichlet_coeffs('chi5real',40)
print(' n  Lambda(n)   b_DH(n)   b_chi5(n)')
for n in range(2,31): print(f"{n:2d} {bz[n]:9.4f} {bd[n]:9.4f} {bc[n]:9.4f}")
