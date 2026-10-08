import numpy as np, cvxpy as cp
import joint_thread_check as J
def make(offset):
    J._cache.clear(); J.Nof=lambda lam,o=offset: int(round(16+10*np.log2(lam/3.0)))+o
S=[3.0,4.0,5.0,6.0,8.0]
res={}
for off in [0,8,16,24,32]:
    make(off)
    t,b,st=J.solve(S,cap=10)
    res[off]=b
    # cross-evaluate at larger bases
    cross=[]
    for off2 in [off+8,off+24]:
        make(off2); cross.append(min(J.margin_of(l,b) for l in S))
    make(off)
    print(f"offset {off:2d}: joint t*={t:.4e} ({st}) cross-margins at +8,+24 more basis: {['%.2e'%c for c in cross]}")
    print("    b(n<=12):",' '.join(f"{n}:{b[n]:.3f}" for n in range(2,13) if n in b))
print("Lambda     :",' '.join(f"{n}:{J.vm(n):.3f}" for n in range(2,13)))
