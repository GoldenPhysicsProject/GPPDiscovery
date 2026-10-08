"""Truncation-artifact test for the thread-extension claim: does the joint far-from-Lambda margin survive basis refinement?
Galerkin compression: lambda_min over a larger N-dim space can only decrease, so a margin that is positive at N(lam) may
vanish at larger N. Fix the joint optimum b* found at N(lam)+0, then (a) evaluate its margin at N(lam)+k, (b) re-optimise at +k."""
import numpy as np, cvxpy as cp, sys
import joint_thread_check as J
from dh_control import build, Gsym_matrix
def make(offset):
    J._cache.clear(); J.Nof=lambda lam,o=offset: int(round(16+10*np.log2(lam/3.0)))+o
S=[3.0,4.0,5.0,6.0,8.0,10.0,12.0]
make(0)
t0,b0,_=J.solve(S,cap=10)
print(f"base (offset 0): joint t*={t0:.3e}")
for off in [0,4,8,12,16]:
    make(off)
    ms=[J.margin_of(l,b0) for l in S]
    print(f"  frozen b*, basis offset +{off}: margins min={min(ms):.3e} max={max(ms):.3e}  per-lam={['%.2e'%m for m in ms]}")
for off in [4,8,12]:
    make(off)
    t,b,st=J.solve(S,cap=10)
    print(f"  re-optimised at offset +{off}: joint t*={'None' if t is None else '%.3e'%t} ({st}); dist(n<9)={None if b is None else round(J.dist(b,9),3)}")
