import mpmath as mp, pickle
from mp_weil import Gsym, vm
mp.mp.dps=50
out=pickle.load(open('mp_N14.pkl','rb')); N=14
track=[6,8,12,15,24]
print("max one-sided reach of a single coefficient b(n) (either sign), as lambda grows")
print("lam   "+"".join(f"   b({n})    " for n in track))
for lam in [2.5,3.0,4.0,5.0,6.0]:
    d=out[lam]; L=mp.log(lam); A=mp.matrix(d['A']); E,Q=mp.eigsy(A)
    Ais=Q*mp.diag([1/mp.sqrt(e) for e in E])*Q.T
    row=[]
    for n in track:
        if mp.log(n)>=2*L: row.append("   (absent) "); continue
        Gn=mp.matrix(N,N)
        for j in range(1,N+1):
            for k in range(j,N+1):
                Gn[j-1,k-1]=Gn[k-1,j-1]=Gsym(j,k,mp.log(n),L)
        M=Ais*(-2*Gn/mp.sqrt(n))*Ais; ev,_=mp.eigsy((M+M.T)/2)
        r=max(min([-1/e for e in ev if e<0],default=mp.inf),min([1/e for e in ev if e>0],default=mp.inf))
        row.append(f"  {mp.nstr(r,2):>9} ")
    print(f"{lam:4.1f} "+"".join(row))
