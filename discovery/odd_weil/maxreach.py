"""Try to break uniqueness: one-sided reach of F_lambda from Lambda along directions that
improve all near-null eigenvalues at first order (J d >= 0), plus axes and scalings."""
import mpmath as mp, pickle, sys, random
from mp_weil import Gsym, vm
mp.mp.dps=50
out=pickle.load(open('mp_N14.pkl','rb')); N=14
random.seed(1)
for lam in [float(x) for x in sys.argv[1:]]:
    d=out[lam]; L=mp.log(lam); A=mp.matrix(d['A'])
    E,Q=mp.eigsy(A)
    Ais=Q*mp.diag([1/mp.sqrt(e) for e in E])*Q.T
    ns=[n for n in range(2,int(lam**2)+1) if mp.log(n)<2*L]; m=len(ns)
    Es=[]
    for n in ns:
        Gn=mp.matrix(N,N)
        for j in range(1,N+1):
            for k in range(j,N+1):
                Gn[j-1,k-1]=Gn[k-1,j-1]=Gsym(j,k,mp.log(n),L)
        Es.append(Ais*(-2*Gn/mp.sqrt(n))*Ais)
    def reach(dv):
        nrm=mp.sqrt(sum(x*x for x in dv)); dv=[x/nrm for x in dv]
        M=sum((dv[i]*Es[i] for i in range(m)),mp.zeros(N,N)); ev,_=mp.eigsy((M+M.T)/2)
        return min([-1/e for e in ev if e<0],default=mp.inf)
    k=sum(1 for e in E if e/L<1e-8)
    J=mp.matrix(k,m)
    for i in range(k):
        for c in range(m):
            # q_i^T (A^{1/2} Es A^{1/2}) q_i  -> in whitened coords: e_i^T Q^T Es Q e_i * E_i
            qi=Q[:,i]; J[i,c]=(qi.T*Es[c]*qi)[0]*E[i]
    lam_vec=[vm(n) if vm(n) else mp.mpf(0) for n in ns]
    results=[]
    for c in range(m):
        for sgn in (1,-1):
            dv=[mp.mpf(0)]*m; dv[c]=sgn; results.append((reach(dv),f"axis {'+' if sgn>0 else '-'}b({ns[c]})"))
    results.append((reach(lam_vec),"scale Lambda up"))
    results.append((reach([-x for x in lam_vec]),"scale Lambda down"))
    # cone directions: minimize over random combos; maximize reach via random search + local refinement
    Jl=[[J[i,c] for c in range(m)] for i in range(k)]
    best=(mp.mpf(0),None)
    for trial in range(300):
        dv=[mp.mpf(random.gauss(0,1)) for _ in range(m)]
        # push into cone J d >= 0 by adding multiples of J rows where violated
        for it in range(50):
            viol=[i for i in range(k) if sum(Jl[i][c]*dv[c] for c in range(m))<0]
            if not viol: break
            for i in viol:
                s=sum(Jl[i][c]*dv[c] for c in range(m)); nn=sum(x*x for x in Jl[i])
                dv=[dv[c]-1.05*s/nn*Jl[i][c] for c in range(m)]
        r=reach(dv)
        if r>best[0]: best=(r,dv)
    results.append((best[0],"best random cone direction"))
    results.sort(key=lambda x:-x[0])
    print(f"\nlam={lam}: m={m} coefficients, k={k} near-null constraints")
    for r,lab in results[:6]: print(f"   reach {mp.nstr(r,4):>10}   {lab}")
    if best[1]:
        nb=mp.sqrt(sum(x*x for x in best[1])); print("   best cone dir (normalized):",[ (ns[c],round(float(best[1][c]/nb),3)) for c in range(m)])
