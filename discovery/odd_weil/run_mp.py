import mpmath as mp, numpy as np, pickle, sys
from mp_weil import build
N=int(sys.argv[1]); lams=[float(x) for x in sys.argv[2:]]
out={}
for lam in lams:
    d=build(lam,N)
    E,Q=mp.eigsy(d['A']/d['L'])
    Em,_=mp.eigsy(d['Mass']/d['L'])
    print(f"lam={lam} N={N} minA={mp.nstr(E[0],8)} 2nd={mp.nstr(E[1],6)} minMass={mp.nstr(Em[0],6)}",flush=True)
    out[lam]=dict(E=[float(x) for x in E],v=[float(Q[i,0]) for i in range(N)],
                  A=d['A'].tolist(),Lap=d['Lap'].tolist(),Mass=d['Mass'].tolist(),L=float(d['L']))
pickle.dump(out,open(f'mp_N{N}.pkl','wb'))
