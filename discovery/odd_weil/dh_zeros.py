import mpmath as mp, numpy as np, pickle
mp.mp.dps=30
kap=(mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
co=[1,kap,-kap,-1]
f=lambda s: sum(co[a-1]*mp.zeta(s,mp.mpf(a)/5) for a in range(1,5))*mp.power(5,-s)
Lam=lambda s: mp.power(5/mp.pi,s/2)*mp.gamma((s+1)/2)*f(s)
ts=np.arange(0.5,100,0.05); vals=[float(mp.re(Lam(mp.mpf(0.5)+1j*t))) for t in ts]
print('imag check', float(mp.im(Lam(mp.mpf(0.5)+20j))/abs(Lam(mp.mpf(0.5)+20j))))
zs=[]
for i in range(len(ts)-1):
    if vals[i]*vals[i+1]<0:
        zs.append(float(mp.findroot(lambda t: mp.re(Lam(0.5+1j*t)), (ts[i],ts[i+1]),solver='bisect')))
off=mp.findroot(f, mp.mpc(0.808517,85.699348))
print('on-line zeros <100:',len(zs)); print([round(z,4) for z in zs[:12]])
print('off-line zero:',off)
pickle.dump(dict(online=zs,off=[complex(off)]),open('dh_zeros.pkl','wb'))
