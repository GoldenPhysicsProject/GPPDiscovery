"""Two-box energy E(t) = W((f-T_t f)*(f-T_t f)~) = 2(C(0) - C(t)) for zeta, l = log 2
(Daniel's fixed-window paper, Thm onesided / eq. twobox). Prime side:
C(t) = A e^{t/2} - S(t) - D(t), D(t) = sum_k H(2k+1/2) e^{-(2k+1/2)t};
C(0) from the paper's eq. Czero; cross-check C(0) = sum_gamma H(i gamma) with zeros + tail.
Evidence only."""
import json, math, os, sys
import numpy as np, mpmath as mp
sys.path.insert(0, os.path.dirname(__file__))
from fixed_window import mangoldt, window

mp.mp.dps = 30
L = mp.log(2)
H = lambda z: L * (mp.sinh(L * z / 2) / (L * z / 2)) ** 2
A = H(mp.mpf(1) / 2)
C0 = 2 * A - mp.log(4 * mp.pi) - mp.euler + mp.log(3) - mp.quad(lambda u: (mp.exp(u / 2) * (1 - u / L) - 1) / mp.sinh(u), [0, L])
# zero side: 2 * sum_{gamma>0} H(i gamma) = 2 * sum l (sin(l g/2)/(l g/2))^2
K = int(sys.argv[1]) if len(sys.argv) > 1 else 400
gs = [mp.im(mp.zetazero(k)) for k in range(1, K + 1)]
zs = 2 * sum(L * (mp.sin(L * g / 2) / (L * g / 2)) ** 2 for g in gs)
T = gs[-1]
# tail: density log(g/2pi)/(2pi); average sin^2 = 1/2
tail = 2 * mp.quad(lambda g: L * (mp.mpf(1) / 2) / (L * g / 2) ** 2 * mp.log(g / (2 * mp.pi)) / (2 * mp.pi), [T, mp.inf])
N = int(2e7)
lam = mangoldt(N)
xs = np.arange(2, int(N / 2) - 1, dtype=np.float64)
S = window(lam, math.log(2), xs)
t = np.log(xs)
D = sum(float(H(2 * k + mp.mpf(1) / 2)) * np.exp(-(2 * k + 0.5) * t) for k in range(1, 40))
C = float(A) * np.sqrt(xs) - S - D
E = 2 * (float(C0) - C)
i = int(np.argmin(E))
blocks = []; k = 2
while 2 * k <= len(xs):
    seg = E[k - 2:2 * k - 2]; blocks.append([k, float(seg.min()), float(seg.mean()), float(seg.max())]); k *= 2
res = {"A": float(A), "C0_exact": float(C0), "C0_zero_side": float(zs + tail), "zeros_used": K, "tail": float(tail),
       "E_min": [float(xs[i]), float(E[i])], "E_over_2C0_min": float(E[i] / (2 * C0)),
       "E_mean_over_2C0": float(E.mean() / (2 * C0)), "C_max": float(C.max()), "dyadic_min_mean_max": blocks}
os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
json.dump(res, open(os.path.join(os.path.dirname(__file__), "results", "two_box.json"), "w"), indent=1)
print(json.dumps(res)[:2000])
