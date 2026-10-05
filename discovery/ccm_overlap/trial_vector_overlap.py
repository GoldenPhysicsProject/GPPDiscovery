"""Trial-vector overlap test (go/no-go for Layer C of the CCM closure; spec: CLAUDE_CODE_TASK_trial_vector_overlap.md).

For (lambda, N): build the truncated Weil matrix QW (ccm.build), T = QW - eps0 I, eigendecompose,
build the Poisson-prolate trial vector k_lambda = E(h_lambda) (CCM 2511.22755 eq. 7.6: E(h)(u) = u^{1/2} sum_{n>=1} h(nu),
h_lambda = the combination of prolate h_{0,lambda}, h_{4,lambda} with vanishing integral), Galerkin-project onto V_n
(c_n = <V_n, k_lambda>, V_n(u) = L^{-1/2} exp(2 pi i n log(lambda u)/L), L = 2 log lambda), normalize eta^T k = 1,
and record overlaps with the near-null eigenvectors of T.  Numerical evidence only.
"""
import argparse, json, math, os, sys, time
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ccm import build
from prolate import make_hlam, eval_H

ZS = [mp.mpc(1, 1), mp.mpc("0.3", -2), mp.mpc(5, "0.1")]


def parse_lam(s):
    return mp.sqrt(mp.mpf(s[5:-1])) if s.startswith("sqrt(") else mp.mpf(s)


def trial_vector(lam, N, dps, M_deg=7):
    """Galerkin coefficients c_n = <V_n, k_lambda>, |n| <= N (list of mpc), plus info."""
    lam = mp.mpf(lam); L = 2 * mp.log(lam)
    g, K, info = make_hlam(lam, dps)
    Msum = int(mp.floor(lam ** 2 + mp.mpf("1e-30")))
    gl = mp.calculus.quadrature.GaussLegendre(mp.mp)
    nodes = gl.calc_nodes(M_deg, mp.mp.prec)           # on [-1,1]
    # breakpoints x = L - log m, m = Msum..1  (segment between L-log(m+1) and L-log m has terms m' <= m)
    edges = [mp.mpf(0)] + [L - mp.log(m) for m in range(Msum, 0, -1)]
    segs = []
    for a, b in zip(edges[:-1], edges[1:]):
        if b > a:
            # term count on (a,b): m' <= floor(lam^2 e^{-a}) evaluated at midpoint
            mid = (a + b) / 2
            segs.append((a, b, int(mp.floor(lam ** 2 * mp.exp(-mid)))))
    pts = []  # (x, w, kval)
    for a, b, mm in segs:
        for (t, w) in nodes:
            x = (b - a) / 2 * t + (a + b) / 2
            u = mp.exp(x) / lam
            kv = mp.sqrt(u) * sum(eval_H(g, lam, m * u) for m in range(1, mm + 1))
            pts.append((x, (b - a) / 2 * w, kv))
    cs = []
    for n in range(-N, N + 1):
        acc = mp.mpc(0)
        for x, w, kv in pts:
            acc += w * kv * mp.expj(-2 * mp.pi * n * x / L)
        cs.append(acc / mp.sqrt(L))
    # symmetry defect (xi is even <=> c_n = c_{-n}); k(x) vs k(L-x)
    sym = max(abs(cs[i] - cs[-1 - i]) for i in range(len(cs)))
    return cs, dict(info={k: mp.nstr(v, 12) for k, v in info.items()}, segments=len(segs), nodes=len(pts), sym_defect=mp.nstr(sym, 5),
                    max_abs_c=mp.nstr(max(abs(c) for c in cs), 8))


def analyse(lamstr, N, dps_in, skip_trial=False):
    t0 = time.time()
    dps = dps_in
    while True:
        mp.mp.dps = dps
        lam = parse_lam(lamstr)
        T, L, idx, ks = build(lam, N, dps)
        n = len(idx); d = [mp.mpf(j) for j in idx]
        E, Q = mp.eigsy(T)
        order = sorted(range(n), key=lambda i: E[i])
        eps0 = E[order[0]]
        need = int(2 * abs(mp.log10(abs(eps0))) + 30) if eps0 != 0 else dps
        if need <= dps:
            break
        dps = need + 10
    # shift
    for r in range(n):
        T[r, r] -= eps0
    lam_r = [E[order[i]] - eps0 for i in range(n)]
    e_vec = [[Q[row, order[i]] for row in range(n)] for i in range(n)]       # e_vec[r][j]
    xi = e_vec[0][:]; s = sum(xi); xi = [x / s for x in xi]                    # eta^T xi = 1
    i0 = idx.index(0)
    beta = [d[r] * T[r, i0] for r in range(n)]
    eta = [mp.mpf(1)] * n
    dot = lambda a, b: sum(x * y for x, y in zip(a, b))
    nxi = mp.sqrt(dot(xi, xi)); xh = [x / nxi for x in xi]
    Dxi = [d[j] * xi[j] for j in range(n)]
    PDxi = [Dxi[j] - dot(xh, Dxi) * xh[j] for j in range(n)]
    res = dict(lam=lamstr, N=N, dps=dps, eps0=mp.nstr(eps0, 8), L=mp.nstr(L, 10), prime_powers=ks)
    # --- controls
    xihat = lambda z: 2 * mp.sin(z * L / 2) * sum(x / (z - 2 * mp.pi * j / L) for x, j in zip(xi, idx)) / mp.sqrt(L)
    try:
        z1 = mp.findroot(xihat, mp.mpf("14.134725141734693"))
        res["first_zero_error"] = mp.nstr(abs(z1 - mp.mpf("14.134725141734693790457251983562470270784257115699")), 5)
    except Exception as ex:
        res["first_zero_error"] = "fail: %s" % ex
    tp_eta = sum(dot(e_vec[r], eta) ** 2 / lam_r[r] for r in range(1, n))
    res["eta_Tplus_eta"] = mp.nstr(tp_eta, 8)
    Tp_beta = [sum(dot(e_vec[r], beta) * e_vec[r][j] / lam_r[r] for r in range(1, n)) for j in range(n)]
    nrm = mp.sqrt(dot(PDxi, PDxi))
    res["Tplus_beta_plus_PDxi_rel"] = mp.nstr(mp.sqrt(sum((Tp_beta[j] + PDxi[j]) ** 2 for j in range(n))) / nrm, 5)
    res["beta_dot_xi"] = mp.nstr(dot(beta, xi), 5)
    def Wz(e, z):
        return sum(e[j] / (d[j] - z) for j in range(n)) - sum(xi[j] / (d[j] - z) for j in range(n)) * sum(e)
    R = lambda z: [1 / (d[j] - z) for j in range(n)]
    wr = []
    for z in ZS:
        Rz = R(z)
        t1 = dot(beta, [Rz[j] * xi[j] for j in range(n)]); t2 = 1 - dot(eta, [Rz[j] * Dxi[j] for j in range(n)])
        t3 = dot(eta, [Rz[j] * xi[j] for j in range(n)]); t4 = dot(beta, [Rz[j] * Dxi[j] for j in range(n)])
        Dz = t1 * t2 + t3 * t4
        wr.append(dict(z=str(z), W_e1=mp.nstr(abs(Wz(e_vec[1], z)), 6), wronskian_rel=mp.nstr(abs(Dz) / (abs(t1 * t2) + abs(t3 * t4)), 4)))
    res["controls_z"] = wr
    # --- spectral data
    R6 = min(7, n)
    res["low_eigs"] = [mp.nstr(lam_r[r], 8) for r in range(1, R6)]
    res["eig_ratios"] = [mp.nstr(lam_r[r] / lam_r[r + 1], 4) for r in range(1, R6 - 1)]
    MC = {}
    for z in ZS:
        MC[str(z)] = sum(abs(Wz(e_vec[r], z)) ** 2 / lam_r[r] for r in range(1, n))
    res["M_C"] = mp.nstr(max(MC.values()), 8)
    res["W_e_r"] = {str(z): [mp.nstr(abs(Wz(e_vec[r], z)), 6) for r in range(1, R6)] for z in ZS}
    res["t_spectral"] = round(time.time() - t0)
    if skip_trial:
        return res
    # --- trial vector
    cs, tinfo = trial_vector(lam, N, dps)
    res["trial"] = tinfo
    res["trial"]["l2_c"] = mp.nstr(mp.sqrt(sum(abs(c) ** 2 for c in cs)), 12)
    s_k = sum(cs)
    res["eta_dot_k_raw"] = mp.nstr(abs(s_k), 8)
    k = [c / s_k for c in cs]                                                   # eta^T k = 1
    ov = [sum(e_vec[r][j] * k[j] for j in range(n)) for r in range(n)]          # <e_r,k>
    Tk = [sum(T[a, b] * k[b] for b in range(n)) for a in range(n)]
    kTk = sum(mp.conj(k[a]) * Tk[a] for a in range(n)).real
    res["kTk"] = mp.nstr(kTk, 8)
    split = sum(lam_r[r] * abs(ov[r]) ** 2 for r in range(1, R6))
    res["kTk_split_r<=6"] = mp.nstr(split, 8)
    res["kTk_split_remainder"] = mp.nstr(kTk - split, 8)
    res["overlaps_abs"] = [mp.nstr(abs(ov[r]), 8) for r in range(0, R6)]        # r = 0 is the radical
    res["per_mode_pi"] = {str(z): [mp.nstr(abs(ov[r]) ** 2 * abs(Wz(e_vec[r], z)) ** 2, 6) for r in range(1, R6)] for z in ZS}
    # cluster overlaps for near-degenerate neighbours (lam_r/lam_{r+1} > 0.1)
    clusters, cur = [], [1]
    for r in range(1, R6 - 1):
        if lam_r[r] / lam_r[r + 1] > mp.mpf("0.1"):
            cur.append(r + 1)
        else:
            clusters.append(cur); cur = [r + 1]
    clusters.append(cur)
    res["clusters"] = [dict(modes=c, overlap=mp.nstr(mp.sqrt(sum(abs(ov[r]) ** 2 for r in c)), 8)) for c in clusters]
    res["product_kTk_MC"] = mp.nstr(kTk * max(MC.values()), 8)
    pk = [k[j] - (sum(xh[i] * k[i] for i in range(n))) * xh[j] for j in range(n)]
    res["dist_k_from_radical"] = mp.nstr(mp.sqrt(sum(abs(x) ** 2 for x in pk)), 8)
    res["overlap_k_xi_hat"] = mp.nstr(abs(sum(xh[i] * k[i] for i in range(n))), 8)
    nP = mp.sqrt(dot(PDxi, PDxi))
    res["PDxi_dot_k_over_norm"] = mp.nstr(abs(sum(PDxi[j] * k[j] for j in range(n))) / nP, 8)
    res["t_total"] = round(time.time() - t0)
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--lam", required=True); ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--dps", type=int, default=0); ap.add_argument("--out", required=True)
    ap.add_argument("--no-trial", action="store_true")
    a = ap.parse_args()
    dps = a.dps or (60 + 3 * a.N)
    r = analyse(a.lam, a.N, dps, skip_trial=a.no_trial)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(r, open(a.out, "w"), indent=1, default=str)
    print(json.dumps(r, indent=1, default=str)[:3500])
