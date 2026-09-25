"""Claude (chat seat) - shape of the fixed-lambda compensated-Weyl remainder R_N(z) = C_N(z) - G(z).
Signed complex R on z in [0.05,12] (real) and z = x + i, for several N.
(1) rank-one test: R_N(z) / R_Nref(z) constant in z?
(2) least-squares fit of R_Nref against candidate channels.
Finding under test (sandbox, lambda 3 N 18): R ~ kappa_N * z * sin^2(pi z), rel residual 0.76%."""
import mpmath as mp, json, argparse, time
from rinf import setup

def pieces(lamstr, N, dps):
    L, bW, T, idx = setup(lamstr, N, dps); n = len(idx); pi = mp.pi
    E, Q = mp.eigsy(T); k = min(range(n), key=lambda i: E[i]); eps = E[k]
    xi = [Q[r, k] for r in range(n)]; s = sum(xi); xi = [x/s for x in xi]
    beta = [idx[a]*T[a, idx.index(0)] for a in range(n)]
    A  = lambda z: sum(xi[j]/(idx[j]-z) for j in range(n))
    Ap = lambda z: sum(xi[j]/(idx[j]-z)**2 for j in range(n))
    B  = lambda z: sum(beta[j]*xi[j]/(idx[j]-z) for j in range(n))
    h  = lambda z: mp.sin(2*pi*z)/(2*pi)
    s2 = lambda z: mp.sin(pi*z)**2
    def row(z):
        a, ap, b = A(z), Ap(z), B(z)
        G = bW(z) - (eps+L)*h(z); C = b/a + (L/pi**2)*s2(z)*ap/a
        return {"z": z, "R": C-G, "G": G, "cand": {
            "h": h(z), "sin2": s2(z), "cos2pz": mp.cos(2*pi*z), "one": mp.mpf(1),
            "sin2_Ap_over_A": s2(z)*ap/a, "B_over_A": b/a, "inv_A": 1/a, "sin2_over_A": s2(z)/a,
            "h_over_A": h(z)/a, "z_sin2": z*s2(z), "z2_sin2": z*z*s2(z)}}
    return eps, L, row

def zgrid():
    zs = [mp.mpf(k)/20 for k in range(1, 241) if k % 10 != 0]
    zs += [mp.mpc(mp.mpf(k)/4, 1) for k in range(0, 49)]
    return zs

def lstsq(rows, keys):
    M = mp.matrix(len(rows), len(keys)); y = mp.matrix(len(rows), 1)
    for i, r in enumerate(rows):
        w = 1/abs(r["G"])
        for j, k in enumerate(keys): M[i, j] = r["cand"][k]*w
        y[i] = r["R"]*w
    MH = M.H; c = mp.lu_solve(MH*M, MH*y); res = M*c - y
    rel = mp.sqrt(sum(abs(x)**2 for x in res))/mp.sqrt(sum(abs(x)**2 for x in y))
    return [c[j] for j in range(len(keys))], rel

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--lam", required=True)
    ap.add_argument("--Ns", default="24,36,48"); ap.add_argument("--out", required=True)
    a = ap.parse_args(); Ns = [int(x) for x in a.Ns.split(",")]; Nref = Ns[-1]; t0 = time.time()
    zs = zgrid(); data = {}; kappa = []
    for N in Ns:
        eps, L, row = pieces(a.lam, N, 60+2*N)
        data[N] = [row(z) for z in zs]
        c1, r1 = lstsq(data[N], ["z_sin2"]); c3, r3 = lstsq(data[N], ["z_sin2", "sin2", "z2_sin2"])
        kappa.append({"N": N, "eps": mp.nstr(eps, 8), "kappa": mp.nstr(c1[0], 10), "rel_resid_1": mp.nstr(r1, 4),
                      "coef_3": [mp.nstr(x, 8) for x in c3], "rel_resid_3": mp.nstr(r3, 4)})
        print(kappa[-1], "t", round(time.time()-t0), flush=True)
    ref = data[Nref]
    rank1 = {}
    for N in Ns[:-1]:
        q = [data[N][i]["R"]/ref[i]["R"] for i in range(len(zs)) if abs(ref[i]["R"]) > 1e-30]
        mean = sum(q)/len(q); spread = max(abs(x-mean) for x in q)/abs(mean)
        rank1[str(N)] = {"mean_ratio": mp.nstr(mean, 8), "rel_spread": mp.nstr(spread, 4)}
    sets = [["h"], ["sin2"], ["cos2pz"], ["one"], ["sin2_Ap_over_A"], ["B_over_A"], ["inv_A"], ["sin2_over_A"], ["h_over_A"],
            ["h", "sin2"], ["sin2_Ap_over_A", "B_over_A"], ["inv_A", "sin2_over_A"], ["sin2_over_A", "h_over_A"],
            ["h", "sin2", "cos2pz", "one"], ["inv_A", "sin2_over_A", "h_over_A"],
            ["z_sin2"], ["z_sin2", "sin2"], ["z_sin2", "sin2", "z2_sin2"],
            ["h", "sin2", "cos2pz", "one", "inv_A", "sin2_over_A", "h_over_A", "sin2_Ap_over_A", "B_over_A", "z_sin2"]]
    fits = []
    for S in sets:
        c, rel = lstsq(ref, S)
        fits.append({"set": S, "coef": [mp.nstr(x, 10) for x in c], "rel_residual": mp.nstr(rel, 5)})
        print(fits[-1], flush=True)
    samples = [{"z": mp.nstr(r["z"], 6), "R": mp.nstr(r["R"], 10), "G": mp.nstr(r["G"], 10)} for r in ref]
    json.dump({"lam": a.lam, "Ns": Ns, "kappa_track": kappa, "rank1_vs_Nref": rank1, "fits": fits, "samples_Nref": samples,
               "sec": round(time.time()-t0)}, open(a.out, "w"), indent=1)
    print(json.dumps(rank1))
