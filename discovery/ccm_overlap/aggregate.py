"""Aggregate trial-vector overlap runs into TRIAL_VECTOR_OVERLAP_RESULTS.md.
Evidence only; nothing here is a proof. Reads results/tv_*.json."""
import glob, json, math, os
here = os.path.dirname(os.path.abspath(__file__))
rows = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(here, "results", "tv_*.json")))]

def fl(x):
    try: return float(x)
    except Exception: return float("nan")

def lamval(s): return eval(s.replace("sqrt", "math.sqrt"))

out = ["# CCM trial-vector overlap (Poisson-prolate k_lambda)", "",
       "Evidence only. Each row is one (lambda, N). eta-normalised (eta^T k = 1).", "",
       "| lambda | N | dps | eps0 | first-zero err | kTk | M_C | kTk*M_C | eta.k raw | max overlap r>=1 | even-mode overlaps (>1e-12) |",
       "|---|---|---|---|---|---|---|---|---|---|---|"]
for d in sorted(rows, key=lambda d: (lamval(d["lam"]), d["N"])):
    ov = [fl(x) for x in d.get("overlaps_abs", [])]
    ev = [f"{o:.2e}" for o in ov[1:] if o > 1e-12]
    mx = max(ov[1:]) if len(ov) > 1 else float("nan")
    out.append("| %s | %d | %s | %s | %s | %s | %s | %s | %s | %.3e | %s |" % (
        d["lam"], d["N"], d["dps"], d["eps0"], d.get("first_zero_error"), d.get("kTk"),
        d.get("M_C"), d.get("product_kTk_MC"), d.get("eta_dot_k_raw"), mx, ", ".join(ev)))

out += ["", "## Low eigenvalues and clusters", ""]
for d in rows:
    out.append("- lam=%s N=%d: low_eigs=%s; clusters=%s" % (
        d["lam"], d["N"], d.get("low_eigs"), [(c["modes"], c["overlap"]) for c in d.get("clusters", [])]))

out += ["", "## Convergence in N (first even overlap, kTk)", ""]
bylam = {}
for d in rows: bylam.setdefault(d["lam"], []).append(d)
for lam, ds in bylam.items():
    ds.sort(key=lambda d: d["N"]); seq = []
    for d in ds:
        ev = [o for o in (fl(x) for x in d.get("overlaps_abs", [])[1:]) if o > 1e-12]
        seq.append((d["N"], ev[0] if ev else float("nan"), fl(d.get("kTk"))))
    out.append("- lam=%s: %s" % (lam, "; ".join("N=%d: ov1=%.3e kTk=%.4f" % s for s in seq)))

out += ["", "## Scale-free reading (eta-normalisation is ill-posed: eta.k_raw -> 0 as N grows)", "",
        "kTk_sf = kTk*(eta.k_raw)^2/max|c|^2 ; cos_r <= |ov_r|*|eta.k_raw|/max|c| (upper bound on the cosine with mode r).", "",
        "| lambda | N | kTk_sf | M_C | kTk_sf*M_C | cos(e2) bound | cos(e4) bound | lambda_even1 |", "|---|---|---|---|---|---|---|---|"]
for d in sorted(rows, key=lambda d: (lamval(d["lam"]), d["N"])):
    s_ = fl(d.get("eta_dot_k_raw")); mx_ = fl(d["trial"]["max_abs_c"])
    ov = [fl(x) for x in d.get("overlaps_abs", [])]; le = [fl(x) for x in d.get("low_eigs", [])]
    evi = [i for i in range(1, len(ov)) if ov[i] > 1e-12]
    ksf = fl(d.get("kTk")) * s_ * s_ / mx_ ** 2
    cs_ = [ov[i] * s_ / mx_ for i in evi]
    out.append("| %s | %d | %.3e | %.3e | %.3e | %s | %s | %s |" % (d["lam"], d["N"], ksf, fl(d.get("M_C")), ksf * fl(d.get("M_C")),
               "%.2e" % cs_[0] if cs_ else "-", "%.2e" % cs_[1] if len(cs_) > 1 else "-", "%.1e" % le[evi[0] - 1] if evi else "-"))
out += ["", "Verdict rule: GO if kTk*M_C stays bounded while the first even overlap decays no faster",
        "than the eigenvalue (alpha <= 1/2); NO-GO if kTk*M_C grows with N; else INCONCLUSIVE.",
        "A verdict only counts once the N-sequence has converged.", ""]
open(os.path.join(here, "TRIAL_VECTOR_OVERLAP_RESULTS.md"), "w").write("\n".join(out))
print("\n".join(out))
