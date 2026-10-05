"""Aggregate trial-vector overlap runs into TRIAL_VECTOR_OVERLAP_RESULTS.md.
Evidence only; nothing here is a proof. Reads results/tv_*.json."""
import glob, json, math, os, sys
here = os.path.dirname(os.path.abspath(__file__))
rows = []
for f in sorted(glob.glob(os.path.join(here, "results", "tv_*.json"))):
    d = json.load(open(f))
    rows.append(d)
def fl(x):
    try: return float(x)
    except Exception: return float("nan")
out = ["# CCM trial-vector overlap (Poisson-prolate k_lambda)", "",
       "Evidence only. Each row: one (lambda, N). eta-normalised (eta^T k = 1).", "",
       "| lambda | N | dps | eps0 | first-zero err | kTk | M_C | kTk*M_C | eta.k raw | max |<e_r,k>| (r>=1) | even-mode overlaps (>1e-12) |",
       "|---|---|---|---|---|---|---|---|---|---|---|"]
for d in sorted(rows, key=lambda d: (fl(eval(d["lam"].replace("sqrt","math.sqrt"))) if True else 0, d["N"])):
    ov = [fl(x) for x in d.get("overlaps_abs", [])]
    ev = [f"{o:.2e}" for o in ov[1:] if o > 1e-12]
    mx = max(ov[1:]) if len(ov) > 1 else float("nan")
    out.append("| %s | %d | %s | %s | %s | %s | %s | %s | %s | %.3e | %s |" % (
        d["lam"], d["N"], d["dps"], d["eps0"], d.get("first_zero_error"), d.get("kTk"), d.get("M_C"),
        d.get("product_kTk_MC"), d.get("eta_dot_k_raw"), mx, ", ".join(ev)))
out += ["", "## Low eigenvalues and clusters", ""]
for d in rows:
    out.append("- lam=%s N=%d: low_eigs=%s; clusters=%s" % (d["lam"], d["N"], d.get("low_eigs"),
               [(c["modes"], c["overlap"]) for c in d.get("clusters", [])]))
# convergence in N per lambda: first even overlap and kTk
out += ["", "## Convergence in N (first even overlap, kTk)", ""]
bylam = {}
for d in rows: bylam.setdefault(d["lam"], []).append(d)
for lam, ds in bylam.items():
    ds.sort(key=lambda d: d["N"])
    seq = []
    for d in ds:
        ov = [fl(x) for x in d.get("overlaps_abs", [])][1:]
        ev = [o for o in ov if o > 1e-12]
        seq.append((d["N"], ev[0] if ev else float("nan"), fl(d.get("kTk"))))
    out.append("- lam=%s: %s" % (lam, "; ".join("N=%d: ov1=%.3e kTk=%.4f" % s for s in seq)))
out += ["", "Verdict rule: GO if kTk*M_C stays bounded while the first even overlap decays no faster than the eigenvalue",
        "(alpha = slope of log overlap vs log eigenvalue <= 1/2); NO-GO if kTk*M_C grows with N; else INCONCLUSIVE.",
        "A verdict is only meaningful once the N-sequence has converged (watch the convergence block).", ""]
open(os.path.join(here, "TRIAL_VECTOR_OVERLAP_RESULTS.md"), "w").write("\n".join(out))
print("\n".join(out))
