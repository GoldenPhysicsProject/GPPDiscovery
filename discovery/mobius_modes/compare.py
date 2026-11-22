"""Task 3(ii): kappa (zeros) next to the Q(L)/L slope (sieve).  Usage: python compare.py direct_summary.json kappa_summary.json OUT.md"""
import sys, json
d = json.load(open(sys.argv[1])); k = json.load(open(sys.argv[2]))
loc = d['local_slopes_decade']
with open(sys.argv[3], 'w') as f:
    f.write("# kappa (zero side) against the Q(L)/L plateau (Moebius side) (NUMERICAL)\n\n")
    f.write("Two independent computations of the same constant, the mean square of B(u) per unit u: the sieve gives the slope of "
            "Q(L) = integral_0^L B^2; the zeros give kappa = 2 sum_(gamma>0) 1/(gamma^2 |zeta'(rho)|^2). Agreement is evidence for "
            "the almost-periodic (explicit-formula) picture and says nothing about proving RH.\n\n")
    f.write("| quantity | value | note |\n|---|---|---|\n")
    f.write(f"| kappa (zeros, {k['zeros']:,} of them, gamma <= {k['gamma_max']:.0f}) | **{k['kappa']:.6f}** +- {k['kappa_err']:.1e} | partial {k['kappa_partial']:.6f} + tail {k['tail_estimate']:.2e} |\n")
    f.write(f"| Q(L) least-squares slope (sieve, N from 1e4 to {d['N_max']:,}) | **{d['ls_slope_from_1e4']:.6f}** | residual rms {d['ls_resid_rms']:.4f} (almost-periodic wobble of Q around the line) |\n")
    f.write(f"| Q/L at N = {d['N_max']:,} | {d['Q_over_L_at_Nmax']:.6f} | approaches the slope like intercept/L, intercept {d['ls_intercept']:.4f} |\n")
    f.write(f"| Q/L^2 at N = {d['N_max']:,} | {d['Q_over_L2_at_Nmax']:.6f} | goes to 0: linear, not quadratic, growth |\n")
    diff = d['ls_slope_from_1e4'] - k['kappa']
    f.write(f"\nDifference slope - kappa = {diff:+.6f} ({diff / k['kappa'] * 100:+.2f}% of kappa).\n\n## Local slopes of Q(L) over one decade of N (sieve)\n\n| N from | N to | slope |\n|---|---|---|\n")
    for a, b, s in loc: f.write(f"| {a:,} | {b:,} | {s:.5f} |\n")
print(f"kappa={k['kappa']:.6f}  slope={d['ls_slope_from_1e4']:.6f}")
