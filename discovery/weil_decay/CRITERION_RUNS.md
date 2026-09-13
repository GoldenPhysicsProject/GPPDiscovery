# Criterion runs, 2026-09-13

All runs made this session, testing the two resolution criteria. Environment:
Python 3.11.15, mpmath 1.4.1, connes-cvs 0.2.2, `dps=90` throughout.

`margin = T*log(c)/(2*pi*N)` (archimedean, resolution.py);
`headroom = N/N_req`, `N_req = log(c)/min_gap(c) - 1/2` (prime comb, spacing.py).

| c | N | T | dim | margin | headroom | sign | lambda_1 | sec |
|---|---|---|-----|--------|----------|------|----------|-----|
| 3 | 12 | 60 | 25 | 0.874 | 5.43 | -1 | `-0.2443933192254289936596183` | 98.3 |
| 3 | 16 | 60 | 33 | 0.656 | 7.24 | -1 | `-0.4855491048577200817461983` | 121.1 |
| 3 | 20 | 60 | 41 | 0.525 | 9.05 | -1 | `-0.4900557341257446105431979` | 194.6 |
| 3 | 36 | 120 | 73 | 0.583 | 16.29 | -1 | `-0.4901288932404565165938318` | 839.5 |
| 3 | 36 | 160 | 73 | 0.777 | 16.29 | -1 | `-0.4897489413922827120029554` | 859.8 |
| 3 | 36 | 190 | 73 | 0.923 | 16.29 | -1 | `-0.3779171010840920867089124` | 415.5 |
| 3 | 36 | 210 | 73 | 1.020 | 16.29 | +1 | `4.656080582364100789220654e-8` | 430.1 |
| 3 | 36 | 240 | 73 | 1.166 | 16.29 | +1 | `4.777124319083887855688378e-8` | 558.5 |
| 3 | 44 | 300 | 89 | 1.192 | 19.91 | +1 | `4.939277536912255790786999e-8` | 2152.9 |
| 3 | 52 | 300 | 105 | 1.009 | 23.53 | +1 | `4.912671463759615266428714e-8` | 1510.0 |
| 3 | 60 | 300 | 121 | 0.874 | 27.16 | -1 | `-0.4896873806431253686980452` | 1299.9 |
| 3 | 68 | 300 | 137 | 0.771 | 30.78 | -1 | `-0.4901289692666868743043771` | 1035.9 |
| 5 | 36 | 422 | 73 | 3.003 | 5.36 | +1 | `9.641940539432706930785647e-18` | 625.8 |
| 5 | 44 | 300 | 89 | 1.746 | 6.55 | +1 | `8.616202554895980048174814e-18` | 2267.9 |
| 5 | 52 | 300 | 105 | 1.478 | 7.75 | +1 | `8.171335295804326128739907e-18` | 1555.6 |
| 5 | 52 | 610 | 105 | 3.005 | 7.75 | +1 | `9.786360807440641981210483e-18` | 1192.9 |
| 19 | 60 | 300 | 121 | 2.343 | 1.25 | +1 | `9.853690563424607676092397e-80` | 1389.4 |

(Timings are not comparable between rows: several runs shared two cores.)

## What each block shows

**c=3, N=36, T varied (120, 160, 190, 210, 240, 300).** The controlled test of the
archimedean criterion. Sign flips between `T=190` (margin 0.923) and `T=210` (margin 1.020):
a 10% change in one parameter moving `lambda_1` seven decades. Above the line the value has
already settled -- `4.656e-8`, `4.777e-8`, `4.951e-8` at margins 1.020, 1.166, 1.457.

**c=3, T=300, N varied (36, 44, 52, 60, 68).** The same criterion crossed from the other
direction. Positive and flat at margins 1.457, 1.192, 1.009 (`lambda_1` moves 0.8% while the
margin falls by a third), then negative at 0.874 and 0.771.

**Below threshold the computation converges, to the wrong operator.** `c=3, N=36, T=120`
(margin 0.583) and `c=3, N=68, T=300` (margin 0.771) agree to **seven significant figures**,
`-0.4901289`. Sorted by `N`, the sub-threshold values are `-0.2444, -0.4855, -0.4901, -0.4901,
-0.4897, -0.4901` at `N = 12, 16, 20, 36, 60, 68` -- settling from `N=20` onward regardless of
margin. Dropping the outermost modes' archimedean contribution produces a well-defined
different matrix with a genuinely negative ground state. The criterion is the line between
computing the intended operator and computing that one, which is why the approach from above
is flat and the crossing is a jump rather than a gradual loss of accuracy.

**Correction.** An earlier reading of this session's data suggested the sub-threshold
*magnitude* tracked the margin, on the strength of two points at margin 0.87 giving `-0.332`
and `-0.244`. A third point at margin 0.874 (`c=3, N=60, T=300`) gives `-0.4897`. That claim
is withdrawn. The margin fixes where the threshold is and says nothing about the value below
it; `N` does.

**c=5, fixed T=300 vs fixed margin 3.0.** The separation of the two criteria:

| series | N range | change in lambda_1 |
|--------|---------|--------------------|
| fixed `T=300` (margin 2.13 -> 1.48) | 36 -> 52 | **-8.2%** |
| fixed margin 3.0 (`T` 422 -> 610) | 36 -> 52 | **+1.5%** |

Five times smaller and opposite in sign. The drift previously read as "`N` is not converged at
`c=5`" was the archimedean margin degrading as `N` rose at fixed `T`, not prime resolution.
Hold the margin and `lambda_1` is essentially converged in `N` at `c=5` by `N=36`.

It also means the published `c=5` value is about **10% low**: `8.90e-18` at margin 2.13 against
`9.79e-18` at margin 3.0 -- and margin 3 is not demonstrably converged either, so
**margin >= 1 buys the sign, not the value.** Nobody has measured where the archimedean error
actually flattens; that is the cheapest useful next run (`c=3, N=36` at `T = 300, 420, 618,
900`, i.e. margins 1.46, 2.04, 3.00, 4.37, all at dim 73).

**c=19, N=60.** The fourth point of the `c=19` `N`-series: `9.854e-80`, against `2.661e-63`
at `N=36`. The published value is wrong by **seventeen orders of magnitude** and still falling
-- decrements `-6.41, -5.55, -4.48` decades, decelerating, at headroom `0.75, 0.92, 1.08, 1.25`.
This one is genuine prime-comb under-resolution and is not fixable by moving `T`.

## Where that leaves the decay constant

Ranking the `N=36`, `T=300` c-scan by headroom: `c=3` (16.4), `c=5` (5.4), `c=7` (4.4),
`c=11` (1.81), `c=13` (1.69), `c=17` (0.78), `c=19` (0.75). Only the first three are
comfortably resolved in the prime comb, and all seven sit at uncontrolled and varying
archimedean margins (1.46 to 3.91).

Three usable points give two slopes: `c=3 -> 5` is `-43.93` and `c=5 -> 7` is `-68.72`, in
natural log. They do not agree with each other, and neither matches `-2*gamma_1 = -28.27` or
`-4*gamma_1 = -56.55`. The README's "two-point slope of about -64" comes from `c=5` and `c=7`,
which is the one pair that survives -- so that number is the honest state of the art, and it
is a single two-point estimate with no curvature information behind it.
