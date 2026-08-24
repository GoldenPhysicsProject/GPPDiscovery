# Analytic mass-regulator limit of the celestial-dispersion scalar box

Starting from `regulated_box_dilog.py`, write

\[
J(-S,-U;m)=\frac{2}{\kappa S U}B,\qquad m=\mu^2>0,
\]

with

\[
R=\sqrt{\frac{U}{U+4m}},\quad
\kappa^2=\frac{S(U+4m)-4m^2}{SU},\quad
y=\frac{1+R}{1-R},\quad
\alpha=\frac{1-\kappa}{1+\kappa}.
\]

For the Euclidean small-`m` regime, `κ>1`, so put

\[
b=-\alpha=\frac{\kappa-1}{\kappa+1}>0,\qquad p=by.
\]

The original principal-branch dilogarithm bracket is

\[
B=\log y\,[\log(1-y/b)-\log(1-by)]
+\operatorname{Li}_2(y/b)-\operatorname{Li}_2(by)
-\operatorname{Li}_2(1/b)+\operatorname{Li}_2(b).
\]

Here `y/b>1`, `1/b>1`, and `0<p=by<1`. Apply the principal-branch inversion identity

\[
\operatorname{Li}_2(x)+\operatorname{Li}_2(x^{-1})
=\frac{\pi^2}{3}-\frac12\log^2x-i\pi\log x,
\qquad x>1,
\]

and

\[
\log(1-x)=\log(x-1)+i\pi,\qquad x>1.
\]

The imaginary terms cancel exactly. After elementary logarithmic simplification one gets the real identity

\[
\boxed{
B=(\log p-\log b)\log\frac{p-b^2}{1-p}
-\operatorname{Li}_2\!\left(\frac{b^2}{p}\right)
+2\operatorname{Li}_2(b)-\operatorname{Li}_2(p)
-\frac12\log^2p+\frac12\log^2b .
}
\]

This form makes the regulator limit transparent. Direct expansion of the algebraic variables gives

\[
R=1-\frac{2m}{U}+O(m^2),
\]

\[
\kappa=1+\frac{2m}{U}+O(m^2),
\]

\[
b=\frac{m}{U}+O(m^2),
\]

and, importantly,

\[
p=by=1-\frac{m}{S}+O(m^2).
\]

Therefore, with `L=log m`,

\[
\log b=L-\log U+o(1),\qquad
\log p=o(1),
\]

\[
\log\frac{p-b^2}{1-p}=-L+\log S+o(1),
\]

while

\[
\operatorname{Li}_2(b^2/p),\operatorname{Li}_2(b)\to0,
\qquad
\operatorname{Li}_2(p)\to\frac{\pi^2}{6}.
\]

Hence

\[
B=\frac32L^2-(\log S+2\log U)L
+\log S\log U+\frac12\log^2U-\frac{\pi^2}{6}+o(1).
\]

Since `κ→1`, multiplying by `2/(κSU)` gives

\[
\boxed{
J(-S,-U;m)=\frac1{SU}\left[
3\log^2m-(2\log S+4\log U)\log m
+\log^2U+2\log S\log U-\frac{\pi^2}{3}
\right]+o(1).
}
\]

This is an analytic consequence of the celestial cut + dispersion + dilogarithm closure. It is not obtained by fitting. The earlier numerical coefficient extraction merely anticipated the exact expression.

## Scope

This is the mass-regulator asymptotic of the specific internally regulated scalar box used in `Loops_from_Cuts_in_Celestial_Holography.tex`. It should not be silently identified with a dimensional-regularization finite part; matching schemes requires an explicit regulator dictionary.
