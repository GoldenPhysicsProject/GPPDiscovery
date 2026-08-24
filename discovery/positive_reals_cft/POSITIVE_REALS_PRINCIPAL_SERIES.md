# Positive reals as a one-dimensional principal-series system

This note isolates an exact representation-theoretic structure behind the working analogy that `(R_+, x)` behaves like a one-dimensional conformal/dilation system.

## 1. Unitary dilation representation on L^2(R_+, dx)

For `a>0`, define

`(U_a f)(x) = a^(1/2) f(a x)`.

Then `U_a` is unitary on `L^2(R_+,dx)`:

`int_0^inf |U_a f(x)|^2 dx = int_0^inf a |f(ax)|^2 dx = int_0^inf |f(y)|^2 dy`.

For the generalized Mellin mode

`f_s(x)=x^(-s)`, `s in C`,

one has

`U_a f_s = a^(1/2-s) f_s`.

Hence, for any fixed nontrivial `a != 1`,

`|a^(1/2-s)| = a^(1/2-Re(s))`,

so

`|a^(1/2-s)|=1  <=>  Re(s)=1/2`.

Thus the Riemann critical line is exactly the unitary spectral line of the ordinary Lebesgue-half-density dilation representation.

## 2. Unitary inversion and the functional-equation reflection

Define

`(J f)(x)=x^(-1) f(1/x)`.

Then `J` is unitary on `L^2(R_+,dx)` and involutive:

`J^2=1`.

Indeed, with `y=1/x`, `dx=dy/y^2`, the factor `x^-2=y^2` cancels the inversion Jacobian.

On Mellin modes,

`J f_s(x)=x^-1 (x^-1)^(-s)=x^(s-1)=f_(1-s)(x)`.

Therefore inversion acts spectrally as

`s -> 1-s`.

On `Re(s)=1/2`, this becomes

`1-s = conjugate(s)`.

So the familiar Riemann functional-equation reflection is exactly the spectral action of the unitary geometric inversion of the positive real line in this half-density representation.

## 3. Relation to multiplicative Haar measure

On `L^2(R_+, dx/x)`, the bare multiplicative characters `x^(it)` are unitary and inversion acts by `t -> -t`; the natural spectral axis is centered at `Re(s)=0`.

Passing from Haar functions to Lebesgue half-densities multiplies by `x^-1/2`, shifting the spectral coordinate by `1/2`. The same unitary axis is therefore represented as

`Re(s)=1/2`.

This explains the critical-line half shift as a measure/half-density effect rather than an arbitrary numerical choice.

## 4. Celestial dictionary

Under the project dictionary `Delta=2s`,

`Re(s)=1/2 <=> Re(Delta)=1`,

and

`s -> 1-s <=> Delta -> 2-Delta`.

Thus the positive-real inversion representation and the celestial scalar shadow involution are the same affine reflection after `Delta=2s`.

This is an exact mathematical dictionary. It does not by itself constrain the zeros of zeta; additional global arithmetic input is still required.

## Lean targets

1. Prove the dilation eigenvalue modulus criterion algebraically for `a>0`, `a!=1`.
2. Define the inversion half-density map and prove `J^2=id` pointwise.
3. Formalize the L^2 norm preservation under inversion once the required change-of-variables API is identified.
4. Bridge `Delta=2s` to the existing celestial shadow/reflection modules.
