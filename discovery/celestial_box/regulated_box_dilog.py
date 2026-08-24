"""Closed-form regulated scalar box from the celestial cut/dispersion representation.

Starting point (Toupin, Loops from Cuts in Celestial Holography):

  C(s',u,m) = (1/8pi) * 4/sqrt(d(d+4c)) * artanh(sqrt(d/(d+4c)))
  d = s' U, c = m (s' + m), U=-u>0, m=mu^2>0

and in the Euclidean region s=-S<0, u=-U<0,

  J(-S,-U) = 8pi int_0^inf C(s',-U,m)/(s'+S) ds'.

The substitution

  r^2 = U s' / ((U+4m)s' + 4m^2)

reduces the dispersion integral exactly to

  J = 8/(S U) int_0^R artanh(r)/(1-kappa^2 r^2) dr,
  R^2 = U/(U+4m),
  kappa^2 = (S(U+4m)-4m^2)/(S U).

Writing y=(1+r)/(1-r), alpha=(1-kappa)/(1+kappa), the primitive is
a difference of ordinary dilogarithms.  This script verifies the closed form
against both the transformed dispersion integral and an independent direct
Feynman-parameter evaluation.
"""

import mpmath as mp

mp.mp.dps = 60


def G(c, d):
    w = d + 4*c
    return 8*mp.atanh(mp.sqrt(d/w))/(mp.sqrt(d)*w**mp.mpf('1.5')) + 2/(c*w)


def J_direct(S, U, m):
    """Independent Euclidean Feynman-parameter box; S=-s>0, U=-u>0, m=mu^2."""
    def inner(x2, x4):
        rho = 1-x2-x4
        if rho <= 0:
            return mp.mpf('0')
        c0 = m*(x2+x4) + U*x2*x4
        d0 = S*rho**2
        return rho*G(c0, d0)
    return mp.quad(
        lambda x2: mp.quad(lambda v: (1-x2)*inner(x2, (1-x2)*v), [0, .5, 1]),
        [0, .01, .5, .99, 1]
    )


def J_reduced(S, U, m):
    R = mp.sqrt(U/(U+4*m))
    k2 = (S*(U+4*m)-4*m*m)/(S*U)
    return 8/(S*U)*mp.quad(lambda r: mp.atanh(r)/(1-k2*r*r), [0, R])


def F(q, y):
    return mp.log(y)*mp.log(1+y/q) + mp.polylog(2, -y/q)


def J_dilog(S, U, m):
    R = mp.sqrt(U/(U+4*m))
    k2 = (S*(U+4*m)-4*m*m)/(S*U)
    k = mp.sqrt(k2)
    if abs(k) < mp.mpf('1e-40'):
        # Degenerate case: denominator in reduced integral is 1.
        return 8/(S*U)*mp.quad(lambda r: mp.atanh(r), [0, R])
    y = (1+R)/(1-R)
    alpha = (1-k)/(1+k)
    bracket = (F(alpha, y)-F(1/alpha, y)
               - mp.polylog(2, -1/alpha) + mp.polylog(2, -alpha))
    return 2*bracket/(k*S*U)


def close(a, b, tol=mp.mpf('1e-45')):
    return abs(a-b) <= tol*max(1, abs(a), abs(b))


def check_point(S, U, m):
    jd = J_direct(S, U, m)
    jr = J_reduced(S, U, m)
    jc = J_dilog(S, U, m)
    assert close(jd, jr, mp.mpf('1e-40'))
    assert close(jd, jc, mp.mpf('1e-40'))
    print(f'S={S}, U={U}, mu^2={m}')
    print('  direct =', mp.nstr(jd, 50))
    print('  reduced=', mp.nstr(jr, 50))
    print('  dilog  =', mp.nstr(jc, 50))


if __name__ == '__main__':
    # Paper anchor plus two generic Euclidean points.
    check_point(mp.mpf(3), mp.mpf(2), mp.mpf('.5'))
    check_point(mp.mpf('4.7'), mp.mpf('1.6'), mp.mpf('.3'))
    check_point(mp.mpf('2.2'), mp.mpf('3.1'), mp.mpf('.7'))
    print('PASS: celestial cut dispersion -> reduced artanh integral -> closed dilogarithm box')
