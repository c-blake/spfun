## This module contains various techniques for special function evaluation, e.g.
## power series, continued fractions, polynomial & Padé approximation, etc.

from fpUt import safeDivisor

template lentz*(F: type,                # FP type
                num, den: untyped,      # numerator, denominator coef series
                f, n, est: untyped,     # var F: running val, it & estim REL.ERR
                tol=1e-7,               # minimum convergence error tolerance
                minIt=2,                # minimum number of coefs
                maxIt=5000,             # maximum number of coefs
                eps=1e-30,              # replacement for 0.0 in x/0.0 contexts
                den0=den(0)) =          # if this is simpler than B(n) (eg. 0.0)
  ## Continued fraction evaluator by modified Lentz method.  I.e., evaluate
  ##               num1   num2
  ##   f = den0 + ------ ------ .. = den0 + num1/(den1 + num2/(den2 + ..))
  ##              den1 + den2 +
  ## This is a template in terms of templates `num(n)` ,`den(n)` since continued
  ## fraction series have variadic mixed type (int/float) inputs.  Leaving all
  ## that type/param count fixing to a parent proc/func seems an ok way to
  ## abstract over that.  NOTE: num(0) is unused.
  var A,B,C,D,scl: F # num,den,running C=PN(n)/PN(n-1),D=PD(n)/PD(n-1),updateFac
  f = den0
  if f.abs < F(eps): f = eps # Could copySign; Shouldn't matter since it's eps
  C = f
  D = F(0)
  n = 1              # iteration number
  while n < maxIt:
    A = num(n)
    B = den(n)
    C = B + A / safeDivisor(C, F(eps))
    D = B + A * D
    D = F(1) / safeDivisor(D, F(eps))
    scl = C * D
    f *= scl
    est = abs(scl - F(1))
    if est < F(tol) and n >= minIt:
      break
    inc n
