## This module contains various techniques for special function evaluation, e.g.
## power series, continued fractions, polynomial & Padé approximation, etc.

from fpUt import safeDivisor

template lentz*(F: type,                # FP type
                num, den: untyped,      # numerator, denominator coef series
                f, n, est: untyped,     # var F: running val, it & estim REL.ERR
                tol=1e-7,               # minimum convergence error tolerance
                maxIt=5000,             # maximum number of coefs
                eps=1e-30,              # replacement for 0.0 in x/0.0 contexts
                den0=den(0)) =          # if this is simpler than B(n) (eg. 0.0)
  ## Continued fraction evaluator by modified Lentz method.  This is a template
  ## in terms of templates `num(n)` ,`den(n)` since continued fraction series
  ## have variadix mixed type (int/float) inputs.  Leaving all that type/param
  ## count fixing to a parent proc/func seems an ok way to abstract over that.
  ## NOTE: num(n) and den(n) are indexed starting from 1 as this is the common
  ## convention for series coefficient formulae.
  n = 1              # iteration number
  var A,B,C,D,scl: F # num,den,running C=PN(n)/PN(n-1),D=PD(n)/PD(n-1),updateFac
  C = den0
  if (f = den0; f < F(eps)):            # Handle 1st step specially
    B = F(1) / safeDivisor(den(1), F(eps))
    C = F(0.5) / F(eps)                 # Suppress num(2) for C; could be `inf`
    D = B
    f = num(1) * D
    n = 2
  while n < maxIt:
    A = num(n)
    B = den(n)
    C = B + A / safeDivisor(C, F(eps))
    D = B + A * D
    D = F(1) / safeDivisor(D, F(eps))
    scl = C * D
    f *= scl
    est = abs(scl - F(1))
    if est < F(tol):
      break
    inc n
