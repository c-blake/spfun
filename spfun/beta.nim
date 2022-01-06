## Provide complete & normalized/regularized incomplete beta functions.
## The latter is done via evaluation of the usual continued frac series:
##          x^a (1-x)^b   1 c1 c2  WHERE numerator coefs c_i are:
##  Ix(a,b)=----------- [-- -- --..  c_{2n+1} = -x(a+n)(a+b+n)/((a+2n)(a+2n+1))
##            a B(a,b)   1+ 1+ 1+    c_{2n}   =  x   n(b-n)   /((a+2n)(a+2n-1))

import math, eval

var doNotUse: float64 # Only for default value; Callers must provide `est`.

func lnBeta*[F](a, b: F): F {.inline.} = lgamma(a) + lgamma(b) - lgamma(a + b)
  ## Natural log of the complete beta function; often written `ln B(a,b)`.

func betaCF[F](x, a, b, err: F; est: var float64=doNotUse): F {.inline.} =
  template den(n): untyped = F(1)
  template num(n): untyped =
    let m = F((n - 1) div 2)
    if n == 1: F(1)
    elif (n and 1) != 0: x * m * (b - m) / ((a + F(2)*m) * (a + F(2)*m - F(1)))
    else   : -(x * (a + m) * (a + b + m) / ((a + F(2)*m) * (a + F(2)*m + F(1))))
  var val: F
  var it: int
  lentz(F, num, den, val, it, est, err, den0=F(0))
# if err > 0.005: ({.cast(noSideEffect).}: echo "its: ", it) # track cost
  if est < err: val else: F(0)          # NAN on convergence failure?

func betaI*[F](x, a, b: F; err: F=F(1e-6), est: var float64=doNotUse): F =
  ## Regularized incomplete beta function.  `B(x,a,b)/B(a,b)` where `B(x,a,b) =
  ## integral(0, x, dt*t^(a-1)*(1-t)^(b-1)`, sometimes written `B_x(a,b)`.
  if x == F(0) or x == F(1):
    result = x
    est = 0f64
  elif x < F(0) or x > F(1):
    result = NaN
    est = 0f64
  else:
    let xC = F(1) - x
    let B = exp(a * ln(x) + b * ln(xC) - lnBeta(a, b))
    var est: float64
    result =
      if x*(a + b + F(2)) < a + F(1):        B * betaCF(x , a, b, err, est) / a
      else                          : F(1) - B * betaCF(xC, b, a, err, est) / b
    est *= result                       # relative -> absolute error estimate

when isMainModule:
  import fpUt
  assert almostEqual(lnBeta(2.0, 3.0), -2.484906649788, 1e-7)
  assert almostEqual(lnBeta(2.0f32, 3.0f32), -2.484906649788f32, 1e-7)
  assert almostEqual(betaI(0.5, 2.0, 3.0), 0.6875, 1e-7)
  assert almostEqual(betaI(0.5f32, 2.0f32, 3.0f32), 0.6875f32, 1e-7)
