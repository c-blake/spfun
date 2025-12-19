## Normalized/regularized incomplete gamma function.
import math, eval

var doNotUse: float64 # Only for default value; Callers must provide `est`.

func lnGamma*[F](x: F; err: F=F(1e-6), est: var float64=doNotUse): F {.inline.}=
  ## Natural log of the complete gamma function; Alias for `lgamma`.
  lgamma(x) # NOTE: Could add some fast, less accurate methods here
  #XXX Should put something for `estp`.

func series[F](a, x, err: F; est: var float64=doNotUse): F {.inline.} =
  var a = a                             # Editable param
  var t = F(1) / a                      # Initial running series term
  template init: untyped = t            # g(a, x)*e^x*x^-a =
  template next(n): untyped =           #     sum(n=0.., G(a)*x^n/G(a + n + 1))
    a += F(1); t *= x / a; t            #.. with G(z + 1) = z*G(z)
  var val: F; var it: int               #.. and G(a,x) = g(a,x)/G(a).
  powSeries(F, init, next, val, it, est, err)
  if est < err*val.abs: val else: F(0)  # NaN on convergence failure?

func conFrac[F](a, x, err: F; est: var float64=doNotUse): F {.inline.} =
  template num(n): untyped =
    if n == 1: F(1) else: F(n-1) * (a - F(n-1)) #   1     1*(a-1) 2*(a-2)
  template den(n): untyped = x - a + F(2*n - 1) # ------- ------- ------- ..
  var val: F; var it: int                       # x-a+1 + x-a+3 + x-a+5 +
  lentz0 F, num, den, F(0), val, it, est, err
  if est < err: val else: F(0)          # NaN on convergence failure?

proc gammaI*[F](a, x: F; err: F=F(1e-6), est: var float64=doNotUse,
                norm: ptr F=nil): F =
  ## Normalized incomplete gamma function: `integ(0,x, dt*exp(-t)*t^(a-1)/G(a)`.
  ## Sometimes this is called `P(a,x)`.  If given, `norm[] = lnGamma(a)`.
  let lnG = lnGamma(a, err)             #XXX assume calculated exactly
  let fac = exp(-x + a*ln(x) - lnG)
  if x < F(0) or a < F(0): result = NaN # invalid
  elif x < a + F(1):
    result = series(a, x, err, est)
    result *= fac; est *= fac           # track absolute error
  else:
    result = conFrac(a, x, err, est)
    result *= fac; est *= result        # relative -> absolute error pre compl
    result = F(1) - result
  if not norm.isNil: norm[] = lnG       # optional return

proc χ²*[F](df, x: F; err: F=F(1e-6), est: var float64=doNotUse): F =
  ## Chi-squared CDF `P(χ²<x)` for `df` degrees of freedom.
  gammaI(0.5*df, 0.5*x, err, est)

proc χ²c*[F](df, x: F; err: F=F(1e-6), est: var float64=doNotUse): F =
  ## Chi-squared ComplementaryCDF/SurvivalFunc `P(χ²>x)` for `df` deg.freedom.
  1.F - χ²(df, x, err, est)

when isMainModule:
  import fpUt; when not declared(assert): import std/[assertions, formatFloat]
  assert almostEqual(gammaI(3.0f32, 1.0f32), 0.080301404, 1e-6)
  assert almostEqual(gammaI(3.0f32, 3.0f32), 0.5768099  , 1e-6)
  assert almostEqual(gammaI(3.0f32, 5.0f32), 0.8753480  , 1e-6)
  # These numerical were checked against Wolfram Alpha online evaluator which
  # seems to have a 2X normalization convention { gammaI(a,inf) = 2!= our 1 }.
  assert almostEqual(gammaI(3.0f64, 1.0f64, 1e-16),0.080301397071394193, 1e-16)
  assert almostEqual(gammaI(3.0f64, 3.0f64, 1e-16), 0.5768099188731565, 1e-16)
  assert almostEqual(gammaI(3.0f64, 5.0f64, 1e-16), 0.87534798051691887, 1e-16)
  when defined bench: # ~20 <ns> (~100 clkCyc)/eval on 5.2GHz i7-1370P Perf-Core
    import std/times; var s = 0f32; let t0 = epochTime()
    for df in 1..100:
      var x = 0.1f32; while x < 10: s += χ²(df.float32, x); x += 0.01
    echo 1e9/100/1000*(epochTime() - t0)," <ns> per eval ",s
