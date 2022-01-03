import math

func pdf*[F: SomeFloat](x: F): F =
  ## PDF for unit Gaussian/normal distribution, N(0,1)
  exp(F(-0.5)*x*x)/sqrt(F(2.0*PI))

func cdf*[F: SomeFloat](x: F): F = F(0.5) * (F(1) + erf(x/sqrt(F(2))))
  ## CDF of the unit Gaussian/normal distribution somtimes known as "Phi(x)".

func qtl*[F: SomeFloat](p: F): F =
  ## Acklam's rational approx for inv. Gaussian CDF. |RelErr|<1.15e-9 globally.
  ## Ie., qtl(0.25) to qtl(0.75) gives the interquartile range of a unit normal.
  let a = [F(-3.969683028665376e+01), 2.209460984245205e+02,
             -2.759285104469687e+02, 1.383577518672690e+02,
             -3.066479806614716e+01, 2.506628277459239e+00]
  let b = [F(-5.447609879822406e+01), 1.615858368580409e+02,
           -1.556989798598866e+02, 6.680131188771972e+01,
           -1.328068155288572e+01]
  let c = [F(-7.784894002430293e-03), -3.223964580411365e-01,
           -2.400758277161838e+00, -2.549732539343734e+00,
           4.374664141464968e+00, 2.938163982698783e+00]
  let d = [F(7.784695709041462e-03), 3.224671290700398e-01,
           2.445134137142996e+00, 3.754408661907416e+00]
  var q: F
  if p < 0.0 or p > 1.0: return NaN     # Undefined
  elif p == 0.0: return -Inf            # Negative infinity
  elif p == 1.0: return +Inf            # Positive infinity
  elif p < F(0.02425):                  # Lower Region
    q = sqrt(-(F(2)*ln(p)))
    return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
            ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + F(1))
  elif p > F(0.97575):                  # Upper Region
    q = sqrt(-(F(2)*ln(F(1) - p)))
    return -((((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
             ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + F(1)))
  q = p - 0.5                           # Central Region
  let r = q*q
  return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5])*q /
          (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + F(1))
