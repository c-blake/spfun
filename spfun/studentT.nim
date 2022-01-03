## Provide various func/sprocs related to the Student's T distribution

import beta, math

proc cdf*[F](df, t: F): F =
  ## CDF of the student's T statistic of `df` degrees of freedom
  return betaI(abs(df / (df + t * t)), F(0.5) * df, F(0.5), err=1e-7)

proc corrP*[F](r: F; n: int): F =
  ## P(Pearson linear corr coef > r|null hypothesis of linear independence).
  if abs(r) >= 1 or n < 2:
    F(0)
  else:
    studentT.cdf(F(n - 2), r * sqrt(F(n - 2) / ((F(1) - r) * (F(1) + r))))
