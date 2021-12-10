## Provide various func/sprocs related to the Student's T distribution

import beta

func T_cdf(df, t: float): float =
  ## CDF of the student's T statistic of `df` degrees of freedom
  return betaI(abs(df / (df + t * t)), 0.5 * df, 0.5, err=1e-7)
