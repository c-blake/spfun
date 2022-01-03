from math import PI, arctan, tan

proc pdf*[F: SomeFloat](x: F): F =
  ## PDF for unit Cauchy distribution.
  F(1)/F(PI)/(F(1) + x*x)

proc cdf*[F: SomeFloat](x: F): F =
  ## CDF for unit Cauchy distribution.
  arctan(x)/F(PI) + F(0.5)

proc qtl*[F: SomeFloat](p: F): F =
  ## Inverse CDF aka quantile for unit Cauchy distribution.
  tan(F(PI)*(p - 0.5))
