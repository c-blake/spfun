## Module for various special functions of classical probability theory (that
## are inspecific to any particular probabilility density/distribution).

import math

func hoeffding*[F](n: int; t: F, err=F(0), est: ptr F=nil): F =
  ## https://en.wikipedia.org/wiki/Hoeffding%27s_inequality
  if not est.isNil: est[] = F(0)  # XXX should really be epsilon[F]
  F(2) * exp(F(-2) * F(n) * t * t)

when isMainModule:
  assert almostEqual(hoeffding(2, 0.5), 0.7357588823428847)
