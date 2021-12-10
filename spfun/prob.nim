## Module for various special functions of classical probability theory (that
## are inspecific to any particular probabilility density/distribution).

import math, fpUt

func hoeffding*[F](n: int; t: F, err=F(0), estp: ptr F=nil): F =
  ## https://en.wikipedia.org/wiki/Hoeffding%27s_inequality
  result = F(2) * exp(F(-2) * F(n) * t * t)
  safeSet estp, F(0) * result           # XXX should * epsilon[F]-ish

when isMainModule:
  assert almostEqual(hoeffding(2, 0.5), 0.7357588823428847)
