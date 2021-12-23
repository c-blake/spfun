## Module for various special functions of classical probability theory (that
## are inspecific to any particular probabilility density/distribution).

import math, fpUt

var doNotUse: float64 # Only for default value; Callers must provide `est`.

func hoeffding*[F](n: int; t: F, err=F(0), est: var float64=doNotUse): F =
  ## https://en.wikipedia.org/wiki/Hoeffding%27s_inequality
  result = F(2) * exp(F(-2) * F(n) * t * t)
  est = 0.0f64                          # XXX should * epsilon[F]-ish * result

when isMainModule:
  assert almostEqual(hoeffding(2, 0.5), 0.7357588823428847)
