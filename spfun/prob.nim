## Module for various special functions of classical probability theory (that
## are inspecific to any particular probabilility density/distribution).

import math

proc hoeffding*[T](n: int; t: T): T = T(2) * exp(T(-2) * T(n) * t * t)
  ## https://en.wikipedia.org/wiki/Hoeffding%27s_inequality
