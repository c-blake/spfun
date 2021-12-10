## This module contains various floating point utility funcs/procs useful in
## either evaluating or testing special functions.

from math import copySign

func almostEqual*[F](x, y, thresh: F): bool = abs(x - y) <= thresh * abs(x + y)
  ## True iff `|x-y| <= thresh * |x+y|`.  More intuitive than stdlib interface.

func safeDivisor*[F](x, eps: F): F {.inline.} =
  ## Replace `x` with `eps` if `|x| < eps`.
  if abs(x) > eps: x else: copySign(x, eps)

template safeSet*(p, v) =
  ## Assign through pointer `p` the value `v`, but do nothing is `p.isNil`.
  if not p.isNil: p[] = v
