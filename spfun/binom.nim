## Provide various func/sprocs related to the binomial distribution
import beta

var doNotUse: float64 # Only for default value; Callers must provide `est`.

func binom_cdf*[F](p: F; n,k: int, err=1e-7,
                   est: var float64=doNotUse): F {.inline.} =
  ## CDF of a (p,n) binomial random variable
  F(1) - betaI(p, F(k + 1), F(n - k), err, est)

proc binom_qtl*[F](p: F, n: int, qtl: F): int =
  ## Least `k` such that `binom_cdf(p,n,k) <= qtl` (by binary search).
  let e = F(0.01) / n.F                 # only compute betaI to needed err
  var lo = 0
  var hi = n
  while hi > lo + 1:
    let k = (hi + lo) shr 1
    if binom_cdf(p, n, k, e) > qtl: hi = k  # cdf@midPt > target; lower hi
    else: lo = k                            # cdf@midPt <= target; raise lo
  lo

when isMainModule:
  import math                           # test with a 10 coin-flip example
  assert almostEqual(binom_cdf(0.5, 10, 3), 0.171875)
  assert almostEqual(binom_cdf(0.5, 10, 6), 0.828125) # 1 - 0.171875
  assert binom_qtl(0.5, 10, 0.1718) == 2
  assert binom_qtl(0.5, 10, 0.1719) == 3
  assert binom_qtl(0.5, 10, 0.828) == 5
  assert binom_qtl(0.5, 10, 0.829) == 6
  var est: float64
  echo binom_cdf(0.5f32, 10, 3, 0.01, est), " +- ", est
