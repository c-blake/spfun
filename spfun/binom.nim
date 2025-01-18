## Provide various func/sprocs related to the binomial distribution
from beta        import betaI
from std/math    import sqrt, almostEqual
from spfun/gauss import qtl

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

type
  BinomP* = object                      ## Estimator for binary proportions.
    s, n: int                           # n successes aka n hits, n trials
  BinomPAlgo* =enum Wald,Wilson,Agresti ## Binary proportion est CI Method.

proc nHit*(bp: BinomP): int = bp.s      ## Get number of hits registered.

proc nTry*(bp: BinomP): int = bp.n      ## Get number of trials registered.

proc clear*(bp: var BinomP) = bp.s = 0; bp.n = 0 ## Zero all counters.

proc inc*(bp: var BinomP, x=0, y=1) = inc bp.s, x; inc bp.n, y
  ## Register a trial with a possible hit.  `y > 1` enables a batch interface.

proc est*[F](bp: BinomP; pLo, pHi: var F; ci=0.95, algo=Wilson) =
  ## Return a `ci` confidence interval for population hit ratio `p`.
  case algo # en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
  of Wald:              # Laplace did it first, but it's known by this name
    let z  = -gauss.qtl (1 - ci)*0.5
    let nI = 1.0/bp.n.float
    let p  = bp.s.float*nI
    let s  = sqrt(p*(1-p)*nI)
    pLo = p - z*s
    pHi = p + z*s
  of Wilson:
    let z  = -gauss.qtl (1 - ci)*0.5
    let z2 = z*z
    let s  = bp.s.float
    let n  = bp.n.float
    let c  = (s + 0.5*z2) / (n + z2)
    let r  = z/(n + z2)*sqrt(s*(n - s)/n + 0.25*z2)
    pLo = max(0.0, c - r)
    pHi = min(1.0, c + r)
  of Agresti:
    let z  = -gauss.qtl (1 - ci)*0.5
    let z2 = z*z
    let nT = bp.n.float + z2
    let nI = 1.0/nT
    let pT = (bp.s.float + 0.5*z2)*nI
    let sT = z*sqrt(pT*(1-pT)*nI)
    pLo = max(0.0, pT - sT)
    pHi = min(1.0, pT + sT)

proc initBinomP*(hits, tries: int): BinomP = result.s = hits; result.n = tries

proc est*(bp: BinomP; ci=0.95, algo=Wilson): tuple[pLo, pHi: float] =
  bp.est(result.pLo, result.pHi, ci, algo)

when isMainModule:
  when defined test:                    # Test with 10 coin-flip examples
    when not declared(assert): import std/[formatfloat, assertions]
    assert almostEqual(binom_cdf(0.5, 10, 3), 0.171875)
    assert almostEqual(binom_cdf(0.5, 10, 6), 0.828125) # 1 - 0.171875
    assert binom_qtl(0.5, 10, 0.1718) == 2
    assert binom_qtl(0.5, 10, 0.1719) == 3
    assert binom_qtl(0.5, 10, 0.828) == 5
    assert binom_qtl(0.5, 10, 0.829) == 6
    var e: float64
    echo binom_cdf(0.5f32, 10, 3, 0.01, e), " +- ", e
  else:
    when defined(nimPreviewSlimSystem): import std/formatFloat
    proc p(sn: seq[int]; ci=0.95, algo=Wilson): tuple[pLo, pHi: float] =
      ## Estimate population binary proportion from successes & total trials.
      if sn.len != 2: quit "Need success & totalN", 1
      var bp = BinomP(s: sn[0], n: sn[1])
      bp.est(result.pLo, result.pHi, ci, algo)
    import cligen; dispatch p, cmdName="binom", echoResult=true, help={
      "sn"  : "success & total trials",
      "ci"  : "confidence interval size",
      "algo": "method: Wald, Wilson, Agresti"}
