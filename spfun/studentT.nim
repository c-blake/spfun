## Provide various func/sprocs related to the Student's T distro & correlations.

from beta          import betaI
from std/math      import sqrt, sum, fac
from std/algorithm import sort, nextPermutation
from std/random    import shuffle, randomize
from std/strutils  import strip, parseFloat

proc cdf*[F](df, t: F): F =
  ## CDF of the student's T statistic of `df` degrees of freedom
  return betaI(abs(df / (df + t * t)), F(0.5) * df, F(0.5), err=F(1e-7))

proc corrP*[F](r: F; n: int): F =
  ## P(Pearson linear corr coef > r|null hypothesis of linear independence).
  if abs(r) >= 1 or n < 2:
    F(0)
  else:
    studentT.cdf(F(n - 2), r * sqrt(F(n - 2) / ((F(1) - r) * (F(1) + r))))

## Studentized permutation test of Pearson & Spearman Correlation Coefs.  Core
## idea is that under H0=no relationship, relative order of elements makes no
## difference => sampling over stat(shuffles) yields P(stat|H0).  Yu shows test
## better controls false positives for association under general scenarios like
## small sample size & dependent-but-uncorrelated variables.  Currently, test
## only supports H0=0, but alternative hypothesis can be 1|2-sided.  Optimized
## version of https://github.com/hyu-ub/perk by Han Yu&Alan Hutson.  This code
## lives here since it & CLI util directly check asymptotic approximation of
## Student's T (ultimately marketed by Kendall & Stuart 1979).

proc toRanks*[F](xs: openArray[F]): seq[F] =
    ## Convert xs[i] to 1-origin ranks, mid-ranking exact ties.
    if xs.len < 1: return
    var xi = newSeq[(F, int)](xs.len)
    for i, x in xs: xi[i] = (x, i)
    xi.sort                             # Create a list of sorted (x,i) pairs
    var i = 0                           # Move indices (i,j) over them where
    while i < xs.len:                   #..j != i only if there are exact ties.
        var j = i + 1
        while j < xs.len:               # Move j to 1-past any block of ties
            if xi[j][0] != xi[i][0]:    # Do float-rounding insensitive compare?
                break
            j += 1
        for k in i..<j:                 # Mid-rank for tie Else .5*(i+i+1+1)=i+1
            xi[k][0] = 0.5*(F(i) + F(j) + 1.0)  # + 0.0 for 0-origin ranks
        i = j
    result.setLen xs.len                # Allocate result list
    for (x, i) in xi: result[i] = x     # Copy just ranks to result

proc ssq[F](xs: openArray[F]): F = (for x in xs: result += x*x) # sum of squares

proc center[F](xs: var openArray[F]) =  # Center about the mean
  let mu = xs.sum/xs.len.F
  for x in mitems(xs): x -= mu

proc scale[F](xs: var openArray[F]) =   # Scale to unit variance
  let vI = sqrt(F(xs.len.F) / xs.ssq)
  for x in mitems(xs): x *= vI

proc cor[F](xs, ys: openArray[F]): F =  # Pearson Corr Coef of Already Unitized
  for i in 0 ..< xs.len: result += xs[i]*ys[i]
  result /= xs.len.F

proc tau[F](xs, ys: openArray[F]): F =  # Constant for Studentization
  for i in 0 ..< xs.len: result += xs[i]*xs[i] * ys[i]*ys[i]
  sqrt result/xs.len.F

template forPerm[F](xs: var openArray[F]; B: int; body) =
  let nF = if xs.len < 12: xs.len.fac else: 12.fac
  let Bp = min(B, nF)                   # B' = either exactly n! or user value
  if Bp < nF:                           # Randomly sample likely unique perms
    for b in 1..Bp: xs.shuffle; body
# When B << n!, above repeats averaging work & creates unneeded sampling noise.
  else: # So, in that case, do whole set of permutations (must repeat `body`).
    xs.sort                             # Nim perm gen is *lexicographic* & with
    block: body                         #..dup vals needs in-order start to work
    while xs.nextPermutation: (block: body)

type Corr* = enum linear, rank                  ## Linear | Rank correlation
type AltH* = enum less="-", greater="+", twoSide, form ## Kind of altern.hypoth.

proc ccPv*[F](xs, ys: openArray[F]; B=999, cc=linear, altH=twoSide):
       tuple[cc, pVal: F] =
  ## Pearson Linear|Spearman Rank Correlation Coefficient with p-Value options
  var xs: seq[F] = if cc == rank: xs.toRanks else: xs[0..^1]
  var ys: seq[F] = if cc == rank: ys.toRanks else: ys[0..^1]
  xs.center; xs.scale; ys.center; ys.scale
  result.cc = cor(xs, ys)
  if altH == form: result.pVal = corrP(result.cc, xs.len); return
  let rS0 = result.cc / tau(xs, ys)     # Studentized to compare to
  var n, nB: int                        # count of true side-conditions
  forPerm(xs, B):                       # Sample(PERM) preserves mean=0, sd=1;
    let rSS = cor(xs, ys) / tau(xs, ys) #..So, can standardize xs & ys ONCE..
    case altH                           #..optimizing the rS sampling process.
      of twoSide: inc n, (rSS.abs > rS0.abs).int
      of less:    inc n, (rSS < rS0.abs).int
      of greater: inc n, (rSS > rS0).int
      else: discard
    inc nB      # Dups may make saturated case have < n! arrangements of xs[i].
  result.pVal = n.F / nB.F

when isMainModule:
  when defined danger: randomize()
  proc toN(p:string):seq[float]=(for x in p.lines:result.add x.strip.parseFloat)
  proc p(xy:seq[string]; B=999,corr=linear,altH=twoSide): tuple[cc,pVal:float]=
    ## Pearson Linear|Spearman Rank Correlation Coefficient with p-Value options
    let x = xy[0].toN; let y = xy[1].toN
    if y.len != x.len: quit "x.len != y.len\n", 1
    ccPv(x, y, B, corr, altH)
  import cligen; dispatch p, cmdName="percc", echoResult=true, help={
    "xy"  : "paths to 2 ASCII number-vector files",
    "B"   : "permutations to sample",
    "corr": "CorrCoeff to test: linear, rank",
    "altH": "Altern. Hypothesis: - + twoSide form"}
