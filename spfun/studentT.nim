## Provide various func/sprocs related to the Student's T distro & correlations.

from beta          import betaI
from binom         import BinomP, inc, est, nHit, nTry
from std/math      import sqrt, sum, fac
from std/algorithm import sort
from std/random    import shuffle, randomize
from std/strutils  import strip, parseFloat
when defined(nimPreviewSlimSystem): import std/syncio

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
## difference => stat(shuffle) sampling gives P(stat|H0).  Yu&Hutson2020 shows
## test better controls falsePos for association under general scenarios like
## small sample size & dependent-but-uncorrelated variables.  Currently, test
## supports only H0=0, but alternative hypoth can be 1|2-sided.  This is like
## github.com/hyu-ub/perk but optimized.  Code is here since it & CLI directly
## check asymptotic approximation of Student's T (marketed by Kendall 1979),
## although there's also a strong argument to be made for `spfun/binom.nim`.
##
## When nSamp<n! this repeats cor work (more so w/ties) & adds sampling noise
## relative to `nextPermutation`.  OTOH, full averages: 1) don't easily fit into
## target pValue-optimized approaches, 2) raise questions about re-sampling from
## data vs. smoothed data, 3) Drawing statistical conclusions from <6 points is
## too rarely sound to optimize for, 4) 6!=720 is often "good enough" anyway.

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

proc tau[F](xs, ys: openArray[F]): F =  # Constant for Studentization
  for i in 0 ..< xs.len: result += xs[i]*xs[i] * ys[i]*ys[i]
  result.sqrt

type Corr* = enum linear, rank                  ## Linear | Rank correlation
type AltH* = enum less="-", greater="+", twoSide, form ## Kind of altern.hypoth.
var nEvalCcPv = 0
proc ccPvEvalCount*(): int = nEvalCcPv ## Global time evaluation counter

proc ccPv*[F](xs,ys: openArray[F]; minTry=9, maxTry=5000, ci=0.95, pVthr=0.05,
              cc=linear, altH=twoSide, verbose=false): tuple[cc, pLo, pHi: F] =
  ## Pearson Linear|Spearman Rank Correlation Coefficient with p-Value options
  if xs.len != ys.len or xs.len == 0: return               # Shuffle preserves..
  var xs: seq[F] = if cc==rank: xs.toRanks else: xs[0..^1] #..mean=0,sd=1 =>can
  var ys: seq[F] = if cc==rank: ys.toRanks else: ys[0..^1] #..standardize xs&ys
  xs.center; xs.scale; ys.center; ys.scale                 #..ONCE to optimize
  let nIRt = 1.0/sqrt(xs.len.F)                            #..sample re-calc to
  result.cc = cor(xs, ys)/xs.len.F                         #..just dot&dot2.rt.
  if altH == form:
    result.pLo = corrP(result.cc, xs.len); result.pHi = result.pLo; return
  let rS0 = result.cc/(nIRt*tau(xs,ys)) # Studentized to compare to
  var bp: BinomP                        # Count side-conditions that are true
  template count[F](bp: var BinomP; xs, ys: openArray[F]; rS0: F; altH: AltH) =
    let rSS = cor(xs, ys)*nIRt/tau(xs, ys)
    case altH
      of twoSide: bp.inc (rSS.abs > rS0.abs).int
      of less:    bp.inc (rSS < rS0.abs).int
      of greater: bp.inc (rSS > rS0).int
      else: discard
  let nE0 = nEvalCcPv
  for it in 1..minTry: xs.shuffle; bp.count(xs, ys, rS0, altH)
  while(bp.est(result.pLo,result.pHi,ci);result.pHi>pVthr and result.pLo<pVthr):
    xs.shuffle                          # Repeat while not sure (@`ci`)..
    bp.count(xs, ys, rS0, altH)         #..if pV is lower or higher.
    inc nEvalCcPv                       # MT-unsafe BUT race only loses counts
    if nEvalCcPv > nE0 + maxTry: break
  if verbose: stderr.write "ccPv: ",bp.nHit," / ",bp.nTry," : ",bp.est,"\n"

when isMainModule:
  when defined(nimPreviewSlimSystem): import std/formatFloat
  when defined danger: randomize()
  proc toN(p:string):seq[float]=(for x in p.lines:result.add x.strip.parseFloat)
  proc p(xy: seq[string]; minTry=9, maxTry=9, ci=0.95, pVthr=0.05, corr=linear,
         altH=twoSide, verbose=false): tuple[cc, pLo, pHi: float] =
    ## Pearson Linear|Spearman Rank Correlation Coefficient with p-Value options
    let x = xy[0].toN; let y = xy[1].toN
    if y.len != x.len: quit "x.len != y.len\n", 1
    ccPv x, y, minTry, maxTry, ci, pVthr, corr, altH, verbose
  import cligen; dispatch p, cmdName="studentT", echoResult=true, help={
    "xy"     : "paths to 2 ASCII number-vector files",
    "minTry" : "min permutations to sample",
    "maxTry" : "max permutations to sample",
    "ci"     : "conf.ival on pValue",
    "pVthr"  : "pValue threshold",
    "corr"   : "CorrCoeff to test: linear, rank",
    "altH"   : "Altern. Hypothesis: - + twoSide form",
    "verbose": "Altern. Hypothesis: - + twoSide form"}
