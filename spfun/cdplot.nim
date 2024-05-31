when not declared(stdin): import std/[syncio, formatfloat]
import std/[strutils, strformat, algorithm]
from spfun/binom     import initBinomP, est
from cligen/osUt     import mkdirOpen
from cligen/colorScl import rgb, hex
type ConfBand = enum pointWise, simultaneous
type Floats = seq[float] #^en.wikipedia.org/wiki/CDF-based_nonparametric_confidence_interval

iterator edf*[T](ts: seq[T]): (T, int) =
  ## (x, CumCnt(x)) unnormalized EDF with ties of already sorted seq `ts`; Unix
  ## uniq -c-esque but cumulative.  E.g.: 1 2 2 2 3 4 4 -> 1,1 2,4 3,5 4,7.
  var tLast: T
  var c: int
  for t in ts:
    if c > 0 and t != tLast:    # Not first point
      yield (tLast, c)
    tLast = t
    inc c
  if c > 0:                     # Could be empty
    yield (tLast, c)

proc cdplot*(tp="/tmp/cd/", gplot="", nCI=20, band=pointWise, wvls: Floats= @[],
             vals: Floats = @[], alphas: Floats = @[], xlabel="Sample Value",
             inputs: seq[string]) =
  ## Generate files & gnuplot script to render CDFs with confidence band blurs.
  ## If `.len<inputs.len` the final value of `wvls`, `vals`, or `alphas` is
  ## re-used for subsequent inputs, otherwise they match pair-wise.
  let inputs = if inputs.len > 0: inputs else: @[""]
  template setup(id, arg, default) =    # Ensure ok (wvls|vals|alphas)[i] ..
    var id = arg                        #.. for each inputs[i].
    if id.len == 0: id.add default
    for i in 1 .. inputs.len - id.len: id.add id[^1]
  setup wvls, wvls, 0.87; setup vals, vals, 1.0; setup alphas, alphas, 0.5
  for p in inputs:
    var xs: seq[float]
    for f in lines(if p.len>0: p.open else: stdin): xs.add f.strip.parseFloat
    xs.sort; let n = xs.len             # Make, then emit EDF, C.Band files
    let e = mkdirOpen(&"{tp}/{p}EDF", fmWrite)
    var cifs: seq[File]
    for j in 1..nCI:
      cifs.add open(&"{tp}/{p}ci-{j}", fmWrite)
      cifs.add open(&"{tp}/{p}ci+{j}", fmWrite)
    for (f, c) in edf(xs):
      e.write f," ",c.float/n.float,"\n"
      for j in 1..nCI:
        let ci = j.float/float(nCI + 1)
        let (lo, hi) = initBinomP(c, n).est(ci)
        cifs[2*j-2].write f," ",lo,"\n"
        cifs[2*j-1].write f," ",hi,"\n"
    e.close; for cif in cifs: cif.close
  let g = if gplot.len > 0: open(gplot, fmWrite) else: stdout
  g.write &"""#!/usr/bin/gnuplot
# set terminal png size 1920,1080 font "Helvetica,10"; set output "cbands.png"
set key top left noautotitle    # EDFs go bot left->up right;Dot keys crowd plot
set style data steps
set xlabel "{xlabel}"
set ylabel "Probability"
plot """
  for i, p in inputs:
    for j in 1..nCI:
      let cCI = rgb(wvls[i], sat=1.0 - j.float/float(nCI + 1), val=vals[i]).hex
      let alph = if alphas[i] < 1.0: int(alphas[i]*256).toHex[^2..^1] else: ""
      let s = if i==0 and j==1: "" else: ",\\\n     "
      g.write s, &"'{tp}/{p}ci-{j}' lw 2 lc rgb '#{alph}{cCI}'"
      g.write &",\\\n     '{tp}/{p}ci+{j}' lw 2 lc rgb '#{alph}{cCI}'"
    let cEDF = rgb(wvls[i], sat=1.0, val=vals[i]).hex
    let lab = if p.len>0: p else: "stdin"
    g.write &",\\\n     '{tp}/{p}EDF'  lw 3 lc rgb '#{cEDF}' t '{lab}'"
  g.write "\n"; g.close

when isMainModule:
  import cligen; include cligen/mergeCfgEnv; dispatch cdplot, help={
    "inputs": "input paths or \"\" for stdin",
    "tp"    : "temp path prefix for numbered CI levels",
    "gplot" : "gnuplot script or \"\" for stdout",
    "nCI"   : "number CI levels for conf.band shading",
    "band"  : "pointWise|simultaneous confidence bands",
    "wvls"  : "cligen/colorScl HSV-based wvlens; 0.6",
    "vals"  : "values (V) of HSV fame; 0.8",
    "alphas": "alpha channel transparencies; 0.5"}
