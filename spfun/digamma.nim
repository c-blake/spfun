var psis = @[0.0, -0.5772156649015328606]   # psis[i] <- DiGammaψ[i]

proc psi*(n: int): float = ## Array-based DiGamma ψ(psi) for integral arguments.
  let n0 = psis.len         # Just use ψ[i+1] - ψ[i] = 1/i recursion to..
  if n0 < n + 2:            #.. lazily populate an array of cached answers.
    psis.setLen n + 2
    for i in n0 - 1 .. n:
      psis[i + 1] = psis[i] + 1/i.float
  psis[n]
