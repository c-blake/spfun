var psis = @[0.0, -0.5772156649015328606]   # psis[i] <- DiGammaψ[i]

proc psi*(n: int): float =  # Just use ψ[i+1] - ψ[i] = 1/i recursion to..
  ## Array-based DiGamma ψ(psi) function for integral arguments.
  let n0 = psis.len     #.. lazily populate an array of cached answers.
  if n0 < n + 2:
    psis.setLen n + 2
    for i in n0 - 1 .. n:
      psis[i + 1] = psis[i] + 1/i.float
  psis[n]
