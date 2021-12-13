This is a library of special functions of statistics, mathematical physics, or
other applied mathematics and numerical methods that support their evalution.

Style
-----

Since often only 4 byte IEEE single precision is needed and since this uses both
1/2 the space and, in a SIMD vectorized context, 1/2 the time, these routines
are also usually generic over a floating point type `F` which may be `float32`
or `float64`.  If the language ever supports it this might be extended to
`float16` or `float128`.  `std/rationals.Rational[T]` is not crazy, but would
likely require re-implementing many base transcendental functions (sin,exp,..).

If you read pretty much any material on numerical methods then you will soon
know that almost all methods have run-time performance which can be very
sensitive to error tolerances & input parameters.  Even the simplest example of
how many terms to sum in a power series exhibits this.  Lack of caution in
"nesting" evaluations, e.g. as in numerical root finding or numerical
integration, can really multiply the time cost of wanton precision.

Therefore, most library entry points have a request-reply protocol for error.
I.e., they accept an optional requested error or convergence tolerance and they
reply with an estimate on the achieved error.  If the achieved error is too big
then the method has failed, and the value is not to be trusted naively.

In Nim, this is all modeled by two extra default parameters to generics for the
calls: `err=1e-7, est: var float64=privateDummy`.  To receive the reply estimate
you must declare and pass a float64.  If you do not care about these details
then you can simply not pass the extra parameters and hope for the best.

The requested error is relative since the caller does not know the answer ahead
of time.  The reply estimate is absolute since additive adjustments seem more
common than multiplicative and these do not alter the error. { This may be up
for debate/change in early stages.  Perhaps relative error returns are best.
E.g., since a caller need only multiply to convert, not divide. }

Scope
-----

It'd be good to add hypergeometric functions (which often subsume or are used by
so many others), spherical harmonics, bessel functions, elliptic, exponential,
and trigonometric integrals, as well as the density/cdf/quantile functions for
the more common probability distributions such as uniform, exponential, normal,
Chi-square, F, logistic, gamma & beta, Cauchy, Weibull (maybe GEV), Poisson,
Kolmogorov-Smirnov, maybe log Normal.  PRs are very welcome.
