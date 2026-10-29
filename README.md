# zmctp

Distribution functions and likelihood fitting for the complex triparametric
Pearson parent and its zero-inflated mixture.

```r
library(zmctp)
dctp(0:5, a=1, b=0.5, gama=6)
x <- rep(0:4, c(40,12,6,2,1))
fit <- zictp.fit(x)
fit$estimates
fit$diagnostics
plot(fit)
```

Version 0.1.4 replaces support-dependent normalisation and capped quantiles,
profiles the mixing probability including zero, and retains transformed
optimiser coordinates. Public names, `gama`, fitted classes and S3 methods
are retained. See NEWS and `vignette("introduction", package="zmctp")` for
behaviour changes, boundary diagnostics and runtime limits.

The implemented mixture adds zeros for positive omega. It does not implement
zero deflation. Bootstrap inference remains in the research supplement.
Published application datasets are not included in this release candidate.

Maintainer: Rasheedat Oladoja <roladoja@ttu.edu>.
Developers: Rasheedat Oladoja and Kabir Opeyemi Olorede.
License: GPL-3. Original source notices and authorship are retained.
