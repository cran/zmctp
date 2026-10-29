# zmctp 0.1.4

* Declare the stats::setNames and utils::head/tail namespace imports, and
  clarify the DESCRIPTION wording following CRAN's incoming pre-tests.
* Replace requested-support normalisation with a shifted log-gamma normalising
  constant and log recurrence. Density/CDF calls no longer depend on the largest
  count passed in the same call.
* Expand quantile searches beyond `max_q` (now an initial endpoint). Add optional
  `max_count`; exceeding the work limit signals an error instead of capping a
  count. Quantiles at probability one use the distribution support endpoint.
* Support omega zero and one in distribution functions. Correct the explanation
  of zero inflation; positive omega increases zero probability for the same
  nondegenerate parent.
* Fit with multiple starts, retain transformed coordinates, and profile omega
  analytically including the exact parent boundary. Preserve method defaults;
  requested methods supplement Nelder-Mead and nlminb. Retain omega_start as a
  validated compatibility argument; profiling makes its value unnecessary.
* Preserve fitted classes, existing fields and S3 methods. Add theta, attempts
  and diagnostics. raw_optim$par now uses a/log_b/log_margin coordinates.
* Reject all-zero fits as unidentified. Validate probability and moment domains.
  Return zero probabilities for invalid support points and propagate missing
  vector entries. Random-generation n must be an integer, rather than truncated.
* Withhold Wald uncertainty near boundaries or for an unsuitable Hessian; remove
  Hessian jitter and negative covariance clipping.
* Include unobserved cells and an upper-tail bin in fitted_freq. Clarify that
  Pearson and historically named Wald discrepancies are descriptive.
* Add numerical regression tests and replace static vignette results with
  executable constructed examples. Add Kabir Opeyemi Olorede as a package author;
  Rasheedat Oladoja remains the registered maintainer.
* Keep bootstrap inference in the research supplement for this release.

# zmctp 0.1.3

## Bug fixes

* Fixed a starting-value inconsistency in `ctp.fit()` that could cause
  `optim()` to silently converge to a degenerate local optimum for
  moderately overdispersed data. `a_start` and `gama_start` are derived
  independently from the sample mean/variance, with nothing guaranteeing
  `gama_start > 2*a_start + 2`. When that margin was non-positive,
  `eta_start` was clamped to `log(1e-6)`, pinning the optimizer's start
  right on the variance-existence singularity (`gama = 2a + 2`), where
  the log-density has pathological gradients. `ctp.fit()` now falls back
  to a safe starting margin (with a `warning()`) instead of clamping to
  the boundary, and additionally tries a small grid of alternative
  starting points, keeping whichever converges to the best
  log-likelihood.

* Applied the same safe-margin fallback to the starting-value grid in
  `zictp.fit()`, which shared the same underlying clamp pattern in each
  of its multi-start candidates.

* `ctp.fit()` and `zictp.fit()` now flag fits whose estimated `gama` is
  very close to the variance-existence boundary (`gama - 2*a - 2` near
  zero) with a `warning()`, a `near_boundary` element in the returned
  object, and a note in `print()`/`summary()` output -- instead of
  silently reporting `Convergence: YES` on an unreliable fit.

* Fixed a severe performance bug in `dctp()`: its normalization-tail loop
  grew a vector one element at a time via `c()`, which is O(n^2) in R.
  Near the `gama = 2*a + 2` boundary the tail decays very slowly and can
  need close to the full 10000-term cap, making individual `dctp()` calls
  take seconds rather than milliseconds -- and this compounds badly
  across an optimizer's many function evaluations during `ctp.fit()`/
  `zictp.fit()`, especially for the exact overdispersed, near-boundary
  data most likely to trigger the starting-value bug above. The tail
  buffer is now preallocated; a `dctp()` call near the boundary that
  previously took ~0.15s now takes ~0.002s, with identical numerical
  output.

# zmctp 0.1.2

## Bug fixes

* Fixed an incorrect mode formula in `mode_ctp()` and `mode_zictp()`.
  The previous implementation used a naive "round to nearest integer"
  rule that could return the wrong mode (e.g., returning 1 when the
  true mode was 0). Both functions now use the correct threshold-based
  rule.

* Fixed incorrect variance and skewness formulas in `var_ctp()` and
  `skew_ctp()`. The previous implementation used an incorrect
  algebraic combination of terms. `var_zictp()` and `skew_zictp()`
  required no changes and inherit the fix automatically, since they
  are built directly on `var_ctp()` and `skew_ctp()`.

* Added an explicit validity check (`gama > 0`) to `ctp.fit()` and
  `zictp.fit()`. This closes a narrow edge case in the internal
  parameter reparameterization where sufficiently extreme starting
  conditions could otherwise allow an invalid (non-positive) `gama`
  during optimization. Typical fits to well-behaved data are
  unaffected.

## Other changes

* `tests/testthat/test-moments.R` now contains a full suite of
  regression tests for `mean_ctp()`, `var_ctp()`, `skew_ctp()`,
  `mode_ctp()`, and their zero-modified counterparts, including a
  Monte Carlo cross-check of the skewness formula against simulated
  data from `rzictp()`.

## Note for existing users

If you used `var_ctp()`, `skew_ctp()`, `var_zictp()`, `skew_zictp()`,
`mode_ctp()`, or `mode_zictp()` in version 0.1.1 or earlier, results
from those functions will differ under 0.1.2 -- the new results are
correct. We recommend re-running any analysis that relied on their
numeric output.

# zmctp 0.1.1

* Added exported functions for theoretical moments and mode of the CTP
  and Zero-Modified CTP distributions: `mean_ctp()`, `var_ctp()`,
  `skew_ctp()`, `mode_ctp()`, `mean_zictp()`, `var_zictp()`,
  `skew_zictp()`, `mode_zictp()`. These were previously computed only
  internally within `summary.ctpfit()` and `summary.zictpfit()`.

# zmctp 0.1.0

* Initial CRAN release.
