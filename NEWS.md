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
