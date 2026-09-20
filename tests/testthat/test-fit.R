# =============================================================================
# tests/testthat/test-fit.R
# Regression tests for the ctp.fit()/zictp.fit() starting-value bug: a_start
# and gama_start are derived independently from the sample mean/variance, so
# nothing guarantees gama_start > 2*a_start + 2. Previously, when that margin
# was <= 0, eta_start was silently clamped to log(1e-6), pinning the optimizer
# right on the variance-existence singularity (gama = 2a + 2).
# =============================================================================

test_that("ctp.fit() warns instead of silently starting on the singularity", {
  set.seed(1)
  x <- rctp(200, a = 1, b = 0.5, gama = 8)

  # a_start = 5 forces gama_start (computed from x) - 2*5 - 2 to be very
  # negative, reproducing the inconsistent-heuristic trigger condition.
  expect_warning(
    fit <- ctp.fit(x, a_start = 5),
    "Default gama_start is inconsistent with a_start"
  )
  expect_true(is.finite(fit$logLik))
})

test_that("ctp.fit() flags a fit whose gama lands on the variance-existence boundary", {
  set.seed(42)
  x <- rnbinom(1550, mu = 5, size = 0.5)

  fit <- suppressWarnings(ctp.fit(x))
  expect_true(fit$converged)
  margin <- unname(fit$estimates["gama"] - 2 * fit$estimates["a"] - 2)
  expect_true(margin >= 0)
  if (margin < 1e-3) {
    expect_true(fit$near_boundary)
  }
})

test_that("ctp.fit() finds a fit at least as good as the truth on well-behaved data", {
  # See the analogous zictp.fit() test below for why this checks logLik
  # against the truth rather than exact parameter recovery.
  set.seed(7)
  x <- rctp(500, a = 1, b = 0.5, gama = 8)

  fit <- expect_silent(ctp.fit(x))
  expect_true(fit$converged)
  expect_false(fit$near_boundary)

  ll_true <- sum(dctp(x, 1, 0.5, 8, log = TRUE))
  expect_true(fit$logLik >= ll_true - 1e-2)
})

test_that("zictp.fit() flags a fit whose gama lands on the variance-existence boundary", {
  set.seed(42)
  x <- rnbinom(1550, mu = 5, size = 0.5)
  x[sample(seq_along(x), 100)] <- 0

  fit <- suppressWarnings(zictp.fit(x))
  expect_true(fit$converged)
  margin <- unname(fit$estimates["gama"] - 2 * fit$estimates["a"] - 2)
  expect_true(margin >= 0)
  if (margin < 1e-3) {
    expect_true(fit$near_boundary)
  }
})

test_that("zictp.fit() finds a fit at least as good as the truth on well-behaved data", {
  # NOTE: exact parameter recovery is not a fair test of an MLE routine here --
  # the ZM-CTP model has known near-flat-likelihood/identifiability regions at
  # high zero-inflation (omega and a small mean can trade off), so a stochastic
  # sample can legitimately have its MLE sit far from the generating parameters.
  # The correct sanity check is that the fit does at least as well, in
  # likelihood, as the true parameters -- not that it recovers them exactly.
  set.seed(7)
  x <- rzictp(1000, a = 2, b = 1, gama = 10, omega = 0.2)

  fit <- zictp.fit(x)
  expect_true(fit$converged)
  expect_false(fit$near_boundary)

  ll_true <- sum(dzictp(x, 2, 1, 10, 0.2, log = TRUE))
  expect_true(fit$logLik >= ll_true - 1e-2)
})
