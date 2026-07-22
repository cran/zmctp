# =============================================================================
# tests/testthat/test-moments.R
# =============================================================================

test_that("mean_ctp, var_ctp, skew_ctp match known hand-verified values", {
  # a=1, b=0.5, gama=6 -- independently verified via direct numerical
  # summation of the pmf (not just the closed-form formula), see paper
  # Section 3.2 for the full derivation
  expect_equal(mean_ctp(1, 0.5, 6), 0.41667, tolerance = 1e-4)
  expect_equal(var_ctp(1, 0.5, 6), 1.12847, tolerance = 1e-4)
  expect_equal(skew_ctp(1, 0.5, 6), 8.15843, tolerance = 1e-3)
})

test_that("mean_zictp, var_zictp, skew_zictp match known hand-verified values", {
  # Same base parameters, omega = 0.3
  expect_equal(mean_zictp(1, 0.5, 6, 0.3), 0.29167, tolerance = 1e-4)
  expect_equal(var_zictp(1, 0.5, 6, 0.3), 0.82639, tolerance = 1e-4)
  expect_equal(skew_zictp(1, 0.5, 6, 0.3), 9.49929, tolerance = 1e-3)
})

test_that("mean_zictp / var_zictp / skew_zictp reduce correctly to CTP at omega = 0", {
  a <- 2; b <- 1; gama <- 10
  expect_equal(mean_zictp(a, b, gama, 0), mean_ctp(a, b, gama))
  expect_equal(var_zictp(a, b, gama, 0), var_ctp(a, b, gama))
  expect_equal(skew_zictp(a, b, gama, 0), skew_ctp(a, b, gama))
})

test_that("mode_ctp gives the correct mode, verified against direct pmf computation", {
  # These five cases were independently verified by directly computing the
  # true pmf via forward recurrence and confirming the argmax (see paper
  # Section 3.2 for the full derivation and numerical verification)
  
  # Over-dispersion regime: true mode is 0 (naive rounding of the old,
  # incorrect formula gave 1 here -- this is the regression test for
  # the mode formula bug fix)
  expect_equal(mode_ctp(4.0, 0.5, 20), 0L)
  
  # Under-dispersion regime: true mode is 1
  expect_equal(mode_ctp(-1.5, 1.8, 3.8), 1L)
  
  # Vaccine adverse events fitted parameters: true mode is 0
  expect_equal(mode_ctp(4.084763, 0.0004433845, 20.2193974), 0L)
  
  # UK coal mining strikes fitted parameters: true mode is 1
  expect_equal(mode_ctp(-1.542897, 1.8604172351, 3.7795715), 1L)
  
  # bioChemists fitted parameters: true mode is 0
  expect_equal(mode_ctp(2.866173, 2.4829756539, 15.2339793), 0L)
})

test_that("mode_ctp correctly identifies the exact tie (bimodal) case", {
  a <- 2; b <- 1; y0_target <- 3
  gama <- (a^2 + b^2 + y0_target*(2*a - 1)) / (y0_target + 1)
  
  result <- mode_ctp(a, b, gama)
  expect_equal(length(result), 2)
  expect_equal(result, c(y0_target, y0_target + 1))
})

test_that("mode_zictp reduces to mode_ctp at omega = 0", {
  a <- -1.5; b <- 1.8; gama <- 3.8
  expect_equal(mode_zictp(a, b, gama, 0), mode_ctp(a, b, gama))
})

test_that("mode_zictp shifts to 0 as omega increases sufficiently", {
  # Under-dispersion regime: parent mode is 1; mode should shift to 0
  # once omega crosses the threshold derived in the paper's Proposition
  # on the mode of ZI-CTP
  a <- -1.5; b <- 1.8; gama <- 3.8
  expect_equal(mode_zictp(a, b, gama, 0.0), 1L)
  expect_equal(mode_zictp(a, b, gama, 0.1), 1L)
  expect_equal(mode_zictp(a, b, gama, 0.3), 0L)
  expect_equal(mode_zictp(a, b, gama, 0.5), 0L)
})

test_that("skew_zictp matches empirical skewness from simulated data", {
  skip_on_cran()  # slower Monte Carlo check, not run during CRAN submission checks
  
  # NOTE: deliberately using gama=15 here rather than gama=6. At gama=6 the
  # skewness-existence margin (gama > 2*a+3) is only 1, giving a slowly-
  # decaying tail that requires an impractically large simulated sample
  # (confirmed during development to need >50,000 draws and even then
  # converges slowly) to match the theoretical skewness within any
  # reasonable tolerance. gama=15 gives a comfortable margin and a stable,
  # fast-converging empirical check.
  set.seed(2026)
  a <- 1; b <- 0.5; gama <- 15; omega <- 0.3
  x <- rzictp(50000, a, b, gama, omega)
  theoretical <- skew_zictp(a, b, gama, omega)
  empirical <- mean((x - mean(x))^3) / sd(x)^3
  expect_equal(theoretical, empirical, tolerance = 0.15)
})