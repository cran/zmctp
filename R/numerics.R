# Normalisation adapted from the audited research implementation. Parameters
# are scalar; vector arguments are observations, quantiles or probabilities.
.ctp_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
    stop(name, " must be a finite numeric scalar.", call. = FALSE)
}
.ctp_parameters <- function(a, b, gama) {
  .ctp_scalar(a, "a"); .ctp_scalar(b, "b"); .ctp_scalar(gama, "gama")
  if (b < 0 || gama <= max(0, 2 * a))
    stop("Require b >= 0 and gama > max(0, 2*a).", call. = FALSE)
}
.ctp_omega <- function(omega) {
  .ctp_scalar(omega, "omega")
  if (omega < 0 || omega > 1)
    stop("omega must be in [0, 1].", call. = FALSE)
}
.ctp_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
}
.ctp_finite <- function(a, b) {
  if (b == 0 && a <= 0 && a == floor(a)) -a else Inf
}
.ctp_logp0 <- function(a, b, gama) {
  .ctp_parameters(a, b, gama)
  if (a == 0 && b == 0) return(0)
  x <- gama - a
  steps <- max(0, ceiling(20 - min(gama, gama - 2 * a)))
  correction <- 0
  if (steps > 0) for (k in seq_len(steps) - 1) {
    u <- x + k
    correction <- correction + log((gama + k) / u) +
      log((gama - 2 * a + k) / u) - log1p((b / u)^2)
  }
  x <- x + steps
  coefficients <- c(1/12, -1/360, 1/1260, -1/1680, 1/1188,
                    -691/360360, 1/156)
  stirling <- function(z) sum(coefficients * Re((1/z)^seq(1, 13, 2)))
  v <- b / x; r <- a / x
  lp <- correction + (x - 0.5) * log1p(v^2) - 2 * b * atan(v) -
    (x + a - 0.5) * log1p(r) - (x - a - 0.5) * log1p(-r) +
    2 * stirling(x + 1i * b) - stirling(x + a) - stirling(x - a)
  if (!is.finite(lp) || lp > 1e-8)
    stop("Normalising constant could not be evaluated.", call. = FALSE)
  min(lp, 0)
}
.ctp_logpmf <- function(x, a, b, gama) {
  lp0 <- .ctp_logp0(a, b, gama)
  if (!length(x)) return(numeric())
  largest <- max(x)
  if (largest > 1000000)
    stop("Probability recurrence exceeds the 1,000,000-count work limit.", call. = FALSE)
  lp <- numeric(largest + 1L); lp[1L] <- lp0
  if (largest > 0) {
    k <- seq_len(largest) - 1L
    scale <- pmax(abs(a + k), b)
    numerator <- rep(-Inf, length(k)); ok <- scale > 0
    numerator[ok] <- 2 * log(scale[ok]) +
      log(((a + k[ok]) / scale[ok])^2 + (b / scale[ok])^2)
    lp[-1L] <- lp0 + cumsum(numerator - log(gama + k) - log(k + 1))
  }
  lp[x + 1L]
}
.ctp_logadd <- function(x, y) {
  z <- pmax(x, y)
  ans <- z + log1p(exp(-abs(x-y)))
  ans[is.infinite(z) & z < 0] <- -Inf
  ans
}
.ctp_log1mexp <- function(x) {
  out <- rep(NA_real_, length(x))
  i <- which(!is.na(x) & x < -log(2)); j <- which(!is.na(x) & x >= -log(2))
  out[i] <- log1p(-exp(x[i])); out[j] <- log(-expm1(x[j]))
  out
}
.ctp_count_n <- function(n) {
  .ctp_scalar(n, "n")
  if (n < 0 || n != floor(n) || n > .Machine$integer.max)
    stop("n must be a nonnegative integer within the R integer range.", call. = FALSE)
  as.integer(n)
}
