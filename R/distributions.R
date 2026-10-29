#' CTP distribution
#'
#' Density, distribution, quantile and random-generation functions for the
#' complex triparametric Pearson distribution.
#' @param x Numeric vector of counts. Nonintegers and negative values have
#'   probability zero; missing values propagate.
#' @param a,b,gama Finite scalar parent parameters, with \code{b >= 0} and
#'   \code{gama > max(0, 2*a)}. The spelling \code{gama} is retained.
#' @param log Return log probabilities.
#' @details The probability at zero is evaluated from a shifted log-gamma
#'   expansion. Subsequent log probabilities use the recurrence
#'   \deqn{p(k+1)/p(k) = ((a+k)^2+b^2)/((gama+k)(k+1)).}
#'   Probabilities are not renormalised on the requested support. Moment
#'   existence is checked by moment functions, not by density functions.
#'   Recurrences have a work limit of 1,000,000 counts and signal an error
#'   when that limit is insufficient. At \code{b = 0}, nonpositive integer
#'   \code{a} gives a terminating distribution.
#' @return Numeric vector of probabilities, cumulative probabilities,
#'   quantiles, or integer random draws. Quantiles return \code{Inf} at
#'   probability one for an infinite-support distribution.
#' @examples
#' dctp(0:5, 1, 0.5, 6)
#' pctp(c(-Inf, 0, 2, Inf), 1, 0.5, 6)
#' qctp(c(0, 0.5, 0.99, 1), 1, 0.5, 6)
#' set.seed(123)
#' rctp(10, 1, 0.5, 6)
#' @export
dctp <- function(x, a, b, gama, log = FALSE) {
  .ctp_parameters(a, b, gama); .ctp_flag(log, "log")
  if (!is.numeric(x)) stop("x must be numeric.", call. = FALSE)
  lp <- rep(-Inf, length(x)); lp[is.na(x)] <- NA_real_
  ok <- is.finite(x) & x >= 0 & x == floor(x)
  top <- .ctp_finite(a, b)
  ok <- ok & x <= top
  if (any(ok)) lp[ok] <- .ctp_logpmf(x[ok], a, b, gama)
  if (log) lp else exp(lp)
}

#' @rdname dctp
#' @param q Numeric vector of quantiles; values are rounded down.
#' @param lower.tail Return \eqn{P(X \le q)} if true and \eqn{P(X > q)} otherwise.
#' @param log.p Probabilities are on the log scale.
#' @export
pctp <- function(q, a, b, gama, lower.tail = TRUE, log.p = FALSE) {
  .ctp_cdf(q, a, b, gama, 0, lower.tail, log.p)
}
.ctp_cdf <- function(q, a, b, gama, omega, lower.tail, log.p) {
  .ctp_parameters(a, b, gama); .ctp_omega(omega)
  .ctp_flag(lower.tail, "lower.tail"); .ctp_flag(log.p, "log.p")
  if (!is.numeric(q)) stop("q must be numeric.", call. = FALSE)
  q <- floor(q); ans <- rep(-Inf, length(q)); ans[is.na(q)] <- NA_real_
  top <- if (omega == 1) 0 else .ctp_finite(a, b)
  ans[!is.na(q) & (q == Inf | q >= top)] <- 0
  ok <- is.finite(q) & q >= 0 & q < top
  if (any(ok)) {
    lp <- .ctp_logpmf(0:max(q[ok]), a, b, gama) + log1p(-omega)
    if (omega > 0) lp[1L] <- .ctp_logadd(log(omega), lp[1L])
    cdf <- cumsum(exp(lp))
    if (any(cdf > 1 + 1e-8)) stop("CDF evaluation exceeded one.", call. = FALSE)
    lcdf <- log(pmin(cdf, 1)); lcdf[1L] <- lp[1L]
    ans[ok] <- lcdf[q[ok] + 1L]
  }
  if (!lower.tail) ans <- .ctp_log1mexp(ans)
  if (log.p) ans else exp(ans)
}

#' @rdname dctp
#' @param p Numeric vector of probabilities. Missing values propagate.
#' @param max_q Initial search endpoint, retained for compatibility. It is
#'   expanded as needed, rather than returned as a capped quantile.
#' @param max_count Maximum search endpoint, at most 1,000,000. Failure to
#'   bracket a requested quantile raises an error.
#' @export
qctp <- function(p, a, b, gama, lower.tail = TRUE, log.p = FALSE,
                 max_q = 1000, max_count = 1000000L) {
  .ctp_quantile(p, a, b, gama, 0, lower.tail, log.p, max_q, max_count)
}
.ctp_quantile <- function(p, a, b, gama, omega, lower.tail, log.p,
                          max_q, max_count) {
  .ctp_parameters(a, b, gama); .ctp_omega(omega)
  .ctp_flag(lower.tail, "lower.tail"); .ctp_flag(log.p, "log.p")
  if (!is.numeric(p)) stop("p must be numeric.", call. = FALSE)
  if (any(!is.na(p) & if (log.p) p > 0 else (p < 0 | p > 1)))
    stop("Invalid probability.", call. = FALSE)
  .ctp_scalar(max_q, "max_q"); .ctp_scalar(max_count, "max_count")
  if (max_q < 0 || max_q != floor(max_q) || max_count < max_q ||
      max_count != floor(max_count) || max_count > 1000000)
    stop("Require integer 0 <= max_q <= max_count <= 1000000.", call. = FALSE)
  p <- if (log.p) {if (lower.tail) exp(p) else -expm1(p)} else {
    if (lower.tail) p else 1-p
  }
  top <- if (omega == 1) 0 else .ctp_finite(a, b)
  ans <- rep(NA_real_, length(p)); ans[!is.na(p) & p == 0] <- 0
  ans[!is.na(p) & p == 1] <- top
  ok <- !is.na(p) & p > 0 & p < 1
  if (any(ok)) {
    upper <- min(max_q, top)
    repeat {
      lp <- .ctp_logpmf(0:upper, a, b, gama) + log1p(-omega)
      if (omega > 0) lp[1L] <- .ctp_logadd(log(omega), lp[1L])
      cdf <- cumsum(exp(lp))
      if (any(cdf > 1+1e-8)) stop("CDF evaluation exceeded one.", call. = FALSE)
      cdf <- pmin(cdf, 1)
      if (is.finite(top) && upper == top) cdf[length(cdf)] <- 1
      if (tail(cdf, 1L) >= max(p[ok])) break
      if (upper >= max_count)
        stop("Quantile not bracketed before max_count; no count was capped.", call. = FALSE)
      upper <- min(max_count, top, max(1, 2 * upper))
    }
    ans[ok] <- findInterval(p[ok], cdf, left.open = TRUE)
  }
  ans
}

#' @rdname dctp
#' @param n Number of draws, a nonnegative integer.
#' @export
rctp <- function(n, a, b, gama) {
  n <- .ctp_count_n(n); .ctp_parameters(a, b, gama)
  as.integer(qctp(stats::runif(n), a, b, gama))
}

#' Zero-inflated CTP distribution
#'
#' A mixture of a point mass at zero and a CTP parent. The historical package
#' name uses "zero-modified"; this nonnegative mixing-probability model adds
#' zeros relative to the same parent distribution.
#' @inheritParams dctp
#' @param omega Scalar mixing probability in \code{[0, 1]}. Zero gives the
#'   parent distribution; one gives a point mass at zero. Positive values
#'   increase the zero probability whenever the parent is nondegenerate.
#' @details \eqn{P(Y=0)=omega+(1-omega)P(X=0)} and
#'   \eqn{P(Y=k)=(1-omega)P(X=k)} for \eqn{k>0}.
#'   Quantile search limits and endpoint conventions are as in \code{dctp}.
#' @return Numeric vector of probabilities, cumulative probabilities or
#'   quantiles, or integer random draws.
#' @examples
#' dzictp(0:5, 1, 0.5, 6, 0.3)
#' pzictp(0:5, 1, 0.5, 6, 0.3)
#' qzictp(c(0.1, 0.5, 0.9), 1, 0.5, 6, 0.3)
#' rzictp(10, 1, 0.5, 6, 0.3)
#' @export
dzictp <- function(x, a, b, gama, omega, log = FALSE) {
  .ctp_omega(omega)
  lp <- dctp(x, a, b, gama, log = TRUE) + log1p(-omega)
  zero <- !is.na(x) & x == 0
  if (omega > 0) lp[zero] <- .ctp_logadd(log(omega), lp[zero])
  .ctp_flag(log, "log")
  if (log) lp else exp(lp)
}

#' @rdname dzictp
#' @export
pzictp <- function(q, a, b, gama, omega, lower.tail = TRUE, log.p = FALSE) {
  .ctp_cdf(q, a, b, gama, omega, lower.tail, log.p)
}

#' @rdname dzictp
#' @export
qzictp <- function(p, a, b, gama, omega, lower.tail = TRUE, log.p = FALSE,
                   max_q = 1000, max_count = 1000000L) {
  .ctp_quantile(p, a, b, gama, omega, lower.tail, log.p, max_q, max_count)
}

#' @rdname dzictp
#' @export
rzictp <- function(n, a, b, gama, omega) {
  n <- .ctp_count_n(n); .ctp_parameters(a, b, gama); .ctp_omega(omega)
  as.integer(qzictp(stats::runif(n), a, b, gama, omega))
}
