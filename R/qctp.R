# =============================================================================
# FILE: R/qctp.R
# =============================================================================

#' Quantile Function of the CTP Distribution
#'
#' @param p Vector of probabilities.
#' @param a Parameter a.
#' @param b Parameter b.
#' @param gama Parameter gamma.
#' @param lower.tail Logical.
#' @param log.p Logical.
#' @param max_q Maximum x to search.
#'
#' @return Numeric vector of quantiles.
#' @export
#'
#' @examples
#' qctp(c(0.25, 0.5, 0.75), a = 1, b = 0.5, gama = 6)
qctp <- function(p, a, b, gama, lower.tail = TRUE, log.p = FALSE, max_q = 1000) {
  
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  
  if (any(is.na(p))) stop("p contains NA values.")
  if (any(p < 0 | p > 1)) stop("p must be in [0, 1].")
  
  x_vals <- 0:max_q
  probs <- dctp(x_vals, a, b, gama, log = FALSE)
  cdf <- cumsum(probs)
  cdf[cdf > 1] <- 1
  
  out <- vapply(p, function(pi) {
    idx <- which(cdf >= pi)[1]
    if (is.na(idx)) max_q else x_vals[idx]
  }, numeric(1))
  
  out
}