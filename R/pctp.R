# =============================================================================
# FILE: R/pctp.R
# =============================================================================

#' Cumulative Distribution Function of the CTP Distribution
#'
#' @param q Vector of quantiles.
#' @param a Parameter a.
#' @param b Parameter b.
#' @param gama Parameter gamma.
#' @param lower.tail Logical; if TRUE, probabilities are P(X <= q),
#'   otherwise P(X > q).
#' @param log.p Logical; if TRUE, probabilities are given on the log scale.
#'
#' @return Numeric vector of cumulative probabilities.
#' @export
#'
#' @examples
#' pctp(0:5, a = 1, b = 0.5, gama = 6)
pctp <- function(q, a, b, gama, lower.tail = TRUE, log.p = FALSE) {
  
  q <- floor(q)
  out <- numeric(length(q))
  
  for (i in seq_along(q)) {
    if (is.na(q[i])) {
      out[i] <- NA_real_
    } else if (q[i] < 0) {
      out[i] <- 0
    } else {
      x_vals <- 0:q[i]
      out[i] <- sum(dctp(x_vals, a, b, gama, log = FALSE))
    }
  }
  
  out <- pmin(out, 1)
  
  if (!lower.tail) out <- 1 - out
  if (log.p) out <- log(out)
  
  out
}