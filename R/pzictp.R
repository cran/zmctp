# =============================================================================
# FILE: R/pzictp.R
# =============================================================================

#' Cumulative Distribution Function of the ZM-CTP Distribution
#'
#' @param q Vector of quantiles.
#' @param a Parameter a.
#' @param b Parameter b.
#' @param gama Parameter gamma.
#' @param omega Zero-modification parameter.
#' @param lower.tail Logical.
#' @param log.p Logical.
#'
#' @return Numeric vector of cumulative probabilities.
#' @export
#'
#' @examples
#' pzictp(0:5, a = 1, b = 0.5, gama = 6, omega = 0.3)
pzictp <- function(q, a, b, gama, omega, lower.tail = TRUE, log.p = FALSE) {
  
  q <- floor(q)
  out <- numeric(length(q))
  
  for (i in seq_along(q)) {
    if (is.na(q[i])) {
      out[i] <- NA_real_
    } else if (q[i] < 0) {
      out[i] <- 0
    } else {
      x_vals <- 0:q[i]
      out[i] <- sum(dzictp(x_vals, a, b, gama, omega, log = FALSE))
    }
  }
  
  out <- pmin(out, 1)
  
  if (!lower.tail) out <- 1 - out
  if (log.p) out <- log(out)
  
  out
}