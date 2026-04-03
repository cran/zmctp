# =============================================================================
# FILE: R/dzictp.R
# =============================================================================

#' Probability Mass Function of the ZM-CTP Distribution
#'
#' @description
#' The Zero-Modified CTP (ZM-CTP) distribution modifies the zero probability
#' of the baseline CTP distribution.
#'
#' @param x Vector of nonnegative integers.
#' @param a Parameter a.
#' @param b Parameter b.
#' @param gama Parameter gamma.
#' @param omega Zero-modification parameter, with 0 < omega < 1.
#' @param log Logical; if TRUE, returns log-probabilities.
#'
#' @return Numeric vector of probabilities.
#' @export
#'
#' @examples
#' dzictp(0:5, a = 1, b = 0.5, gama = 6, omega = 0.3)
dzictp <- function(x, a, b, gama, omega, log = FALSE) {
  
  if (!is.numeric(omega) || length(omega) != 1 || !is.finite(omega) ||
      omega <= 0 || omega >= 1) {
    stop("omega must be a single finite number in (0, 1).")
  }
  
  x <- as.integer(round(x))
  p_ctp <- dctp(x, a, b, gama, log = FALSE)
  
  prob <- numeric(length(x))
  prob[x == 0] <- omega + (1 - omega) * p_ctp[x == 0]
  prob[x > 0] <- (1 - omega) * p_ctp[x > 0]
  prob[x < 0] <- 0
  
  if (log) return(log(prob))
  prob
}