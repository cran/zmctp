# =============================================================================
# FILE: R/rzictp.R
# =============================================================================

#' Random Generation from the ZM-CTP Distribution
#'
#' @param n Number of observations.
#' @param a Parameter a.
#' @param b Parameter b.
#' @param gama Parameter gamma.
#' @param omega Zero-modification parameter.
#'
#' @return Integer vector of random draws.
#' @export
#'
#' @examples
#' set.seed(123)
#' rzictp(10, a = 1, b = 0.5, gama = 6, omega = 0.3)
rzictp <- function(n, a, b, gama, omega) {
  
  if (length(n) != 1 || !is.numeric(n) || is.na(n) || n < 0) {
    stop("n must be a single nonnegative number.")
  }
  
  n <- as.integer(n)
  if (n == 0) return(integer(0))
  
  u <- stats::runif(n)
  as.integer(qzictp(u, a, b, gama, omega))
}