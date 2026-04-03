# =============================================================================
# FILE: R/rctp.R
# =============================================================================

#' Random Generation from the CTP Distribution
#'
#' @param n Number of observations.
#' @param a Parameter a.
#' @param b Parameter b.
#' @param gama Parameter gamma.
#'
#' @return Integer vector of random draws.
#' @export
#'
#' @examples
#' set.seed(123)
#' rctp(10, a = 1, b = 0.5, gama = 6)
rctp <- function(n, a, b, gama) {
  
  if (length(n) != 1 || !is.numeric(n) || is.na(n) || n < 0) {
    stop("n must be a single nonnegative number.")
  }
  
  n <- as.integer(n)
  if (n == 0) return(integer(0))
  
  u <- stats::runif(n)
  as.integer(qctp(u, a, b, gama))
}