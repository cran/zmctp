# =============================================================================
# FILE: R/moments.R  
# =============================================================================

#' Theoretical moments of the CTP distribution
#'
#' @param a,b,gama CTP parameters.
#'
#' @details
#' These functions do not validate that the parameter constraints required
#' for a given moment to exist are satisfied (e.g. \code{gama > 2*a+2} for
#' the variance, \code{gama > 2*a+3} for the skewness). Outside these
#' regions the returned value may be \code{Inf}, \code{NaN}, or negative;
#' see \code{\link{dctp}} for the corresponding existence warnings applied
#' at the density level. See \code{\link{mode_ctp}} for the mode, including
#' its tie-case handling.
#'
#' @examples
#' mean_ctp(a = 1, b = 0.5, gama = 6)
#' var_ctp(a = 1, b = 0.5, gama = 6)
#' skew_ctp(a = 1, b = 0.5, gama = 6)
#'
#' @export
mean_ctp <- function(a, b, gama) (a^2 + b^2) / (gama - 2*a - 1)

#' @rdname mean_ctp
#' @export
var_ctp <- function(a, b, gama) {
  (a^2 + b^2) * ((gama - 1 - a)^2 + b^2) / ((gama - 2*a - 1)^2 * (gama - 2*a - 2))
}

#' @rdname mean_ctp
#' @export
skew_ctp <- function(a, b, gama) {
  mu3 <- (a^2+b^2) * (4*b^2+(gama-1)^2) * ((gama-1-a)^2+b^2) /
    ((gama-2*a-1)^3 * (gama-2*a-2) * (gama-2*a-3))
  mu3 / var_ctp(a, b, gama)^1.5
}

#' Mode of the CTP distribution
#'
#' @param a alpha parameter
#' @param b b parameter (b >= 0)
#' @param gama gamma parameter
#' @return the mode (a single integer), or, in the exact tie case,
#'   a numeric vector of length 2 giving both consecutive modes
#' @export
mode_ctp <- function(a, b, gama) {
  y0 <- (a^2 + b^2 - gama) / (gama - 2*a + 1)
  
  if (y0 < 0) {
    return(0L)
  }
  
  # Check for the exact tie case (y0 is an integer) using a small
  # numerical tolerance, since y0 is a floating-point calculation
  tol <- 1e-8
  if (abs(y0 - round(y0)) < tol) {
    y0_int <- round(y0)
    return(c(y0_int, y0_int + 1))  # two consecutive modes
  }
  
  # Generic case: mode is the ceiling of y0
  return(as.integer(ceiling(y0)))
}

#' Theoretical moments and mode of the ZI-CTP distribution
#'
#' @param a,b,gama CTP parameters.
#' @param omega Zero-inflation parameter.
#'
#' @details
#' These functions do not validate that the parameter constraints required
#' for a given moment to exist are satisfied (e.g. \code{gama > 2*a+2} for
#' the variance, \code{gama > 2*a+3} for the skewness), the same caveat
#' documented in \code{\link{mean_ctp}}.
#'
#' @examples
#' # Quick illustration (fast, runs during R CMD check)
#' mean_zictp(a = 1, b = 0.5, gama = 6, omega = 0.3)
#' var_zictp(a = 1, b = 0.5, gama = 6, omega = 0.3)
#' skew_ctp(a = 1, b = 0.5, gama = 6)
#' skew_zictp(a = 1, b = 0.5, gama = 6, omega = 0.3)
#' mode_zictp(a = 1, b = 0.5, gama = 6, omega = 0.3)
#'
#' \donttest{
#' # Empirical check against simulated data (slower, skipped on CRAN checks)
#' set.seed(2026)
#' x <- rzictp(50000, a = 1, b = 0.5, gama = 6, omega = 0.3)
#' skew_zictp(1, 0.5, 6, 0.3)                       # theoretical
#' mean((x - mean(x))^3) / sd(x)^3                  # empirical
#' }
#'
#' @export
mean_zictp <- function(a, b, gama, omega) (1 - omega) * mean_ctp(a, b, gama)

#' @rdname mean_zictp
#' @export
var_zictp <- function(a, b, gama, omega) {
  mu <- mean_ctp(a, b, gama); s2 <- var_ctp(a, b, gama)
  (1 - omega) * s2 + omega * (1 - omega) * mu^2
}

#' @rdname mean_zictp
#' @export
skew_zictp <- function(a, b, gama, omega) {
  mu  <- mean_ctp(a, b, gama)
  s2  <- var_ctp(a, b, gama)
  mu3 <- skew_ctp(a, b, gama) * s2^1.5   # recover CTP's raw mu3
  mu3_zm <- (1-omega)*mu3 + 3*omega*(1-omega)*mu*s2 +
    omega*(1-omega)*(2*omega-1)*mu^3
  mu3_zm / var_zictp(a, b, gama, omega)^1.5
}

#' Mode of the ZI-CTP distribution
#'
#' @param a alpha parameter
#' @param b b parameter (b >= 0)
#' @param gama gamma parameter
#' @param omega zero-inflation parameter
#' @return the mode (a single integer)
#' @export
mode_zictp <- function(a, b, gama, omega) {
  y_star <- mode_ctp(a, b, gama)
  
  # Handle the rare tie case from mode_ctp() by using the smaller
  # of the two consecutive modes for this comparison (a reasonable,
  # documented convention -- see package documentation)
  if (length(y_star) > 1) {
    y_star <- y_star[1]
  }
  
  if (y_star == 0) {
    return(0L)  # mode is trivially 0 already
  }
  
  f0 <- dctp(0, a, b, gama)
  f_ystar <- dctp(y_star, a, b, gama)
  
  threshold <- omega / (1 - omega)
  
  if (threshold >= (f_ystar - f0)) {
    return(0L)
  } else {
    return(as.integer(y_star))
  }
}