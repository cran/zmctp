# =============================================================================
# FILE: R/moments.R  
# =============================================================================

#' Theoretical moments and mode of the CTP distribution
#'
#' @param a,b,gama CTP parameters.
#'
#' @details
#' \code{mode_ctp()} returns the unimodal-case mode (rounded). The
#' distribution has two consecutive modes in the knife-edge case where
#' \code{((a-1)^2+b^2)/(gama-2*a+1)} is exactly an integer; this edge case
#' is not separately flagged.
#'
#' These functions do not validate that the parameter constraints required
#' for a given moment to exist are satisfied (e.g. \code{gama > 2*a+2} for
#' the variance, \code{gama > 2*a+3} for the skewness). Outside these
#' regions the returned value may be \code{Inf}, \code{NaN}, or negative;
#' see \code{\link{dctp}} for the corresponding existence warnings applied
#' at the density level.
#'
#' @examples
#' mean_ctp(a = 1, b = 0.5, gama = 6)
#' var_ctp(a = 1, b = 0.5, gama = 6)
#' skew_ctp(a = 1, b = 0.5, gama = 6)
#' mode_ctp(a = 1, b = 0.5, gama = 6)
#'
#' @export
mean_ctp <- function(a, b, gama) (a^2 + b^2) / (gama - 2*a - 1)

#' @rdname mean_ctp
#' @export
var_ctp <- function(a, b, gama) {
  mu <- mean_ctp(a, b, gama)
  mu * (mu + gama - 1) / (gama - 2*a - 2)
}

#' @rdname mean_ctp
#' @export
skew_ctp <- function(a, b, gama) {
  mu3 <- ((a^2+b^2)*(4*b^2+(gama-1)^2) + (b^2+(gama-1-a)^2)) /
    ((gama-2*a-1)^3 * (gama-2*a-2) * (gama-2*a-3))
  mu3 / var_ctp(a, b, gama)^1.5
}

#' @rdname mean_ctp
#' @export
mode_ctp <- function(a, b, gama) {
  round(((a-1)^2 + b^2) / (gama - 2*a + 1))
}

#' Theoretical moments and mode of the ZM-CTP distribution
#'
#' @param a,b,gama CTP parameters.
#' @param omega Zero-modification parameter.
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

#' @rdname mean_zictp
#' @export
mode_zictp <- function(a, b, gama, omega) {
  f0 <- dctp(0, a, b, gama)
  ystar <- mode_ctp(a, b, gama)
  fystar <- dctp(ystar, a, b, gama)
  if (omega / (1 - omega) >= fystar - f0) 0 else ystar
}