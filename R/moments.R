#' Moments of the CTP distribution
#'
#' @param a,b,gama Scalar parent parameters; see \code{dctp}.
#' @details For an infinite-support parent the mean, variance and third
#'   moment require \code{gama > 2*a+1}, \code{gama > 2*a+2}, and
#'   \code{gama > 2*a+3}, respectively. A nonexisting mean or variance
#'   returns \code{Inf} with a warning; nonexisting skewness returns
#'   \code{NaN} with a warning. Terminating parents are handled separately.
#' @return A scalar mean, variance or skewness.
#' @examples
#' mean_ctp(1, 0.5, 6)
#' var_ctp(1, 0.5, 6)
#' skew_ctp(1, 0.5, 6)
#' @export
mean_ctp <- function(a, b, gama) .ctp_moments(a,b,gama,1L)
#' @rdname mean_ctp
#' @export
var_ctp <- function(a, b, gama) .ctp_moments(a,b,gama,2L)
#' @rdname mean_ctp
#' @export
skew_ctp <- function(a, b, gama) .ctp_moments(a,b,gama,3L)

.ctp_moments <- function(a,b,gama,order) {
  .ctp_parameters(a,b,gama)
  top <- .ctp_finite(a,b)
  if (is.finite(top)) {
    x <- 0:top; p <- dctp(x,a,b,gama)
    mu <- sum(x*p); v <- sum((x-mu)^2*p)
    return(switch(order,mu,v,if (v>0) sum((x-mu)^3*p)/v^1.5 else NaN))
  }
  if (gama <= 2*a+order) {
    warning("Requested moment does not exist for these parameters.",call.=FALSE)
    return(if (order == 3L) NaN else Inf)
  }
  mu <- (a^2+b^2)/(gama-2*a-1)
  variance <- (a^2+b^2)*((gama-1-a)^2+b^2)/
    ((gama-2*a-1)^2*(gama-2*a-2))
  if (order == 1L) return(mu)
  if (order == 2L) return(variance)
  mu3 <- (a^2+b^2)*(4*b^2+(gama-1)^2)*((gama-1-a)^2+b^2)/
    ((gama-2*a-1)^3*(gama-2*a-2)*(gama-2*a-3))
  mu3/variance^1.5
}

#' Mode of the CTP distribution
#' @param a,b,gama Scalar parent parameters; see \code{dctp}.
#' @return The modal integer, or both adjacent integers in a numerical tie.
#' @details Ties are detected to relative tolerance \code{1e-10}. The
#'   recurrence minus one changes sign at most once in the valid domain.
#' @examples
#' mode_ctp(-1.5, 1.8, 3.8)
#' @export
mode_ctp <- function(a,b,gama) {
  .ctp_parameters(a,b,gama)
  y <- (a^2+b^2-gama)/(gama-2*a+1)
  if (y < 0) return(0L)
  if (abs(y-round(y)) <= 1e-10*max(1,abs(y))) return(c(round(y),round(y)+1))
  ceiling(y)
}

#' Moments of the zero-inflated CTP distribution
#' @inheritParams mean_ctp
#' @param omega Scalar mixing probability in \code{[0, 1]}.
#' @details Moment domains are the same as for the parent, except that
#'   \code{omega = 1} gives mean and variance zero. Skewness of a degenerate
#'   distribution is undefined and returns \code{NaN}.
#' @return A scalar mean, variance or skewness.
#' @examples
#' mean_zictp(1, 0.5, 6, 0.3)
#' var_zictp(1, 0.5, 6, 0.3)
#' skew_zictp(1, 0.5, 6, 0.3)
#' @export
mean_zictp <- function(a,b,gama,omega) {
  .ctp_parameters(a,b,gama); .ctp_omega(omega)
  if (omega == 1) 0 else (1-omega)*mean_ctp(a,b,gama)
}
#' @rdname mean_zictp
#' @export
var_zictp <- function(a,b,gama,omega) {
  .ctp_parameters(a,b,gama); .ctp_omega(omega)
  if (omega == 1) return(0)
  if (omega == 0) return(var_ctp(a,b,gama))
  mu <- mean_ctp(a,b,gama); v <- var_ctp(a,b,gama)
  (1-omega)*v+omega*(1-omega)*mu^2
}
#' @rdname mean_zictp
#' @export
skew_zictp <- function(a,b,gama,omega) {
  .ctp_parameters(a,b,gama); .ctp_omega(omega)
  if (omega == 1) return(NaN)
  if (omega == 0) return(skew_ctp(a,b,gama))
  if (!is.finite(.ctp_finite(a,b)) && gama <= 2*a+3) {
    warning("Requested moment does not exist for these parameters.",call.=FALSE)
    return(NaN)
  }
  mu <- mean_ctp(a,b,gama); v <- var_ctp(a,b,gama)
  if (v == 0) return(NaN)
  mu3 <- skew_ctp(a,b,gama)*v^1.5
  vzi <- var_zictp(a,b,gama,omega)
  ((1-omega)*mu3+3*omega*(1-omega)*mu*v+
     omega*(1-omega)*(2*omega-1)*mu^3)/vzi^1.5
}

#' Mode of the zero-inflated CTP distribution
#' @inheritParams mean_zictp
#' @return The smallest modal integer; ties are resolved to the smallest value.
#' @examples
#' mode_zictp(-1.5, 1.8, 3.8, 0.3)
#' @export
mode_zictp <- function(a,b,gama,omega) {
  .ctp_parameters(a,b,gama); .ctp_omega(omega)
  if (omega == 1) return(0L)
  m <- min(mode_ctp(a,b,gama))
  if (m == 0) return(0L)
  p <- dzictp(c(0,m),a,b,gama,omega,log=TRUE)
  if (p[1] >= p[2]) 0L else m
}
