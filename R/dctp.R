# =============================================================================
# FILE: R/dctp.R
# =============================================================================

#' Probability Mass Function of the CTP Distribution
#'
#' @description
#' Evaluates the probability mass function of the Complex TriParametric
#' Pearson (CTP) distribution for nonnegative integer values.
#'
#' The pmf is
#' \deqn{p(x) = p(0)\frac{(a+ib)_x(a-ib)_x}{(\gamma)_x x!}, \quad x=0,1,\dots}
#' and satisfies the recurrence
#' \deqn{\frac{p(x+1)}{p(x)} = \frac{(a+x)^2+b^2}{(\gamma+x)(x+1)}.}
#'
#' @param x Vector of nonnegative integers.
#' @param a Parameter a.
#' @param b Parameter b (must satisfy b >= 0).
#' @param gama Parameter gamma.
#' @param log Logical; if TRUE, returns log-probabilities.
#'
#' @return A numeric vector of probabilities.
#' @export
#'
#' @examples
#' dctp(0:5, a = 1, b = 0.5, gama = 6)
dctp <- function(x, a, b, gama, log = FALSE) {
  
  if (!is.numeric(x)) stop("x must be numeric.")
  if (any(!is.finite(x))) stop("x must contain only finite values.")
  if (any(x < 0) || any(abs(x - round(x)) > sqrt(.Machine$double.eps))) {
    stop("x must contain nonnegative integers.")
  }
  if (!is.numeric(a) || length(a) != 1 || !is.finite(a)) {
    stop("a must be a single finite number.")
  }
  if (!is.numeric(b) || length(b) != 1 || !is.finite(b) || b < 0) {
    stop("b must be a single finite number with b >= 0.")
  }
  if (!is.numeric(gama) || length(gama) != 1 || !is.finite(gama) || gama <= 0) {
    stop("gama must be a single positive finite number.")
  }
  
  x <- as.integer(round(x))
  
  # Moment conditions from the paper:
  # mean exists if gama > 2*a + 1
  # variance exists if gama > 2*a + 2
  # third central moment exists if gama > 2*a + 3
  if (gama <= 2 * a + 2) {
    warning("Variance does not exist because gama <= 2*a + 2.")
  }
  if (gama <= 2 * a + 3) {
    warning("Third central moment does not exist because gama <= 2*a + 3.")
  }
  
  max_x <- max(x)
  
  # Build unnormalized weights via recurrence:
  # w(0)=1
  # w(x+1)=w(x)*((a+x)^2+b^2)/((gama+x)(x+1))
  w <- numeric(max_x + 1L)
  w[1L] <- 1
  
  if (max_x >= 1L) {
    for (k in 0:(max_x - 1L)) {
      ratio <- ((a + k)^2 + b^2) / ((gama + k) * (k + 1))
      if (!is.finite(ratio) || ratio < 0) {
        stop("Invalid recurrence ratio encountered.")
      }
      w[k + 2L] <- w[k + 1L] * ratio
    }
  }
  
  # Extend tail until negligible for normalization.
  # Preallocate the extension buffer instead of growing it one element at a
  # time via c(): near the gama = 2*a + 2 boundary the tail decays very
  # slowly and can need close to the full 10000-term cap, and repeated
  # vector growth there is O(n^2) -- easily seconds per call, which then
  # multiplies out to minutes/hours across an optimizer's function
  # evaluations.
  extra_cap <- 10000L
  tail_extra <- numeric(extra_cap)
  last_w <- w[max_x + 1L]
  k <- max_x
  n_extra <- 0L
  repeat {
    ratio <- ((a + k)^2 + b^2) / ((gama + k) * (k + 1))
    next_w <- last_w * ratio

    if (!is.finite(next_w) || next_w < 0) {
      stop("Failed while computing normalization tail.")
    }

    n_extra <- n_extra + 1L
    tail_extra[n_extra] <- next_w
    last_w <- next_w
    k <- k + 1L

    # stop when the added term is tiny
    if (next_w < 1e-12 || n_extra >= extra_cap) break
  }

  denom <- sum(w) + sum(tail_extra[seq_len(n_extra)])
  if (!is.finite(denom) || denom <= 0) {
    stop("Failed to compute a valid normalization constant.")
  }
  
  prob <- w[x + 1L] / denom
  
  if (log) return(log(prob))
  prob
}