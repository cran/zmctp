# =============================================================================
# FILE: R/zictp_fit.R
# =============================================================================

#' Maximum Likelihood Estimation for the Zero-Modified CTP Distribution
#'
#' @description
#' Fits the Zero-Modified Complex TriParametric Pearson (ZM-CTP) distribution
#' to count data using maximum likelihood estimation.
#'
#' A logit reparameterization is used for omega to ensure 0 < omega < 1,
#' and a log reparameterization is used for gamma so that
#' gama = 2*a + 2 + exp(eta), which guarantees variance existence
#' throughout optimization.
#'
#' @param x Numeric vector of nonnegative counts.
#' @param a_start Optional starting value for parameter a.
#' @param b_start Optional starting value for parameter b.
#' @param gama_start Optional starting value for parameter gamma.
#' @param omega_start Optional starting value for omega (0 < omega < 1).
#' @param method Optimization method (default: "BFGS").
#'
#' @return An object of class `"zictpfit"`.
#' @export
#' @examples
#' \donttest{
#' set.seed(123)
#' x <- rzictp(30, a = 1, b = 0.5, gama = 5, omega = 0.3)
#' fit <- zictp.fit(x)
#' fit$estimates
#' }
zictp.fit <- function(x, a_start = NULL, b_start = NULL, gama_start = NULL,
                      omega_start = NULL, method = "BFGS") {
  
  x <- as.numeric(x)
  
  if (length(x) == 0) stop("x must not be empty.")
  if (any(!is.finite(x))) stop("x must contain only finite values.")
  if (any(x < 0)) stop("x must contain nonnegative counts.")
  if (any(abs(x - round(x)) > sqrt(.Machine$double.eps))) {
    stop("x must contain integer-valued counts.")
  }
  
  x <- as.integer(round(x))
  N <- length(x)
  
  mu_emp <- mean(x)
  var_emp <- stats::var(x)
  if (is.na(var_emp)) var_emp <- mu_emp + 1
  
  p_zero_obs <- mean(x == 0)
  
  theta_nb <- ifelse(var_emp > mu_emp, mu_emp^2 / (var_emp - mu_emp), 1)
  theta_nb <- max(theta_nb, 0.5)
  
  if (is.null(gama_start)) gama_start <- max(4.5, theta_nb * 1.5 + 2)
  if (is.null(a_start)) a_start <- max(0.05, mu_emp * (gama_start / (gama_start + mu_emp + 1e-8)))
  if (is.null(b_start)) b_start <- max(0.05, stats::sd(x) / 4)
  if (is.null(omega_start)) omega_start <- min(0.8, max(0.02, p_zero_obs - 0.05))
  
  omega_start <- min(0.95, max(0.001, omega_start))
  
  nll <- function(theta, x) {
    a <- theta[1]
    b <- theta[2]
    eta <- theta[3]
    delta <- theta[4]
    
    gama <- 2 * a + 2 + exp(eta)
    omega <- plogis(delta)
    
    if (!is.finite(a) || !is.finite(b) || !is.finite(eta) || !is.finite(delta)) {
      return(1e10)
    }
    # Note: a can be any real; only b must be >= 0
    if (b < 0) {
      return(1e10)
    }
    
    log_p_ctp <- tryCatch(
      suppressWarnings(dctp(x, a, b, gama, log = TRUE)),
      error = function(e) rep(-Inf, length(x))
    )
    
    if (any(!is.finite(log_p_ctp))) return(1e10)
    
    p_ctp_0 <- tryCatch(
      suppressWarnings(dctp(0, a, b, gama, log = FALSE)),
      error = function(e) NA_real_
    )
    
    if (!is.finite(p_ctp_0) || p_ctp_0 <= 0 || p_ctp_0 >= 1) return(1e10)
    
    n_zero <- sum(x == 0)
    n_pos <- sum(x > 0)
    
    p_zero_zmctp <- omega + (1 - omega) * p_ctp_0
    if (!is.finite(p_zero_zmctp) || p_zero_zmctp <= 0 || p_zero_zmctp >= 1) {
      return(1e10)
    }
    
    ll_zero <- n_zero * log(p_zero_zmctp)
    ll_positive <- if (n_pos > 0) {
      sum(log1p(-omega) + log_p_ctp[x > 0])
    } else {
      0
    }
    
    -(ll_zero + ll_positive)
  }
  
  # base start
  eta_start <- log(max(1e-6, gama_start - 2 * a_start - 2))
  delta_start <- qlogis(omega_start)
  
  start_grid <- list(
    c(a = a_start,          b = b_start,          eta = eta_start,                                      delta = delta_start),
    c(a = max(0.05, 0.5*a_start), b = max(0.05, 0.5*b_start), eta = log(max(1e-6, gama_start - 2*(max(0.05, 0.5*a_start)) - 2)), delta = qlogis(min(0.95, max(0.01, 0.5*omega_start)))),
    c(a = max(0.05, 1.5*a_start), b = max(0.05, 1.5*b_start), eta = log(max(1e-6, (gama_start + 1) - 2*(max(0.05, 1.5*a_start)) - 2)), delta = qlogis(min(0.95, max(0.01, 1.5*omega_start)))),
    c(a = max(0.05, mu_emp/4),    b = max(0.05, stats::sd(x)/5), eta = log(max(1e-6, 6 - 2*(max(0.05, mu_emp/4)) - 2)),         delta = qlogis(0.05)),
    c(a = max(0.05, mu_emp/2),    b = max(0.05, stats::sd(x)/3), eta = log(max(1e-6, 8 - 2*(max(0.05, mu_emp/2)) - 2)),         delta = qlogis(0.20))
  )
  
  fits <- lapply(start_grid, function(st) {
    try(
      stats::optim(
        par = st,
        fn = nll,
        x = x,
        method = method,
        hessian = TRUE,
        control = list(maxit = 10000, reltol = 1e-10)
      ),
      silent = TRUE
    )
  })
  
  ok <- Filter(function(f) !inherits(f, "try-error") && is.finite(f$value), fits)
  if (length(ok) == 0) stop("Optimization failed for all starting values.")
  
  best_idx <- which.min(vapply(ok, function(f) f$value, numeric(1)))
  fit <- ok[[best_idx]]
  
  a_hat <- unname(fit$par["a"])
  b_hat <- unname(fit$par["b"])
  eta_hat <- unname(fit$par["eta"])
  delta_hat <- unname(fit$par["delta"])
  
  gama_hat <- 2 * a_hat + 2 + exp(eta_hat)
  omega_hat <- plogis(delta_hat)
  
  est <- c(a = a_hat, b = b_hat, gama = gama_hat, omega = omega_hat)
  
  vcov_mat <- NULL
  se <- NULL
  
  if (!is.null(fit$hessian) &&
      all(is.finite(fit$hessian)) &&
      nrow(fit$hessian) == 4 &&
      ncol(fit$hessian) == 4) {
    
    H <- fit$hessian
    inv_ok <- FALSE
    jitter <- 0
    
    while (!inv_ok && jitter <= 1e-4) {
      V_raw <- try(solve(H + diag(jitter, 4)), silent = TRUE)
      if (is.matrix(V_raw) && all(is.finite(V_raw))) {
        inv_ok <- TRUE
      } else {
        jitter <- ifelse(jitter == 0, 1e-8, jitter * 10)
      }
    }
    
    if (inv_ok) {
      V_raw <- solve(H + diag(jitter, 4))
      
      J <- diag(4)
      J[3, 1] <- 2
      J[3, 3] <- exp(eta_hat)
      J[4, 4] <- omega_hat * (1 - omega_hat)
      
      V_trans <- J %*% V_raw %*% t(J)
      colnames(V_trans) <- rownames(V_trans) <- c("a", "b", "gama", "omega")
      
      vcov_mat <- V_trans
      se <- sqrt(pmax(diag(V_trans), 0))
      names(se) <- c("a", "b", "gama", "omega")
    } else {
      warning("Hessian inversion failed; SE not computed.")
    }
  }
  
  obs_freq <- as.data.frame(table(x), stringsAsFactors = FALSE)
  obs_freq$x <- as.numeric(obs_freq$x)
  obs_freq$expected <- N * dzictp(obs_freq$x, a_hat, b_hat, gama_hat, omega_hat)
  
  valid_exp <- is.finite(obs_freq$expected) & obs_freq$expected > 0
  pearson_chisq <- if (any(valid_exp)) {
    sum((obs_freq$Freq[valid_exp] - obs_freq$expected[valid_exp])^2 / obs_freq$expected[valid_exp])
  } else {
    NA_real_
  }
  
  mu_fitted <- NA_real_
  var_fitted <- NA_real_
  
  if (gama_hat > 2 * a_hat + 1) {
    mu_ctp <- (a_hat^2 + b_hat^2) / (gama_hat - 2 * a_hat - 1)
    mu_fitted <- (1 - omega_hat) * mu_ctp
  }
  
  if (gama_hat > 2 * a_hat + 2 && is.finite(mu_fitted)) {
    var_ctp <- mu_ctp * (mu_ctp + gama_hat - 1) / (gama_hat - 2 * a_hat - 2)
    var_fitted <- (1 - omega_hat) * (var_ctp + mu_ctp^2) - mu_fitted^2
  }
  
  wald_chisq <- if (is.finite(var_fitted) && var_fitted > 0) {
    sum((x - mu_fitted)^2 / var_fitted)
  } else {
    NA_real_
  }
  
  logLik <- -fit$value
  k <- 4
  AIC <- -2 * logLik + 2 * k
  BIC <- -2 * logLik + k * log(N)
  
  result <- list(
    estimates = est,
    se = se,
    vcov = vcov_mat,
    logLik = logLik,
    AIC = AIC,
    BIC = BIC,
    pearson_chisq = pearson_chisq,
    wald_chisq = wald_chisq,
    fitted_freq = obs_freq,
    data = x,
    converged = (fit$convergence == 0),
    method = method,
    raw_optim = fit
  )
  
  class(result) <- "zictpfit"
  return(result)
}