# =============================================================================
# FILE: R/ctp_fit.R - CORRECTED FITTING FUNCTION
# =============================================================================

#' Maximum Likelihood Estimation for the CTP Distribution
#'
#' @description
#' Fits the Complex Triparametric Pearson (CTP) distribution to count data
#' using maximum likelihood estimation.
#'
#' @param x Numeric vector of nonnegative counts.
#' @param a_start Optional starting value for parameter a.
#' @param b_start Optional starting value for parameter b.
#' @param gama_start Optional starting value for parameter gamma.
#' @param method Optimization method (default: "L-BFGS-B").
#' @param penalty numeric penalty added for numerical stability when b → 0

#'
#' @return An object of class "ctpfit" containing:
#' \item{estimates}{Named vector of MLEs}
#' \item{se}{Standard errors}
#' \item{vcov}{Variance-covariance matrix}
#' \item{logLik}{Log-likelihood}
#' \item{AIC}{Akaike Information Criterion}
#' \item{BIC}{Bayesian Information Criterion}
#' \item{pearson_chisq}{Pearson's chi-squared statistic}
#' \item{wald_chisq}{Wald's chi-squared statistic}
#' \item{fitted_freq}{Data frame of observed vs expected frequencies}
#' \item{data}{Original data}
#' \item{converged}{Convergence status}
#'
#' @export
#' @examples
#' set.seed(123)
#' \donttest{
#' x <- rctp(30, a = 1, b = 0.5, gama = 5)
#' fit <- ctp.fit(x)
#' print(fit)
#' plot(fit)
#' }
ctp.fit <- function(x, a_start = NULL, b_start = NULL, gama_start = NULL,
                    method = "L-BFGS-B", penalty = 1e10) {
  
  x <- as.numeric(x)
  N <- length(x)
  mu_emp <- mean(x)
  var_emp <- var(x)
  
  # sensible fallback heuristics for starts
  theta_nb <- ifelse(var_emp > mu_emp, mu_emp^2 / (var_emp - mu_emp), 1)
  theta_nb <- max(theta_nb, 0.5)
  
  if (is.null(gama_start)) gama_start <- max(4.0, theta_nb * 1.5 + 2)
  if (is.null(a_start)) a_start <- max(1e-3, mu_emp * (gama_start / (gama_start + mu_emp)))
  if (is.null(b_start)) b_start <- 0.05
  
  # reparameterize gama = 2*a + 2 + exp(eta) so variance-existence constraint holds
  eta_start <- log(max(1e-6, gama_start - 2 * a_start - 2))
  
  # parameters we optimize: a, b, eta
  start_vals <- c(a = a_start, b = b_start, eta = eta_start)
  
  # negative log-likelihood in new parametrization
  nll <- function(theta, x) {
    a <- theta[1]
    b <- theta[2]
    eta <- theta[3]
    # re-build gama
    gama <- 2 * a + 2 + exp(eta)
    
    # basic domain checks (return penalty if invalid)
    # Note: a can be any real number (negative a gives under-dispersion)
    if (!is.finite(a) || !is.finite(b) || !is.finite(gama)) return(penalty)
    if (b <= 0 || exp(eta) <= 0) return(penalty)
    
    # call user-provided density; dctp must accept log=TRUE
    log_probs <- tryCatch({
      dctp(x, a, b, gama, log = TRUE)
    }, error = function(e) {
      rep(-Inf, length(x))
    })
    
    if (any(!is.finite(log_probs))) return(penalty)
    val <- -sum(log_probs)
    # guard: if val is absurdly large (penalty region), return penalty
    if (!is.finite(val) || val > penalty) return(penalty)
    return(val)
  }
  
  # bounds: a unrestricted (negative a gives under-dispersion), b>0, eta unrestricted
  lower_bounds <- c(-Inf, 1e-8, -30)   # eta lower ~ exp(-30) ~ 9e-14
  upper_bounds <- c(Inf, Inf,  30)
  
  fit <- optim(
    par = start_vals,
    fn = nll,
    x = x,
    method = method,
    lower = lower_bounds,
    upper = upper_bounds,
    hessian = TRUE,
    control = list(maxit = 10000)
  )
  
  # convert back to natural params
  est_raw <- fit$par
  est <- numeric(3)
  names(est) <- c("a", "b", "gama")
  est["a"] <- est_raw["a"]
  est["b"] <- est_raw["b"]
  est["gama"] <- 2 * est["a"] + 2 + exp(est_raw["eta"])
  
  # logLik handling: if optimizer returned penalty, mark as NA
  if (!is.null(fit$value) && is.finite(fit$value) && fit$value < penalty) {
    logLik <- -fit$value
  } else {
    logLik <- NA_real_
  }
  
  # Hessian -> vcov with jitter to improve inversion robustness
  vcov_mat <- NULL
  se <- NULL
  if (!is.null(fit$hessian) && all(is.finite(fit$hessian))) {
    # need to transform covariance from (a,b,eta) -> (a,b,gama)
    # get vcov in (a,b,eta) space
    H <- fit$hessian
    # try inversion with small jitter
    inv_ok <- FALSE
    jitter <- 0
    while (!inv_ok && jitter <= 1e-4) {
      try({
        V_raw <- try(solve(H + diag(jitter, nrow(H))), silent = TRUE)
        if (is.matrix(V_raw) && all(is.finite(V_raw))) {
          inv_ok <- TRUE
        } else {
          inv_ok <- FALSE
        }
      }, silent = TRUE)
      if (!inv_ok) jitter <- ifelse(jitter == 0, 1e-8, jitter * 10)
    }
    if (inv_ok) {
      V_raw <- solve(H + diag(jitter, nrow(H)))
      # jacobian of transform (a,b,eta) -> (a,b,gama)
      # gama = 2*a + 2 + exp(eta) -> d(gama)/d(a) = 2, d(gama)/d(b)=0, d(gama)/d(eta)=exp(eta)
      a_hat <- est_raw["a"]; b_hat <- est_raw["b"]; eta_hat <- est_raw["eta"]
      J <- matrix(0, nrow = 3, ncol = 3)
      colnames(J) <- rownames(J) <- c("a","b","gama")
      # map order: (a,b,eta) -> (a,b,gama)
      ## ∂a/∂a = 1, ∂a/∂b=0, ∂a/∂eta=0
      J[1,1] <- 1; J[1,2] <- 0; J[1,3] <- 0
      ## ∂b/∂a = 0, ∂b/∂b=1, ∂b/∂eta=0
      J[2,1] <- 0; J[2,2] <- 1; J[2,3] <- 0
      ## ∂gama/∂a = 2, ∂gama/∂b=0, ∂gama/∂eta=exp(eta)
      J[3,1] <- 2; J[3,2] <- 0; J[3,3] <- exp(est_raw["eta"])
      # Note: V_raw is cov for (a,b,eta); transform to (a,b,gama) via J %*% V_raw %*% t(J)
      V_trans <- J %*% V_raw %*% t(J)
      colnames(V_trans) <- rownames(V_trans) <- c("a","b","gama")
      vcov_mat <- V_trans
      # standard errors:
      se <- tryCatch({
        s <- sqrt(pmax(0, diag(vcov_mat)))
        names(s) <- c("a","b","gama")
        s
      }, error = function(e) NULL)
    } else {
      warning("Hessian inversion failed (even after jitter); SE not computed.")
    }
  } else {
    warning("No finite Hessian available; SE not computed.")
  }
  
  # Goodness-of-fit: observed vs expected frequencies (safe)
  obs_freq <- as.data.frame(table(x))
  obs_freq$x <- as.numeric(as.character(obs_freq$x))
  obs_freq$expected <- tryCatch({
    N * dctp(obs_freq$x, est["a"], est["b"], est["gama"])
  }, error = function(e) {
    rep(NA_real_, nrow(obs_freq))
  })
  
  # remove rows with NA expected (e.g. if density fails)
  obs_freq_valid <- obs_freq[is.finite(obs_freq$expected) & obs_freq$expected > 0, , drop = FALSE]
  
  pearson_chisq <- if (nrow(obs_freq_valid) > 0) {
    sum((obs_freq_valid$Freq - obs_freq_valid$expected)^2 / obs_freq_valid$expected)
  } else NA_real_
  
  # theoretical moments: compute only if constraints satisfied
  mu_fitted <- NA_real_; var_fitted <- NA_real_; wald_chisq <- NA_real_
  try({
    if (est["gama"] > 2 * est["a"] + 1) {
      mu_fitted <- (est["a"]^2 + est["b"]^2) / (est["gama"] - 2 * est["a"] - 1)
    }
    if (est["gama"] > 2 * est["a"] + 2) {
      var_fitted <- mu_fitted * (mu_fitted + est["gama"] - 1) / (est["gama"] - 2 * est["a"] - 2)
    }
    if (is.finite(var_fitted) && var_fitted > 0) {
      wald_chisq <- sum((x - mu_fitted)^2 / var_fitted)
    }
  }, silent = TRUE)
  
  # information criteria (safely)
  k <- length(est)
  AIC <- if (!is.na(logLik)) -2 * logLik + 2 * k else NA_real_
  BIC <- if (!is.na(logLik)) -2 * logLik + k * log(N) else NA_real_
  
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
    converged = (is.numeric(fit$convergence) && fit$convergence == 0),
    method = method,
    raw_optim = fit
  )
  class(result) <- "ctpfit"
  return(result)
}