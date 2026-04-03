# =============================================================================
# FILE: R/methods.R - S3 METHODS
# =============================================================================

#' Print method for ctpfit objects
#' @param x A ctpfit object
#' @param ... Additional arguments (ignored)
#' @return
#' Invisibly returns the original fitted model object.
#' The function is called for its side effect of printing a model summary.
#' @export
print.ctpfit <- function(x, ...) {
  cat("\nCTP Distribution - Maximum Likelihood Estimates\n")
  cat("===============================================\n\n")

  cat("Parameter Estimates:\n")
  est_table <- data.frame(
    Estimate = round(x$estimates, 4),
    Std.Error = if (!is.null(x$se)) round(x$se, 4) else rep(NA, length(x$estimates))
  )
  print(est_table)

  cat("\nGoodness-of-Fit Statistics:\n")
  cat(sprintf("  Log-Likelihood: %.4f\n", x$logLik))
  cat(sprintf("  AIC: %.4f\n", x$AIC))
  cat(sprintf("  BIC: %.4f\n", x$BIC))
  cat(sprintf("  Pearson Chi-sq: %.4f\n", x$pearson_chisq))
  cat(sprintf("  Wald Chi-sq: %.4f\n", x$wald_chisq))

  cat("\nConvergence:", ifelse(x$converged, "YES", "NO"), "\n")
  invisible(x)
}
#' Print method for zictpfit objects
#' @param x A zictpfit object
#' @param ... Additional arguments (ignored)
#' @return
#' Invisibly returns the original fitted model object.
#' The function is called for its side effect of printing a model summary.
#' @export
print.zictpfit <- function(x, ...) {
  cat("\nZero-Modified CTP Distribution - Maximum Likelihood Estimates\n")
  cat("=============================================================\n\n")

  cat("Parameter Estimates:\n")
  est_table <- data.frame(
    Estimate = round(x$estimates, 4),
    Std.Error = if (!is.null(x$se)) round(x$se, 4) else rep(NA, length(x$estimates))
  )
  print(est_table)

  cat("\nGoodness-of-Fit Statistics:\n")
  cat(sprintf("  Log-Likelihood: %.4f\n", x$logLik))
  cat(sprintf("  AIC: %.4f\n", x$AIC))
  cat(sprintf("  BIC: %.4f\n", x$BIC))
  cat(sprintf("  Pearson Chi-sq: %.4f\n", x$pearson_chisq))
  cat(sprintf("  Wald Chi-sq: %.4f\n", x$wald_chisq))

  cat("\nConvergence:", ifelse(x$converged, "YES", "NO"), "\n")
  invisible(x)
}
#' Summary method for ctpfit objects
#' @param object A ctpfit object
#' @param ... Additional arguments (ignored)
#' @return
#' Invisibly returns the original fitted model object.
#' The function is called for its side effects, producing
#' a formatted summary of parameter estimates, moments,
#' and goodness-of-fit diagnostics.
#' @export
summary.ctpfit <- function(object, ...) {
  cat("\n=== CTP Distribution Fit Summary ===\n\n")

  print(object)

  cat("\n--- Theoretical Moments ---\n")
  a <- object$estimates["a"]
  b <- object$estimates["b"]
  gama <- object$estimates["gama"]

  mu_theory <- (a^2 + b^2) / (gama - 2*a - 1)
  var_theory <- mu_theory * (mu_theory + gama - 1) / (gama - 2*a - 2)

  cat(sprintf("  Mean: %.4f (empirical: %.4f)\n", mu_theory, mean(object$data)))
  cat(sprintf("  Variance: %.4f (empirical: %.4f)\n", var_theory, var(object$data)))

  cat("\n--- Observed vs Expected Frequencies ---\n")
  print(object$fitted_freq)

  invisible(object)
}
#' Summary method for zictpfit objects
#' @param object A zictpfit object
#' @param ... Additional arguments (ignored)
#' @return
#' Invisibly returns the original fitted model object.
#' The function is called for its side effects, producing
#' a formatted summary of parameter estimates, moments,
#' and goodness-of-fit diagnostics.
#' @export
summary.zictpfit <- function(object, ...) {
  cat("\n=== Zero-Modified CTP Distribution Fit Summary ===\n\n")

  print(object)

  cat("\n--- Moments ---\n")
  a <- object$estimates["a"]
  b <- object$estimates["b"]
  gama <- object$estimates["gama"]
  omega <- object$estimates["omega"]

  mu_ctp <- (a^2 + b^2) / (gama - 2*a - 1)
  var_ctp <- mu_ctp * (mu_ctp + gama - 1) / (gama - 2*a - 2)

  p_ctp_0 <- dctp(0, a, b, gama)

  mu_zmctp <- (1 - omega) * mu_ctp
  var_zmctp <- (1 - omega) * (var_ctp + mu_ctp^2) - mu_zmctp^2

  cat(sprintf("  Mean: %.4f (empirical: %.4f)\n", mu_zmctp, mean(object$data)))
  cat(sprintf("  Variance: %.4f (empirical: %.4f)\n", var_zmctp, var(object$data)))
  cat(sprintf("  P(X=0): %.4f (empirical: %.4f)\n",
              omega + (1-omega)*p_ctp_0, mean(object$data == 0)))

  cat("\n--- Observed vs Expected Frequencies ---\n")
  print(object$fitted_freq)

  invisible(object)
}
