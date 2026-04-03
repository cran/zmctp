#' Plot method for zictpfit objects
#'
#' @description
#' Creates diagnostic plots for Zero-Modified CTP distribution fits,
#' including frequency comparisons, CDF plots, and Q-Q plots.
#'
#' @param x A zictpfit object from zictp.fit()
#' @param type Type of plot: "frequency", "cdf", or "qq"
#' @param ... Additional graphical parameters
#'
#' @return No return value. Called for its side effect of producing a plot.
#' @export
plot.zictpfit <- function(x, type = c("frequency", "cdf", "qq"), ...) {
  type <- match.arg(type)
  
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(try(graphics::par(oldpar), silent = TRUE), add = TRUE)
  graphics::par(mar = c(4.1, 4.1, 2.1, 1.1))
  
  if (type == "frequency") {
    freq_df <- x$fitted_freq
    
    if (is.null(freq_df) || nrow(freq_df) == 0) {
      stop("No frequency data available for plotting.")
    }
    if (!all(c("x", "Freq", "expected") %in% names(freq_df))) {
      stop("Frequency dataframe is missing required columns.")
    }
    
    freq_df <- freq_df[
      is.finite(freq_df$x) &
        is.finite(freq_df$Freq) &
        is.finite(freq_df$expected), ,
      drop = FALSE
    ]
    
    if (nrow(freq_df) == 0) {
      stop("No valid frequency data after removing NA/Inf values.")
    }
    
    counts_mat <- rbind(
      Observed = freq_df$Freq,
      Expected = freq_df$expected
    )
    colnames(counts_mat) <- as.character(freq_df$x)
    
    graphics::barplot(
      counts_mat,
      beside = TRUE,
      col = c("steelblue", "coral"),
      main = "ZM-CTP Fit: Observed vs Expected Frequencies",
      xlab = "Count",
      ylab = "Frequency",
      legend.text = TRUE,
      args.legend = list(x = "topright", bty = "n"),
      ...
    )
    
  } else if (type == "cdf") {
    if (length(x$data) == 0) {
      stop("No data available for CDF plot.")
    }
    
    est <- x$estimates
    a_hat <- unname(est["a"])
    b_hat <- unname(est["b"])
    gama_hat <- unname(est["gama"])
    omega_hat <- unname(est["omega"])
    
    x_range <- 0:max(x$data, na.rm = TRUE)
    ecdf_vals <- vapply(x_range, function(q) mean(x$data <= q, na.rm = TRUE), numeric(1))
    
    tcdf_vals <- tryCatch(
      pzictp(x_range, a_hat, b_hat, gama_hat, omega_hat),
      error = function(e) {
        warning("Error computing theoretical CDF: ", e$message)
        rep(NA_real_, length(x_range))
      }
    )
    
    valid_idx <- is.finite(ecdf_vals) & is.finite(tcdf_vals)
    x_range <- x_range[valid_idx]
    ecdf_vals <- ecdf_vals[valid_idx]
    tcdf_vals <- tcdf_vals[valid_idx]
    
    if (length(x_range) == 0) {
      stop("No valid data for CDF plot.")
    }
    
    graphics::plot(
      x_range, ecdf_vals,
      type = "s",
      lwd = 2,
      col = "steelblue",
      main = "ZM-CTP Fit: Empirical vs Theoretical CDF",
      xlab = "Count",
      ylab = "Cumulative Probability",
      ylim = c(0, 1),
      ...
    )
    graphics::lines(
      x_range, tcdf_vals,
      type = "s",
      lwd = 2,
      col = "coral",
      lty = 2
    )
    graphics::legend(
      "bottomright",
      legend = c("Empirical", "Theoretical"),
      col = c("steelblue", "coral"),
      lwd = 2,
      lty = c(1, 2),
      bty = "n"
    )
    
  } else if (type == "qq") {
    if (length(x$data) == 0) {
      stop("No data available for Q-Q plot.")
    }
    
    est <- x$estimates
    a_hat <- unname(est["a"])
    b_hat <- unname(est["b"])
    gama_hat <- unname(est["gama"])
    omega_hat <- unname(est["omega"])
    
    n <- length(x$data)
    probs <- ((1:n) - 0.5) / n
    empirical_q <- sort(x$data)
    
    theoretical_q <- tryCatch(
      qzictp(probs, a_hat, b_hat, gama_hat, omega_hat),
      error = function(e) {
        warning("Error computing theoretical quantiles: ", e$message)
        rep(NA_real_, length(probs))
      }
    )
    
    valid_idx <- is.finite(theoretical_q)
    empirical_q <- empirical_q[valid_idx]
    theoretical_q <- theoretical_q[valid_idx]
    
    if (length(theoretical_q) == 0) {
      stop("No valid data for Q-Q plot.")
    }
    
    graphics::plot(
      theoretical_q, empirical_q,
      main = "ZM-CTP Fit: Q-Q Plot",
      xlab = "Theoretical Quantiles",
      ylab = "Empirical Quantiles",
      pch = 16,
      col = "steelblue",
      ...
    )
    graphics::abline(0, 1, col = "coral", lwd = 2, lty = 2)
  }
  
  invisible(x)
}