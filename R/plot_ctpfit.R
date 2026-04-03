#' Plot method for ctpfit objects
#'
#' @param x A ctpfit object
#' @param type Type of plot: "frequency", "cdf", or "qq"
#' @param ... Additional graphical parameters
#'
#' @return No return value. Called for its side effect of producing a plot.
#' @export
plot.ctpfit <- function(x, type = c("frequency", "cdf", "qq"), ...) {
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
      stop("Frequency data frame is missing required columns.")
    }
    
    freq_df <- freq_df[
      is.finite(freq_df$x) &
        is.finite(freq_df$Freq) &
        is.finite(freq_df$expected), ,
      drop = FALSE
    ]
    
    if (nrow(freq_df) == 0) {
      stop("No valid frequency data available for plotting.")
    }
    
    graphics::barplot(
      rbind(Observed = freq_df$Freq, Expected = freq_df$expected),
      beside = TRUE,
      col = c("steelblue", "coral"),
      main = "CTP Fit: Observed vs Expected Frequencies",
      xlab = "Count",
      ylab = "Frequency",
      legend.text = TRUE,
      args.legend = list(x = "topright", bty = "n"),
      ...
    )
    
  } else if (type == "cdf") {
    if (is.null(x$data) || length(x$data) == 0) {
      stop("No data available for CDF plot.")
    }
    
    est <- x$estimates
    if (is.null(est) || !all(c("a", "b", "gama") %in% names(est))) {
      stop("Missing parameter estimates needed for CDF plot.")
    }
    
    x_range   <- 0:max(x$data, na.rm = TRUE)
    ecdf_vals <- vapply(x_range, function(q) mean(x$data <= q, na.rm = TRUE), numeric(1))
    tcdf_vals <- pctp(x_range, est["a"], est["b"], est["gama"])
    
    valid_idx  <- is.finite(ecdf_vals) & is.finite(tcdf_vals)
    x_range    <- x_range[valid_idx]
    ecdf_vals  <- ecdf_vals[valid_idx]
    tcdf_vals  <- tcdf_vals[valid_idx]
    
    if (length(x_range) == 0) {
      stop("No valid data for CDF plot.")
    }
    
    graphics::plot(
      x_range, ecdf_vals,
      type = "s", lwd = 2, col = "steelblue",
      main = "CTP Fit: Empirical vs Theoretical CDF",
      xlab = "Count", ylab = "Cumulative Probability",
      ylim = c(0, 1),
      ...
    )
    graphics::lines(
      x_range, tcdf_vals,
      type = "s", lwd = 2, col = "coral", lty = 2
    )
    graphics::legend(
      "bottomright",
      legend = c("Empirical", "Theoretical"),
      col = c("steelblue", "coral"),
      lwd = 2, lty = c(1, 2), bty = "n"
    )
    
  } else if (type == "qq") {
    if (is.null(x$data) || length(x$data) == 0) {
      stop("No data available for Q-Q plot.")
    }
    
    est <- x$estimates
    if (is.null(est) || !all(c("a", "b", "gama") %in% names(est))) {
      stop("Missing parameter estimates needed for Q-Q plot.")
    }
    
    n             <- length(x$data)
    probs         <- ((1:n) - 0.5) / n
    empirical_q   <- sort(x$data)
    theoretical_q <- qctp(probs, est["a"], est["b"], est["gama"])
    
    valid_idx     <- is.finite(theoretical_q) & is.finite(empirical_q)
    theoretical_q <- theoretical_q[valid_idx]
    empirical_q   <- empirical_q[valid_idx]
    
    if (length(theoretical_q) == 0) {
      stop("No valid data for Q-Q plot.")
    }
    
    graphics::plot(
      theoretical_q, empirical_q,
      main = "CTP Fit: Q-Q Plot",
      xlab = "Theoretical Quantiles",
      ylab = "Empirical Quantiles",
      pch = 16, col = "steelblue",
      ...
    )
    graphics::abline(0, 1, col = "coral", lwd = 2, lty = 2)
  }
  
  invisible(x)
}