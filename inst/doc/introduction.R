## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)
library(zmctp)

# Flag: skip slow chunks on CRAN (they time out during R CMD build)
NOT_CRAN <- !identical(Sys.getenv("NOT_CRAN"), "true")

## ----basic-dist---------------------------------------------------------------
# Probability mass function
dctp(0:5, a = 1, b = 0.5, gama = 5)

# Cumulative distribution function
pctp(3, a = 1, b = 0.5, gama = 5)

# Quantile function (inverse CDF)
qctp(c(0.25, 0.5, 0.75), a = 1, b = 0.5, gama = 5)

# Random generation
set.seed(123)
x <- rctp(30, a = 1, b = 0.5, gama = 5)
cat("Sample mean:", mean(x), "\nSample variance:", var(x), "\n")

## ----check-dispersion---------------------------------------------------------
# Generate overdispersed data
set.seed(456)
x_over <- rctp(20, a = 1.2, b = 0.3, gama = 6)
head(x_over)
# Calculate dispersion index
dispersion_index <- var(x_over) / mean(x_over)
cat("Dispersion Index:", dispersion_index, "\n")

if (dispersion_index > 1) {
  cat("Data is OVERDISPERSED\n")
} else if (dispersion_index < 1) {
  cat("Data is UNDERDISPERSED\n")
} else {
  cat("Data is EQUIDISPERSED\n")
}

## ----fit-ctp, eval=FALSE------------------------------------------------------
# # Fit CTP model to moderate-sized data
# set.seed(456)
# x_over <- rctp(200, a = 1.2, b = 0.3, gama = 6)
# fit_ctp <- ctp.fit(x_over)
# print(fit_ctp)

## ----fit-ctp-output, echo=FALSE, comment=""-----------------------------------
cat("CTP Distribution - Maximum Likelihood Estimates
===============================================

Parameter Estimates:
     Estimate Std.Error
a      1.4593    0.3585
b      0.0011    0.0000
gama   8.5912    3.1878

Goodness-of-Fit Statistics:
  Log-Likelihood: -178.1003
  AIC: 362.2006
  BIC: 372.0956
  Pearson Chi-sq: 19.6799
  Wald Chi-sq: 189.8432

Convergence: YES")

## ----plot-ctp, fig.height=6, eval=FALSE---------------------------------------
# plot(fit_ctp, type = "frequency")
# plot(fit_ctp, type = "cdf")
# plot(fit_ctp, type = "qq")

## ----zm-ctp-example, eval=FALSE-----------------------------------------------
# # Generate zero-inflated data
# 
# x_zi <- rzictp(300, a = 1, b = 0.5, gama = 6, omega = 0.3)
# cat("Proportion of zeros:", mean(x_zi == 0), "\n")
# cat("Expected P(X=0) under CTP:", dctp(0, 1, 0.5, 6), "\n")
# fit_ctp <- ctp.fit(x_zi)
# fit_zm <- zictp.fit(x_zi)
# summary(fit_zm)
# cat("Standard CTP AIC:", fit_ctp$AIC, "\n")
# cat("ZM-CTP AIC:", fit_zm$AIC, "\n")
# cat("Omega estimate:", fit_zm$estimates["omega"], "\n")

## ----echo=FALSE, comment=""---------------------------------------------------
cat("Proportion of zeros: 0.7966667 
Expected P(X=0) under CTP: 0.7570218 

=== Zero-Modified CTP Distribution Fit Summary ===


Zero-Modified CTP Distribution - Maximum Likelihood Estimates
=============================================================

Parameter Estimates:
      Estimate Std.Error
a       2.5646        NA
b       0.0017        NA
gama   16.8737        NA
omega   0.4656        NA

Goodness-of-Fit Statistics:
  Log-Likelihood: -216.5545
  AIC: 441.1090
  BIC: 455.9241
  Pearson Chi-sq: 2.7242
  Wald Chi-sq: 287.6122

Convergence: YES 

--- Moments ---
  Mean: 0.3271 (empirical: 0.3267)
  Variance: 0.6467 (empirical: 0.6220)
  P(X=0): 0.7967 (empirical: 0.7967)

--- Observed vs Expected Frequencies ---
  x Freq    expected
1 0  239 239.0001656
2 1   40  38.7105351
3 2   11  13.7597381
4 3    5   5.0633781
5 4    4   1.9722958
6 5    1   0.8143692

Standard CTP AIC: 440.03 
ZM-CTP AIC: 441.109 
Omega estimate: 0.4656324")

## ----cpd-comparison, eval=FALSE-----------------------------------------------
# # Install cpd if needed
# # install.packages("cpd")
# 
# library(zmctp)
# 
# # Generate data where cpd struggles
# set.seed(100)
# x_problem <- rzictp(2000, a = 1, b = 0.001, gama = 8, omega = 0.2)
# 
# # Compare results
# cpd_fit <- cpd::fitctp(x_problem, astart=1, bstart=0.001, gammastart=8)
# zmctp_fit <- zictp.fit(x_problem)
# cat("cpd b estimate:", cpd_fit$coefficients[2], "\n")
# cat("zmctp b estimate:", zmctp_fit$estimates["b"], "\n")
# cat("zmctp omega estimate:", zmctp_fit$estimates["omega"], "\n")

## ----cpd-comparison-result, echo=FALSE, comment=""----------------------------
cat("
cpd b estimate: 8.072796e-07 
zmctp b estimate: 0.0007942037 
zmctp omega estimate: 0.04831814")

## ----session-info-------------------------------------------------------------
sessionInfo()

