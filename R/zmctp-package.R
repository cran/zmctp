#' @importFrom stats plogis qlogis optim runif var sd
#' @importFrom graphics abline barplot legend lines par plot
NULL

#' zmctp: Zero-Modified Complex Triparametric Pearson Distribution
#'
#' @description
#' The zmctp package extends the Complex Triparametric Pearson (CTP)
#' distribution with zero-modified versions for handling overdispersed
#' count data, particularly when the parameter b approaches zero.
#'
#' @details
#' Main functions:
#' \itemize{
#'   \item \code{\link{dctp}}, \code{\link{pctp}}, \code{\link{qctp}},
#'         \code{\link{rctp}} - CTP distribution functions
#'   \item \code{\link{dzictp}}, \code{\link{pzictp}}, \code{\link{qzictp}},
#'         \code{\link{rzictp}} - Zero-Modified CTP distribution functions
#'   \item \code{\link{ctp.fit}} - Fit CTP model
#'   \item \code{\link{zictp.fit}} - Fit Zero-Modified CTP model
#' }
#'
"_PACKAGE"
