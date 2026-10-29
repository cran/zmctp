#' Maximum likelihood fitting for CTP
#'
#' Fit an independent-observation CTP model in the finite-variance domain.
#' @param x Nonempty numeric vector of finite nonnegative integer counts,
#'   with at least one positive count and maximum at most 1,000,000.
#' @param a_start,b_start,gama_start Optional finite scalar starting values.
#'   A supplied triple must have \code{b_start > 0}, \code{gama_start > 0},
#'   and \code{gama_start > 2*a_start+2}. A partial triple supplements safe
#'   defaults. An inconsistent gamma margin is replaced with a safe margin
#'   with a warning, as in version 0.1.3.
#' @param method An \code{optim} method used to supplement the default
#'   multistart Nelder-Mead search and \code{nlminb} polishing. Supported values
#'   are \code{"L-BFGS-B"}, \code{"BFGS"}, \code{"Nelder-Mead"}, and \code{"CG"}.
#' @param penalty Positive finite objective penalty, retained for compatibility.
#'   Values too small to exceed the initial feasible objectives are increased.
#' @details Coordinates are \eqn{(a,\log b,\log(\gamma-2a-2))}.
#'   The final transformed coordinates are retained. Six fixed starting triples
#'   are used, with a user-supplied start added when present. The two best
#'   finite Nelder-Mead solutions are polished with \code{nlminb} and the
#'   requested supplementary method. Convergence requires a successful
#'   optimiser within \code{1e-4} log-likelihood units of the best finite point.
#'
#'   Wald standard errors require convergence, \code{b >= 1e-6}, variance
#'   margin at least \code{1e-3}, and a positive definite observed Hessian
#'   with minimum eigenvalue greater than \code{1e-10} times its maximum.
#'   Missing uncertainty is reported as \code{NULL}, with an explanation in
#'   \code{diagnostics}. No Hessian jitter or negative-variance clipping is used.
#'
#'   The frequency table includes unobserved counts through the sample maximum
#'   and a final upper-tail bin. Pearson and the historically named
#'   \code{wald_chisq} statistics are descriptive discrepancies, not calibrated
#'   tests. The latter is the sum of squared standardised observations.
#' @return An object of class \code{ctpfit}, retaining fields \code{estimates},
#'   \code{se}, \code{vcov}, \code{logLik}, \code{AIC}, \code{BIC},
#'   \code{pearson_chisq}, \code{wald_chisq}, \code{fitted_freq}, \code{data},
#'   \code{converged}, \code{near_boundary}, \code{method}, and \code{raw_optim}.
#'   Additional fields \code{theta}, \code{attempts}, and \code{diagnostics}
#'   record transformed coordinates, search outcomes and numerical flags.
#'   The \code{raw_optim$par} coordinates are now \code{a}, \code{log_b}, and
#'   \code{log_margin}; they should not be interpreted as natural parameters.
#' @examples
#' x <- rep(0:4, c(30, 16, 8, 3, 1))
#' fit <- ctp.fit(x)
#' fit$estimates
#' fit$diagnostics
#' @export
ctp.fit <- function(x, a_start = NULL, b_start = NULL, gama_start = NULL,
                    method = "L-BFGS-B", penalty = 1e10) {
  .ctp_fit_public(x, a_start, b_start, gama_start, method, penalty, FALSE)
}

#' Maximum likelihood fitting for zero-inflated CTP
#'
#' Fits the zero-inflated model while profiling the mixing probability
#' analytically. The exact parent boundary \code{omega = 0} is permitted.
#' @inheritParams ctp.fit
#' @param omega_start Retained for call compatibility. If supplied it must be
#'   a scalar in \code{[0, 1)}. It does not affect the analytic profile.
#' @details For a fixed parent with zero probability \eqn{f_0} and observed zero
#'   fraction \eqn{q}, the profile estimate is zero when \eqn{q\le f_0} and
#'   \eqn{(q-f_0)/(1-f_0)} otherwise. A parent fit is included as a feasible
#'   nested point. All-zero samples are rejected because the parent parameters
#'   are unidentified. In addition to the checks in \code{ctp.fit}, Wald
#'   uncertainty is withheld within \code{1e-5} of an omega endpoint.
#'   Bootstrap tests remain in the research supplement in this release.
#' @return An object of class \code{zictpfit} with the fields described in
#'   \code{ctp.fit}, plus \code{diagnostics$parent_logLik}.
#' @examples
#' x <- rep(0:4, c(40, 12, 6, 2, 1))
#' fit <- zictp.fit(x)
#' fit$estimates
#' @export
zictp.fit <- function(x, a_start = NULL, b_start = NULL, gama_start = NULL,
                      omega_start = NULL, method = "BFGS") {
  if (!is.null(omega_start)) {
    .ctp_omega(omega_start)
    if (omega_start == 1) stop("omega_start must be below one.", call. = FALSE)
  }
  .ctp_fit_public(x, a_start, b_start, gama_start, method, 1e10, TRUE)
}

.ctp_problem <- function(x) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)) ||
      any(x < 0 | x != floor(x)) || max(x) > 1000000)
    stop("x must be finite nonnegative integer counts, maximum <= 1000000.", call. = FALSE)
  if (all(x == 0))
    stop("All-zero samples do not identify the parent parameters.", call. = FALSE)
  freq <- tabulate(as.integer(x) + 1L, nbins = max(x) + 1L)
  support <- seq_along(freq) - 1L; positive <- which(support > 0 & freq > 0)
  n0 <- freq[1L]; n <- length(x)
  evaluate <- function(theta, omega = NULL) tryCatch({
    a <- unname(theta[1L]); b <- exp(unname(theta[2L]))
    gama <- 2*a + 2 + exp(unname(theta[3L]))
    if (any(!is.finite(c(a,b,gama))) || b <= 0 || gama <= max(0,2*a+2))
      stop("Outside fitting domain.")
    lp <- .ctp_logpmf(support,a,b,gama)
    if (lp[1L] >= 0) stop("Degenerate parent.")
    if (is.null(omega)) {
      p0 <- exp(lp[1L]); q <- n0/n
      omega <- if (q <= p0) 0 else (q-p0)/(-expm1(lp[1L]))
    }
    if (!is.finite(omega) || omega < 0 || omega >= 1) stop("Invalid omega.")
    lzero <- if (omega == 0) lp[1L] else
      .ctp_logadd(log(omega),log1p(-omega)+lp[1L])
    ll <- (if (n0 > 0) n0*lzero else 0) +
      sum(freq[positive]*(log1p(-omega)+lp[positive]))
    if (!is.finite(ll)) stop("Nonfinite likelihood.")
    list(logLik=ll,estimates=c(a=a,b=b,gama=gama,omega=unname(omega)))
  }, error=function(e) list(logLik=-Inf,estimates=rep(NA_real_,4)))
  list(evaluate=evaluate,n=n,zero_fraction=n0/n)
}
.ctp_theta <- function(p) {
  c(a=unname(p[1]),log_b=log(unname(p[2])),
    log_margin=log(unname(p[3]-2*p[1]-2)))
}
.ctp_search <- function(problem, omega, starts, method, penalty) {
  starts <- unique(do.call(rbind,starts))
  initial <- apply(starts,1,function(s) -problem$evaluate(s,omega)$logLik)
  if (!any(is.finite(initial))) stop("No feasible initial likelihood.", call. = FALSE)
  penalty <- max(penalty,10*max(1,initial[is.finite(initial)]))
  objective <- function(theta) {
    ll <- problem$evaluate(theta,omega)$logLik
    if (is.finite(ll)) -ll else penalty
  }
  rows <- list()
  add <- function(theta,engine,stage) {
    warnings <- character()
    z <- tryCatch(withCallingHandlers({
      if (engine == "nlminb") stats::nlminb(theta,objective,
        control=list(iter.max=3000L,eval.max=6000L,rel.tol=1e-12,x.tol=1e-9)) else
        stats::optim(theta,objective,method=engine,
          control=if (engine == "L-BFGS-B") list(maxit=5000L,factr=10) else
            list(maxit=5000L,reltol=1e-12))
    },warning=function(w) {
      warnings <<- unique(c(warnings,conditionMessage(w))); invokeRestart("muffleWarning")
    }),error=function(e)e)
    failed <- inherits(z,"error")
    par <- if (failed) rep(NA_real_,3) else unname(z$par)
    e <- problem$evaluate(par,omega)
    rows[[length(rows)+1L]] <<- data.frame(stage=stage,method=engine,
      a=par[1],log_b=par[2],log_margin=par[3],logLik=e$logLik,
      convergence=if (failed) 99L else z$convergence,
      message=if (failed) conditionMessage(z) else paste(z$message,collapse=" | "),
      warnings=paste(warnings,collapse=" | "),stringsAsFactors=FALSE)
  }
  for (i in seq_len(nrow(starts))) add(starts[i,],"Nelder-Mead","fixed starts")
  tab <- do.call(rbind,rows); valid <- which(is.finite(tab$logLik))
  if (!length(valid)) stop("No finite fit from any start.",call.=FALSE)
  top <- head(valid[order(tab$logLik[valid],decreasing=TRUE)],2L)
  for (i in top) {
    s <- as.numeric(tab[i,c("a","log_b","log_margin")])
    add(s,"nlminb","polish")
    if (method != "Nelder-Mead") add(s,method,"supplementary method")
  }
  # Retain feasible nested/user points without claiming optimiser convergence.
  for (i in seq_len(nrow(starts))) if (is.finite(initial[i])) {
    rows[[length(rows)+1L]] <- data.frame(stage="retained feasible point",method="retained",
      a=starts[i,1],log_b=starts[i,2],log_margin=starts[i,3],logLik=-initial[i],
      convergence=NA_integer_,message="Feasible point; no convergence claim",
      warnings="",stringsAsFactors=FALSE)
  }
  tab <- do.call(rbind,rows); best <- which.max(tab$logLik)
  theta <- setNames(as.numeric(tab[best,c("a","log_b","log_margin")]),
                    c("a","log_b","log_margin"))
  e <- problem$evaluate(theta,omega)
  near <- is.finite(tab$logLik) & tab$logLik >= e$logLik-1e-4
  list(theta=theta,estimates=e$estimates,logLik=e$logLik,attempts=tab,
       convergence=tab$convergence[best],
       converged=any(near & !is.na(tab$convergence) & tab$convergence==0L))
}
.ctp_uncertainty <- function(problem,fit,zi) {
  p <- fit$estimates; margin <- p["gama"]-2*p["a"]-2
  unavailable <- function(status) list(se=NULL,vcov=NULL,hessian=NULL,status=status)
  if (!fit$converged) return(unavailable("Fit did not meet convergence criterion"))
  if (p["b"] < 1e-6 || margin < 1e-3 ||
      (zi && min(p["omega"],1-p["omega"]) < 1e-5))
    return(unavailable("Near boundary; Wald uncertainty withheld"))
  theta <- fit$theta
  if (zi) theta <- c(theta,logit_omega=stats::qlogis(p["omega"]))
  objective <- function(t) {
    ll <- problem$evaluate(t[1:3],if (zi) stats::plogis(t[4]) else 0)$logLik
    if (is.finite(ll)) -ll else 1e100
  }
  H <- tryCatch(stats::optimHess(theta,objective),error=function(e)NULL)
  if (is.null(H) || any(!is.finite(H))) return(unavailable("Hessian unavailable"))
  H <- (H+t(H))/2
  values <- eigen(H,symmetric=TRUE,only.values=TRUE)$values
  if (max(values)<=0 || min(values)<=1e-10*max(values))
    return(unavailable("Hessian non-positive or poorly conditioned"))
  J <- diag(c(1,unname(p["b"]),unname(margin),
              if (zi) unname(p["omega"]*(1-p["omega"]))))
  J[3,1] <- 2
  V <- tryCatch(J %*% chol2inv(chol(H)) %*% t(J),error=function(e)NULL)
  if (is.null(V) || any(!is.finite(V)) || any(diag(V)<0))
    return(unavailable("Covariance unavailable"))
  names_p <- c("a","b","gama",if (zi) "omega")
  dimnames(V) <- list(names_p,names_p)
  list(se=setNames(sqrt(diag(V)),names_p),vcov=V,hessian=H,
       status="Local Wald approximation available")
}
.ctp_fit_public <- function(x,a_start,b_start,gama_start,method,penalty,zi) {
  problem <- .ctp_problem(x)
  method <- match.arg(method,c("L-BFGS-B","BFGS","Nelder-Mead","CG"))
  .ctp_scalar(penalty,"penalty")
  if (penalty <= 0) stop("penalty must be positive.",call.=FALSE)
  grid <- rbind(c(-2,1.5,2),c(-1,1.5,4),c(-0.5,1,4),
                c(0.5,0.5,5),c(2,1,10),c(4,0.5,20))
  starts <- lapply(seq_len(nrow(grid)),function(i).ctp_theta(grid[i,]))
  supplied <- list(a_start,b_start,gama_start)
  if (any(!vapply(supplied,is.null,logical(1)))) {
    p <- c(a=0.5,b=0.5,gama=5)
    for (i in seq_len(3L)) if (!is.null(supplied[[i]])) {
      .ctp_scalar(supplied[[i]],c("a_start","b_start","gama_start")[i])
      p[i] <- supplied[[i]]
    }
    if (p[2] <= 0 || p[3] <= 0) stop("Starting b and gama must be positive.",call.=FALSE)
    if (p[3]-2*p[1]-2 <= 1e-4) {
      warning("Default gama_start is inconsistent with a_start; using a safe starting margin.",call.=FALSE)
      p[3] <- max(0,2*p[1]+2)+2
    }
    starts <- c(list(.ctp_theta(p)),starts)
  }
  parent <- .ctp_search(problem,0,starts,method,penalty)
  fit <- if (zi) .ctp_search(problem,NULL,c(starts,list(parent$theta)),method,penalty) else parent
  p <- fit$estimates
  u <- .ctp_uncertainty(problem,fit,zi)
  margin <- unname(p["gama"]-2*p["a"]-2)
  support <- 0:max(x); observed <- tabulate(as.integer(x)+1L,nbins=length(support))
  probabilities <- if (zi) dzictp(support,p[1],p[2],p[3],p[4]) else dctp(support,p[1],p[2],p[3])
  upper <- if (zi) pzictp(max(x),p[1],p[2],p[3],p[4],lower.tail=FALSE) else
    pctp(max(x),p[1],p[2],p[3],lower.tail=FALSE)
  freq <- data.frame(x=c(support,max(x)+1),Freq=c(observed,0L),
    expected=length(x)*c(probabilities,upper),tail=c(rep(FALSE,length(support)),TRUE))
  expected <- freq$expected; valid <- expected > 0
  pearson <- if (any(expected==0 & freq$Freq>0)) Inf else
    sum((freq$Freq[valid]-expected[valid])^2/expected[valid])
  mu <- if (zi) mean_zictp(p[1],p[2],p[3],p[4]) else mean_ctp(p[1],p[2],p[3])
  variance <- if (zi) var_zictp(p[1],p[2],p[3],p[4]) else var_ctp(p[1],p[2],p[3])
  k <- if (zi) 4L else 3L
  near <- margin < 1e-3
  if (near) warning("Estimated gama is very close to the variance-existence boundary; inspect diagnostics.",call.=FALSE)
  result <- list(estimates=p[seq_len(k)],se=u$se,vcov=u$vcov,logLik=fit$logLik,
    AIC=-2*fit$logLik+2*k,BIC=-2*fit$logLik+k*log(length(x)),
    pearson_chisq=pearson,wald_chisq=if (is.finite(variance) && variance>0)
      sum((x-mu)^2/variance) else NA_real_,fitted_freq=freq,data=x,
    converged=fit$converged,near_boundary=near,method=method,
    raw_optim=list(par=fit$theta,value=-fit$logLik,convergence=fit$convergence,
                   hessian=u$hessian,parameterisation="a, log_b, log_margin"),
    theta=fit$theta,attempts=fit$attempts,
    diagnostics=list(uncertainty=u$status,variance_margin=margin,
      b_near_zero=unname(p["b"]<1e-6),gama_near_zero=unname(p["gama"]<1e-6),
      omega_boundary=if (zi) unname(p["omega"]==0) else NULL,
      parent_logLik=parent$logLik,parent_converged=parent$converged,
      fit_domain="b > 0; gama > max(0, 2*a+2)",
      goodness_of_fit="Descriptive discrepancies; no calibrated p-values"))
  class(result) <- if (zi) "zictpfit" else "ctpfit"
  result
}
