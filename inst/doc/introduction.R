## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE,comment="#>",fig.width=6,fig.height=4)
library(zmctp)

## ----distributions------------------------------------------------------------
dctp(0:6,a=1,b=0.5,gama=6)
pctp(c(-Inf,0,2,Inf),1,0.5,6)
qctp(c(0,0.5,0.9,1),1,0.5,6)

## ----zero-inflation-----------------------------------------------------------
c(parent=dctp(0,1,0.5,6),mixture=dzictp(0,1,0.5,6,0.1))
mean_zictp(1,0.5,6,0.3)
var_zictp(1,0.5,6,0.3)
set.seed(123)
rzictp(12,1,0.5,6,0.3)

## ----fitting------------------------------------------------------------------
x <- rep(0:4,c(40,12,6,2,1))
fit <- zictp.fit(x)
fit$estimates
fit$logLik
fit$diagnostics

## ----plot---------------------------------------------------------------------
plot(fit,type="frequency")

## ----expanded-search----------------------------------------------------------
# This parent has F(k) = 1 - 1/(k+2).
qctp(0.9997,a=1,b=0,gama=3,max_q=1000)

