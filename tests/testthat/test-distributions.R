test_that("the analytic heavy-tail example is normalised without truncation", {
  # For a=1,b=0,gama=3, p(k)=1/((k+1)*(k+2)), S(k)=1/(k+2).
  x <- c(0,1,2,100,1000,10000)
  expect_equal(dctp(x,1,0,3),1/((x+1)*(x+2)),tolerance=1e-12)
  expect_equal(pctp(x,1,0,3),1-1/(x+2),tolerance=1e-12)
  expect_equal(pctp(x,1,0,3,lower.tail=FALSE),1/(x+2),tolerance=1e-11)
  expect_equal(dctp(0,1,0,3),dctp(c(0,10000),1,0,3)[1],tolerance=1e-14)
  expect_equal(dctp(1000000,1,0,3,log=TRUE),-log(1000001)-log(1000002),tolerance=1e-10)
})

test_that("quantiles expand, obey inversion, and never silently cap a count", {
  expect_equal(qctp(0.9997,1,0,3),3332)
  expect_equal(qzictp(0.9997,1,0,3,0),3332)
  expect_error(qctp(0.9997,1,0,3,max_count=1000),"not bracketed")
  expect_equal(qctp(c(0,1),1,0,3),c(0,Inf))
  p <- c(.01,.1,.5,.9,.99,.9997)
  for (w in c(0,.3,1)) {
    q <- qzictp(p,1,0,3,w)
    expect_true(all(pzictp(q,1,0,3,w) >= p-1e-12))
    expect_true(all(pzictp(q-1,1,0,3,w) < p+1e-12))
    expect_equal(qzictp(log(p),1,0,3,w,log.p=TRUE),q)
    expect_equal(qzictp(log1p(-p),1,0,3,w,lower.tail=FALSE,log.p=TRUE),q)
  }
})

test_that("mixture endpoints, terminating parents and missing values behave correctly", {
  x <- 0:10
  expect_equal(dzictp(x,1,.5,6,0),dctp(x,1,.5,6))
  expect_equal(dzictp(x,1,.5,6,1),as.numeric(x==0))
  expect_equal(rzictp(10,1,.5,6,1),rep(0L,10))
  expect_equal(pzictp(c(-Inf,-1,0,Inf,NA),1,.5,6,1),c(0,0,1,1,NA))
  expect_equal(qctp(c(0,.5,1),0,0,1),c(0,0,0))
  expect_equal(qctp(1,-2,0,.5),2)
  expect_equal(sum(dctp(0:2,-2,0,.5)),1,tolerance=1e-12)
  expect_true(is.finite(skew_ctp(-1,0,.5)))
  expect_equal(dctp(c(-Inf,-1,.1,Inf,NA),1,.5,6),c(0,0,0,0,NA))
  expect_equal(pctp(c(-Inf,0,Inf,NA),1,.5,6,lower.tail=FALSE),
               c(1,1-dctp(0,1,.5,6),0,NA),tolerance=1e-12)
  expect_equal(qctp(c(NA,0,1),1,.5,6),c(NA,0,Inf))
  expect_length(dctp(numeric(),1,.5,6),0)
  expect_length(rctp(0,1,.5,6),0)
  expect_error(rctp(1.2,1,.5,6),"integer")
  expect_error(dctp(0,2,.5,3),"Require")
  expect_error(dzictp(0,1,.5,6,-.1),"omega")
  expect_warning(expect_equal(mean_ctp(1,.5,2.5),Inf),"moment")
  expect_warning(expect_equal(var_ctp(1,.5,3.5),Inf),"moment")
  expect_warning(expect_true(is.nan(skew_ctp(1,.5,4.5))),"moment")
})

test_that("PMFs agree with independent 75-digit complex gamma references", {
  r <- read.csv(test_path("fixtures_pmf.csv"))
  actual <- mapply(function(a,b,gama,x) dctp(x,a,b,gama,log=TRUE),
                   r$a,r$b,r$gama,r$x)
  expect_equal(actual,r$log_probability,tolerance=2e-9)
})
