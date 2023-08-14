## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----manypar------------------------------------------------------------------
# example using var_dim function from More et al. (or funconstrain)
bfn = function(par) {
  n <- length(par)
  if (n < 1) {stop("Variably Dimensioned: n must be positive")}
  fsum <- 0;  fn1 <- 0;
  for (j in 1:n) {
    fj <- par[j] - 1; fsum <- fsum + fj * fj; fn1 <- fn1 + j * fj
  }
  fn1_fn1 <- fn1 * fn1  # f_n+1 and f_n+2
  fsum <- fsum + fn1_fn1 + fn1_fn1 * fn1_fn1
  fsum
}
bgr = function(par) {
  n <- length(par)
#  if (n < 1) {stop("Variably Dimensioned: n must be positive")  }
  fsum <- 0;  grad <- rep(0, n);  fn1 <- 0
  for (j in 1:n) {
    fj <- par[j] - 1; fn1 <- fn1 + j * fj; grad[j] <- grad[j] + 2 * fj
  }
  fn1_2 <- fn1 * 2; fn13_4 <- fn1_2 * fn1_2 * fn1
  grad <- grad + 1:n * (fn1_2 + fn13_4)
  grad
}
x0<-rep(pi,100)
library(optimx)
defaultW <- getOption("warn")
options(warn = -1)
mm<-c("ncg","Rcgmin","lbfgs","L-BFGS-B","Rtnmin")
res1<-opm(x0, bfn, bgr, method=mm)
options(warn = defaultW)
res1[,100:108]

