## ----listsolvers--------------------------------------------------------------
library(optimx)
checkallsolvers()

## ----sqmod--------------------------------------------------------------------
library(optimx)
sqmod<-function(z, xx){
   nn<-length(z)
   yy<-xx^(1:nn)
   f<-sum((yy-z)^2)
   f
}
sqmod.g <- function(z, xx){
   nn<-length(z)
   yy<-xx^(1:nn)
   gg<- 2*(z - yy)
}

z0 <- rep(0.1, 2) # initial parameters
sol <- optimr(par=z0, fn=sqmod, gr=sqmod.g, method="nvm", xx=1.5)
proptimr(sol)

## ----sqmod2-------------------------------------------------------------------
solm <- opm(par=z0, fn=sqmod, gr=sqmod.g, method="MOST", xx=1.5)
cat("Result of applying opm() to ",attr(solm,"fname"))
print(summary(solm, order=value))

## ----polyoptx1----------------------------------------------------------------
fnR <- function (x, gs=100.0) 
{
    n <- length(x)
    x1 <- x[2:n]
    x2 <- x[1:(n - 1)]
    sum(gs * (x1 - x2^2)^2 + (1 - x2)^2)
}
grR <- function (x, gs=100.0) 
{
    n <- length(x)
    g <- rep(NA, n)
    g[1] <- 2 * (x[1] - 1) + 4*gs * x[1] * (x[1]^2 - x[2])
    if (n > 2) {
        ii <- 2:(n - 1)
        g[ii] <- 2 * (x[ii] - 1) + 4 * gs * x[ii] * (x[ii]^2 - x[ii + 
            1]) + 2 * gs * (x[ii] - x[ii - 1]^2)
    }
    g[n] <- 2 * gs * (x[n] - x[n - 1]^2)
    g
}

x0 <- rep(pi, 4)
mc <- data.frame(method=c("Nelder-Mead","Rvmmin"), maxit=c(1000, 100), maxfeval= c(1000, 1000))

ans <- polyopt(x0, fnR, grR, methcontrol=mc, control=list(trace=0))
ans
mc <- data.frame(method=c("Nelder-Mead","Rvmmin"), maxit=c(100, 100), maxfeval= c(100, 1000))

ans <- polyopt(x0, fnR, grR, methcontrol=mc, control=list(trace=0))
ans

mc <- data.frame(method=c("Nelder-Mead","Rvmmin"), maxit=c(10, 100), maxfeval= c(10, 1000))

ans <- polyopt(x0, fnR, grR, methcontrol=mc, control=list(trace=0))
ans

## ----multix1------------------------------------------------------------------
# edited from inst/doc/examples/StyblinskiTang2308.R
require(optimx, quietly=TRUE)
fnStyblinskiTang <- function(x) 0.5 * sum(x^4 - 16*x^2 + 5*x)
grST <- function(x) {
   2*x^3 - 16*x + 2.5
}
n <- 4
cat("Globalmin ",n," = ",(-39.16599)*length(x0),"\n")
methlist<-c("hjn", "nvm", "Nelder-Mead")
stmat<-matrix(runif(6*n, -5, 5), nrow=6, ncol=n)
cat("Starting parameters (rows of matrix):\n")
print(stmat)
for (mstr in methlist){
  cat("Method ",mstr,":\n")
  mtst <- multistart(stmat, fnStyblinskiTang, grST, method=mstr)
  print(mtst)
}
# Try a single start
tst1<-opm(rep(-5, n), fnStyblinskiTang,grST,method=methlist)
summary(tst1)

## ----pracmanm1----------------------------------------------------------------
library(optimx)
fnR <- function (x, gs=100.0) 
{
  n <- length(x)
  x1 <- x[2:n]
  x2 <- x[1:(n - 1)]
  sum(gs * (x1 - x2^2)^2 + (1 - x2)^2)
}
x0 <- rep(pi, 4)

apnm0<-optimr(x0, fnR, method="pracmanm")
proptimr(apnm0)
apnm1<-optimr(x0, fnR, method="pracmanm", control=list(pracmanmtol=1e-13))
proptimr(apnm1)

## ----diffderivs1--------------------------------------------------------------
mymeth<-"ncg"
# using fnR and grR from earlier chunk
system.time(sola <- optimr(par=x0, fn=fnR, gr=grR, method=mymeth))
proptimr(sola)
system.time(solfwd <- optimr(par=x0, fn=fnR, gr="grfwd", method=mymeth))
proptimr(solfwd)
system.time(solback <- optimr(par=x0, fn=fnR, gr="grback", method=mymeth))
proptimr(solback)
system.time(solctl <- optimr(par=x0, fn=fnR, gr="grcentral", method=mymeth))
proptimr(solctl)
system.time(solnd <- optimr(par=x0, fn=fnR, gr="grnd", method=mymeth))
proptimr(solnd)
system.time(solpra <- optimr(par=x0, fn=fnR, gr="grpracma", method=mymeth))
proptimr(solpra)

## ----bobspec------------------------------------------------------------------
require(optimx)
hobbs.f<- function(x){ # # Hobbs weeds problem -- function
    if (abs(12*x[3]) > 500) { # check computability
       fbad<-.Machine$double.xmax
       return(fbad)
    }
    res<-hobbs.res(x)
    f<-sum(res*res)
}
hobbs.res<-function(x){ # Hobbs weeds problem -- residual
# This variant uses looping
    if(length(x) != 3) stop("hobbs.res -- parameter vector n!=3")
    y<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558, 50.156, 62.948,
         75.995, 91.972)
    t<-1:12
    if(abs(12*x[3])>50) {
       res<-rep(Inf,12)
    } else {
       res<-100*x[1]/(1+10*x[2]*exp(-0.1*x[3]*t)) - y
    }
}
x0 <- c(1,1,1)
cat("default rhobeg and rhoend in bobyqa\n")
tdef <- optimr(x0, hobbs.f, method="bobyqa")
proptimr(tdef)
tsp3 <- optimr(x0, hobbs.f, method="bobyqa", control=list(rhobeg=10, rhoend=1e-4))
proptimr(tsp3)

