## ----code=readLines("SimpHess.R")---------------------------------------------
# SimpHess.R
x0<-c(1,2,3,4)
lo <- x0-0.5
up <- x0+1.0
fnt <- function(x, fscale=10){
  yy <- length(x):1
  val <- sum((yy*x)^2)*fscale
  val
}
attr(fnt,"fname")<-"SimpHess"
grt <- function(x, fscale=10){
  nn <- length(x)
  yy <- nn:1
  #    gg <- rep(NA,nn)
  gg <- 2*(yy^2)*x*fscale
  gg
}
hesst <- function(x, fscale=10){
  nn <- length(x)
  yy <- nn:1
  hh <- diag(2*yy^2*fscale)
  hh
}

require(optimx)
t1 <- optimr(x0, fnt, grt, hesst, method="snewton", control=list(trace=0, usexxxmeth=TRUE), fscale=3.0)
proptimr(t1)
t1m <- optimr(x0, fnt, grt, hesst, method="snewtonm", control=list(trace=0), fscale=3.0)
proptimr(t1m)

t1nlmo <- optimr(x0, fnt, grt, hess=hesst, method="nlm", fscale=3.0, 
                 control=list(trace=0))
proptimr(t1nlmo)

## BUT ... nlminb may not be using a true Newton-type method
tst <- try(t1nlminbo <- optimr(x0, fnt, grt, hess=hesst, method="nlminb", 
                               fscale=3.0, control=list(trace=0)))

# Bounded
tstb <- try(t1nlminbb <- optimr(x0, fnt, grt, hess=hesst, method="nlminb", 
                lower=lo, upper=up, fscale=3.0, control=list(trace=0)))
proptimr(t1nlminbb) 

t1smb <-  optimr(x0, fnt, grt, hess=hesst, method="snewtonm", fscale=3.0, 
                 lower=lo, upper=up, control=list(trace=0))
proptimr(t1smb)

# Masked
lo[1]<-x0[1]
up[1]<-x0[1]
lo[4]<-x0[4]
up[4]<-x0[4]
tstm <- try(t1nlminbm <- optimr(x0, fnt, grt, hess=hesst, method="nlminb", 
                                 lower=lo, upper=up, fscale=3.0, control=list(trace=0)))
proptimr(t1nlminbm) 

t1smm <-  optimr(x0, fnt, grt, hess=hesst, method="snewtonm", fscale=3.0, 
                 lower=lo, upper=up, control=list(trace=0))
proptimr(t1smm)

## ----code=readLines("RosenHess.R")--------------------------------------------
# RosenHess.R
require(optimx)
f <- function(x){ #Rosenbrock banana valley function
  return(100*(x[2] - x[1]*x[1])^2 + (1-x[1])^2)
}
attr(f,"fname")<-"RosenHess"
gr <- function(x){ #gradient
  return(c(-400*x[1]*(x[2] - x[1]*x[1]) - 2*(1-x[1]), 200*(x[2] - x[1]*x[1])))
}
h <- function(x) { #Hessian
  a11 <- 2 - 400*x[2] + 1200*x[1]*x[1]; a21 <- -400*x[1]
  return(matrix(c(a11, a21, a21, 200), 2, 2))
}
x0 <- c(-1.2, 1)
t1 <- snewton(x0, fn=f, gr=gr, hess=h, control=list(trace=0))
proptimr(t1)

# we can also use nlm and nlminb
fght <- function(x){ ## combine f, g and h into single function for nlm
  ff <- f(x)
  gg <- gr(x)
  hh <- h(x)
  attr(ff, "gradient") <- gg
  attr(ff, "hessian") <- hh
  ff
}

t1nlmo <- optimr(x0, f, gr, hess=h, method="nlm", control=list(trace=0))
proptimr(t1nlmo)

t1so <- optimr(x0, f, gr, hess=h, method="snewton", control=list(trace=0))
proptimr(t1so)

t1smo <-  optimr(x0, f, gr, hess=h, method="snewtonm", control=list(trace=0))
proptimr(t1smo)


## nlminb 
tst <- try(t1nlminbo <- optimr(x0, f, gr, hess=h, method="nlminb", 
                               control=list(trace=0)))
if (class(tst) == "try-error"){
  cat("try-error on attempt to run nlminb in optimr()\n")
} else { proptimr(t1nlminbo) }

tstnoh <- try(t1nlminbnoho <- optimr(x0, f, gr, hess=NULL, method="nlminb", 
                               control=list(trace=0)))
if (class(tstnoh) == "try-error"){
  cat("try-error on attempt to run nlminb in optimr()\n")
} else { proptimr(t1nlminbnoho) }


## ----code=readLines("WoodHess.R")---------------------------------------------
#WoodHess.R -- Wood function example
wood.f <- function(x){
  res <- 100*(x[1]^2-x[2])^2+(1-x[1])^2+90*(x[3]^2-x[4])^2+(1-x[3])^2+
    10.1*((1-x[2])^2+(1-x[4])^2)+19.8*(1-x[2])*(1-x[4])
  attr(res,"fname")<-"WoodHess"
  return(res)
}
wood.g <- function(x){ #gradient
  g1 <- 400*x[1]^3-400*x[1]*x[2]+2*x[1]-2
  g2 <- -200*x[1]^2+220.2*x[2]+19.8*x[4]-40
  g3 <- 360*x[3]^3-360*x[3]*x[4]+2*x[3]-2
  g4 <- -180*x[3]^2+200.2*x[4]+19.8*x[2]-40
  return(c(g1,g2,g3,g4))
}
wood.h <- function(x){ #hessian
  h11 <- 1200*x[1]^2-400*x[2]+2;    h12 <- -400*x[1]; h13 <- h14 <- 0
  h22 <- 220.2; h23 <- 0;    h24 <- 19.8
  h33 <- 1080*x[3]^2-360*x[4]+2;    h34 <- -360*x[3]
  h44 <- 200.2
  H <- matrix(c(h11,h12,h13,h14,h12,h22,h23,h24,
                h13,h23,h33,h34,h14,h24,h34,h44),ncol=4)
  return(H)
}
x0 <- c(-3,-1,-3,-1) # Wood standard start

require(optimx) # call methods through optimr() function
# But note that nlm default settings have lower iteration limit 
#  and in 100 "iterations" do not get to solution
t1nlm <- optimr(x0, fn=wood.f, gr=wood.g, hess=wood.h, method="nlm")
proptimr(t1nlm)
# But both optimx Newton approaches do work
wd <- optimr(x0, fn=wood.f, gr=wood.g, hess=wood.h, method="snewton")
proptimr(wd)
wdm <- optimr(x0, fn=wood.f, gr=wood.g, hess=wood.h, method="snewtonm")
proptimr(wdm)
# nlminb 
t1nlminb <- optimr(x0, fn=wood.f, gr=wood.g, hess=wood.h, method="nlminb")
proptimr(t1nlminb)

wood.fgh <- function(x){
  fval <- wood.f(x)
  gval <- wood.g(x)
  hval <- wood.h(x)
  attr(fval,"gradient") <- gval
  attr(fval,"hessian")<- hval
  fval
}

# direct call to nlm
t1nlm <- nlm(wood.fgh, x0, print.level=0)
print(t1nlm)
# Check that optimr gets same result with similar iteration limit of 100
t1nlmo <- optimr(x0, wood.f, wood.g, hess=wood.h, method="nlm", control=list(maxit=100))
proptimr(t1nlmo)
print(wood.g(t1nlmo$par))

# Run for allowed iteration limit in optimr 500*round(sqrt(npar+1)) = 1000
tst<-try(t1nlminbo <- optimr(x0, wood.f, wood.g, hess=wood.h, method="nlminb"))
if (class(tst) == "try-error"){
  cat("try-error on attempt to run nlminb in optimr()\n")
} else { proptimr(t1nlminbo) }


## ----code=readLines("GenRoseHess.R")------------------------------------------
# GenRoseHess.R
# genrosa function code -- attempts to match the rosenbrock at gs=100 and x=c(-1.2,1)
genrosa.f<- function(x, gs=NULL){ # objective function
  ## One generalization of the Rosenbrock banana valley function (n parameters)
  n <- length(x)
  if(is.null(gs)) { gs=100.0 }
  # Note do not at 1.0 so min at 0
  fval<-sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
}
attr(genrosa.f, "fname")<-"genrosa"

genrosa.g <- function(x, gs=NULL){
  # vectorized gradient for genrose.f
  # Ravi Varadhan 2009-04-03
  n <- length(x)
  if(is.null(gs)) { gs=100.0 }
  gg <- as.vector(rep(0, n))
  tn <- 2:n
  tn1 <- tn - 1
  z1 <- x[tn] - x[tn1]^2
  z2 <- 1 - x[tn1]
  # f = gs*z1*z1 + z2*z2
  gg[tn] <- 2 * (gs * z1)
  gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1 - 2 *z2 
  return(gg)
}

genrosa.h <- function(x, gs=NULL) { ## compute Hessian
  if(is.null(gs)) { gs=100.0 }
  n <- length(x)
  hh<-matrix(rep(0, n*n),n,n)
  for (i in 2:n) {
    z1<-x[i]-x[i-1]*x[i-1]
    #		z2<-1.0 - x[i-1]
    hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
    hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
    hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
    hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
  }
  return(hh)
}

require(optimx)
cat("Generalized Rosenbrock tests\n")

cat("original n and x0")

x0 <- c(-1.2, 1)
# solorigs <- snewton(x0, genrosa.f, genrosa.g, genrosa.h) # WORKS OK if optimx loaded
solorig <- optimr(x0, genrosa.f, genrosa.g, genrosa.h, method="snewton", hessian=TRUE)

proptimr(solorig)
print(eigen(solorig$hessian)$values)
solorigm <- optimr(x0, genrosa.f, genrosa.g, genrosa.h, method="snewtonm", hessian=TRUE)
proptimr(solorigm)
print(eigen(solorigm$hessian)$values)

# Start with 50 values of pi and scale factor 10
x0 <- rep(pi, 50)
sol50pi <- optimr(x0, genrosa.f, genrosa.g, genrosa.h, method="snewton", 
                  hessian=TRUE, gs=10)
proptimr(sol50pi)
print(eigen(sol50pi$hessian)$values)
hhi <- genrosa.h(sol50pi$par, gs=10)
print(eigen(hhi)$values)
sol50pim <- optimr(x0, genrosa.f, genrosa.g, genrosa.h, method="snewtonm", 
                   hessian=TRUE, gs=10)
proptimr(sol50pim)
hhm <- genrosa.h(sol50pim$par, gs=10)
print(eigen(hhm)$values)

# Bounds constraints

lo<-rep(3,50)
up<-rep(4,50)
sol50pimb <- optimr(x0, genrosa.f, genrosa.g, genrosa.h, lower=lo, upper=up, method="snewtonm", 
                     hessian=TRUE, gs=10)
proptimr(sol50pimb)

# approximate hessian
solom01 <- optimr(x0, genrosa.f, gr=NULL, hess="approx", method="snewtonm", hessian=TRUE)
proptimr(solom01)
print(eigen(solom01$hessian)$values)
solomg1 <- optimr(x0, genrosa.f, genrosa.g, hess="approx", method="snewtonm", hessian=TRUE)
proptimr(solomg1)
print(eigen(solomg1$hessian)$values)
# Following should fail
solomrr <- try(optimr(x0, genrosa.f, gr=NULL, hess="rubbish", method="snewtonm", hessian=TRUE))


## ----code=readLines("HobbHess.R")---------------------------------------------
# HobbHess.R
## Optimization test function HOBBS
## Nash and Walker-Smith (1987, 1989) ...
require(optimx)

hobbs.f<- function(x){ # # Hobbs weeds problem -- function
  if (abs(12*x[3]) > 500) { # check computability
    fbad<-.Machine$double.xmax
    return(fbad)
  }
  res<-hobbs.res(x)
  f<-sum(res*res)
}
attr(hobbs.f, "fname")<- "Hobbs"

hobbs.res<-function(x){ # Hobbs weeds problem -- residual
  # This variant uses looping
  if(length(x) != 3) stop("hobbs.res -- parameter vector n!=3")
  y<-c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558, 50.156, 62.948,
       75.995, 91.972)
  t<-1:12
  if(abs(12*x[3])>50) {
    res<-rep(Inf,12)
  } else {
    res<-x[1]/(1+x[2]*exp(-x[3]*t)) - y
  }
}

hobbs.jac<-function(x){ # Jacobian of Hobbs weeds problem
  jj<-matrix(0.0, 12, 3)
  t<-1:12
  yy<-exp(-x[3]*t)
  zz<-1.0/(1+x[2]*yy)
  jj[t,1] <- zz
  jj[t,2] <- -x[1]*zz*zz*yy
  jj[t,3] <- x[1]*zz*zz*yy*x[2]*t
  return(jj)
}

hobbs.g<-function(x){ # gradient of Hobbs weeds problem
  # NOT EFFICIENT TO CALL AGAIN
  jj<-hobbs.jac(x)
  res<-hobbs.res(x)
  gg<-as.vector(2.*t(jj) %*% res)
  return(gg)
}


hobbs.rsd<-function(x) { # Jacobian second derivative
  rsd<-array(0.0, c(12,3,3))
  t<-1:12
  yy<-exp(-x[3]*t)
  zz<-1.0/(1+x[2]*yy)
  rsd[t,1,1]<- 0.0
  rsd[t,2,1]<- -yy*zz*zz
  rsd[t,1,2]<- -yy*zz*zz
  rsd[t,2,2]<- 2.0*x[1]*yy*yy*zz*zz*zz
  rsd[t,3,1]<- t*x[2]*yy*zz*zz
  rsd[t,1,3]<- t*x[2]*yy*zz*zz
  rsd[t,3,2]<- t*x[1]*yy*zz*zz*(1-2*x[2]*yy*zz)
  rsd[t,2,3]<- t*x[1]*yy*zz*zz*(1-2*x[2]*yy*zz)
  ##    rsd[t,3,3]<- 2*t*t*x[1]*x[2]*x[2]*yy*yy*zz*zz*zz
  rsd[t,3,3]<- -t*t*x[1]*x[2]*yy*zz*zz*(1-2*yy*zz*x[2])
  return(rsd)
}

hobbs.h <- function(x) { ## compute Hessian
  #   cat("Hessian not yet available\n")
  #   return(NULL)
  H<-matrix(0,3,3)
  res<-hobbs.res(x)
  jj<-hobbs.jac(x)
  rsd<-hobbs.rsd(x)
  ##    H<-2.0*(t(res) %*% rsd + t(jj) %*% jj)
  for (j in 1:3) {
    for (k in 1:3) {
      for (i in 1:12) {
        H[j,k]<-H[j,k]+res[i]*rsd[i,j,k]
      }
    }
  }
  H<-2*(H + t(jj) %*% jj)
  return(H)
}

x0 <- c(200, 50, .3)
cat("Good start for Hobbs:")
print(x0)
solx0 <- optimr(x0, hobbs.f, hobbs.g, hobbs.h, method="snewton", hessian=TRUE)
## Note that we exceed count limit, but have answer
proptimr(solx0)
print(eigen(solx0$hessian)$values)
## Note that we exceed count limit, but have answer

## Setting relative check offset larger gets quicker convergence
solx0a <- optimr(x0, hobbs.f, hobbs.g, hobbs.h, method="snewton", 
                  control=list(offset=1000.))
proptimr(solx0a)


x1s <- c(100, 10, .1)
cat("Scaled start for Hobbs:")
print(x1s)
solx1s <- optimr(x1s, hobbs.f, hobbs.g, hobbs.h, method="snewton", hessian=TRUE , control=list(trace=0))
proptimr(solx1s)
print(eigen(solx1s$hessian)$values)
solx1m <- optimr(x1s, hobbs.f, hobbs.g, hobbs.h, method="snewtonm", hessian=TRUE , control=list(trace=0))
proptimr(solx1m)
print(eigen(solx1m$hessian)$values)

cat("Following test fails in snewton with ERROR \n
     -- Not run as function infinite.\n")
x3 <- c(1, 1, 1)
# solx3 <- try(optimr(x3, hobbs.f, hobbs.g, hobbs.h, method="snewton", control=list(trace=4)))
# if ((solx3$convergence != 0) || class(solx3) != "try-error") {
#   proptimr(solx3)
#   print(eigen(solx3$hessian)$values)
# }
# dirx3 <- try(snewton(x3, hobbs.f, hobbs.g, hobbs.h, control=list(trace=4)))
# if ((dirx3$convergence != 0) || class(dirx3) != "try-error") {
#   proptimr(dirx3)
#   print(eigen(dirx3$Hess)$values)
# }
cat("But Marquardt variant succeeds\n")
solx3m <- optimr(x3, hobbs.f, hobbs.g, hobbs.h, method="snewtonm", 
                  hessian=TRUE, control=list(trace=0))
proptimr(solx3m)
print(eigen(solx3m$hessian)$values)
# we could also use nlm and nlminb and call them from optimr

solx3 <- try(optimr(x3, hobbs.f, hobbs.g, hobbs.h, method="snewton", control=list(trace=0)))
if ((class(solx3) != "try-error") && (solx3$convergence == 0)) {
  proptimr(solx3)
  print(eigen(solx3$hessian)$values)
} else cat("solx3 failed!\n")
dirx3 <- try(snewton(x3, hobbs.f, hobbs.g, hobbs.h, control=list(trace=0)))
if ((class(dirx3) != "try-error") && (dirx3$convcode == 0)) {
  proptimr(dirx3)
  print(eigen(dirx3$Hess)$values)
} else cat("dirx3 failed!\n")


