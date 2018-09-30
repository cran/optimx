snewtonm<-function(par,fn,gr,hess,control=list(trace=0, maxit=500),...) {
  ## Safeguarded Newton minimizer 
  ##
  ##Input
  ##       - fn is the function we wish to minimize
  ##       - gr is the gradient function
  ##       - hess is the hessian function
  ##       - ... are exogenous data
  ##
  ##?? fixup documentation here??
  ##       - par is the initial value
  ##       - ... is data used in the function fn
  ##Output (list) -- need to match optim() output!! ???
  ##       - xs is the value at the minimum
  ##       - fv is the fn evaluated at xs
  ##       - grd is the gradient
  ##       - Hess is the Hessian
  ##       - niter is the number of interations needed (gradient and Hessian evals).
  ##       - add fevals??, other reports

npar <- length(par)


ctrldefault <- list(
  trace = 0,
  maxit = 500,
  maxfeval = npar*500,
  acctol = 0.0001,
  epstol = .Machine$double.eps,
  stepdec = 0.2, 
  stepmax = 5,
  stepmin = 0,
  offset = 100.0,
  defstep=1,
  bigval = .Machine$double.xmax*0.01,
  watch = FALSE
)  

ncontrol <- names(control)
nctrld <- names(ctrldefault)
for (onename in nctrld) {
  if (! (onename %in% ncontrol)) {
    control[onename]<-ctrldefault[onename]
  }
}
trace <- control$trace # convenience

#  cat(" Start snewtonm ")
  convcode <- 0
  nfn <- 0
  ngr <- 0
  nhess <- 0
  eps0<-.Machine$double.eps
  eps <- 10*eps0
  fval <- fn(par, ...)
  nfn <- nfn + 1
  if (control$trace > 0) {
   cat("  f0=",fval,"  at  ")
   print(par)
  }
  fbest <- abs(fval)*1.1 + 100.
  itn <- 1
  lambdamin<-(eps0^(1/4)) ## ?? do better
  laminc <- 10.
  lamdec <- 0.15
  lambda<- lambdamin/lamdec## ?? do better
  while ( (itn <= control$maxit) && (fval < fbest) ) { ## main loop
    fbest <- fval
    lambda <- lambda * lamdec
    if (control$trace>0) {cat(itn," ",nfn," ",ngr," fbest=",fbest,"\n")}
    grd<-gr(par,...)
    ngr <- ngr + 1
    if ( max(abs(grd)) < eps ) {
        cat("Small gradient\n")
        break
    }
    H<-hess(par,...)
    nhess <- nhess + 1
    while (fval >= fbest){
       itn <- itn+1
       Haug<-H + diag(npar)*lambda # To avoid singularity
       stp<-solve(Haug, -grd)
       xn <- par + stp # try unit step
       if (identical(par,xn)) {
#       if (max(abs(stp)) <= 1000*eps*max(abs(par))){
             break
       }
       #    if (control$trace) {cat(" step =", gvl,"  fval=", attr(gvl,"Fval"),"\n")}
       fval <- fn(xn, ...)
       nfn <- nfn + 1
       if (control$trace) {cat(" lambda =", lambda,"  fval=", fval,"\n")}
       if (fval >= fbest) {lambda <- max(lambdamin,lambda)*laminc} # increase lambda
    }
#    if (control$trace>0) {cat(" Success at lambda =", lambda,"  fval=", fval,"\n")}
     par<-xn
     if (control$watch>0) { tmp <- readline("end iteration") }

  }
  if (itn >= control$maxit) {
     msg <- "snewtonm: Too many iterations!"
     if(control$trace > 0) cat(msg,"\n")
     convcode <- 1
  } else { msg <- "snewtonm: Normal exit" }
  out<-NULL
  out$par<-xn
  out$value<-fn(xn,...)
  out$grad<-grd
  out$Hess<-H
  out$counts <- list(niter=itn, nfn=nfn, ngr=ngr, nhess=nhess)
  out$convcode <- convcode
  out$message <- msg
  out 
}
