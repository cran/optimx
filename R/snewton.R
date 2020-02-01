snewton<-function(par, fn, gr, hess, control=list(trace=0, maxit=500),...) {
## Safeguarded Newton minimizer with backtrack line search
##
##Input
##       - fn is the function we wish to minimize
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
nf <- ng <- nh <- niter <- 0 # counters

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
  msg <- "Snewton - no msg"
  xb <- par # best so far
  fbest <- fn(xb, ...)
  nf <- nf + 1 
  repeat { # MAIN LOOP
    if (trace > 0) cat(niter," ",nf," ",ng,"  fbest=",fbest,"\n")
    if (trace > 1) print(xb)
    niter <- niter + 1
    grd<-gr(xb,...) # compute gradient
    ng <- ng + 1
    halt <- FALSE # default is keep going
    # tests on too many counts??
    if (niter > control$maxit) {
      if (trace > 0) cat("Too many (",niter,") iterations\n")
      halt <- TRUE
      convcode <- 1
      break
    }
    if (nf > control$maxfeval){
      msg <- paste("Too many (",nf," function evaluations")
      if (trace > 0) cat(msg,"\n")
      halt <- TRUE
      convcode <- 91 # ?? value
      break
    }
    gmax <- max(abs(grd))
    if (trace > 1) cat("current gradient norm =",gmax,"\n")
    if (gmax <= control$epstol) {
      msg <- paste("Small gradient norm")
      if (trace > 0) cat(msg,"\n")
      halt <- TRUE
      convcode <- 0 # OK
      break
    }
    # Note if we get here, 
#    if (trace > 0) {cat("Iteration ",niter,":")}
    H<-hess(xb,...)
    nh <- nh + 1
    d<-try(solve(H, -grd))
    if (inherits(d, "try-error")) {
          stop("Failure of default solve of Newton equations")
    }
    if (trace > 2) {
         cat("Search vector:")
         print(d)
    }
    gprj <- as.numeric(crossprod(d, grd))
    if (trace > 0) cat("Gradient projection = ",gprj,"\n")
# ?? Do we want a test on size of gprj
    st <- control$defstep
    if (gprj > 0) st <- -st # added 180330 to allow downhill
    xnew <- xb + st*d # new point
    if (all((control$offset+xnew) == (control$offset+xb))) {
        convcode <- 92 # no progress
        msg <- "No progress before linesearch!"
        if (trace > 0) cat(msg,"\n")
        break # finished        
    }
    fval <- try(fn(xnew, ...))
    nf <- nf + 1
    if (inherits(fval, "try-error")) stop("snewton: function evaluation error")
    if (trace > 1) {
       cat("f(xnew)=",fval," at ")
       print(xnew)
    }
    if (fval > control$bigval) {
       msg <- "snewton: New function value infinite"
       if (trace > 1) cat(msg,"\n")
#       fval <- control$bigval
       convcode <- 9999
       break
    }
    while ((fval > fbest + control$acctol*abs(st)*gprj) # 180330 Note abs(st)
           && (all((control$offset+xnew) != (control$offset+xb)))) { 
        # continue until satisfied
        st <- st * control$stepdec
        if (trace > 1) cat("Stepsize now =",st,"\n")
        xnew <- xb + st*d # new point
        fval <- fn(xnew, ...)    
        nf <- nf + 1
        if (trace > 1) cat("* f(xnew)=",fval,"\n")
        if (fval > control$bigval) {
           stop("snewton: Function value infinite in backtrack")
        }
    } # end while
    if (all((control$offset+xnew) == (control$offset+xb))) {
        convcode <- 93 # no progress in linesearch
        msg <- "No progress in linesearch!"
        if (trace > 0) cat(msg,"\n")
        break
    }
    if (trace > 1) cat("end major loop\n")  
    xb <- xnew
    fbest <- fval
    if (control$watch) { tmp <- readline("end iteration") }
  } # end repeat
  out<-NULL
  out$par<-xb
  out$value<-fbest
  out$grad<-grd
  out$Hess<-H
  out$counts <- list(niter=niter,  nfn=nf, ngr=ng, nhess=nh)
  out$convcode <- convcode
  out$message <- msg
  out
}
