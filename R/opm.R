opm <- function(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=c("Nelder-Mead","BFGS"), hessian=FALSE,
            control=list(),
             ...) {

  npar <- length(par)
  pstring<-names(par)
  npar <- length(par)
  ctrl <- ctrldefault(npar)
  ncontrol <- names(control)
  nctrl <- names(ctrl)
  for (onename in ncontrol) {
     if (onename %in% nctrl) {
       if (! is.null(control[onename]) || ! is.na(control[onename]) )
       ctrl[onename]<-control[onename]
     }
  }
  control <- ctrl # note the copy back! control now has a FULL set of values
  ## 180706: Should we try to streamline?
  if(control$trace > 0) cat("opm: wrapper to call optimr to run multiple optimizers\n")

  fnscale <- 1 # default to ensure defined
  if (is.null(control$fnscale)) {
     if (! is.null(control$maximize) && control$maximize ) {fnscale <- -1}
  else if (! is.null(control$maximize)) {
          if ( (control$fnscale < 0) && control$maximize) {fnscale <- -1} # this is OK
          else stop("control$fnscale and control$maximize conflict")
       } # end ifelse
  } # end else
  control$fnscale <- fnscale # to ensure set again

  allmeth <- control$allmeth
  # 160628: uobyqa removed as it fails hobbs from 1,1,1 unscaled

  bdmeth <- control$bdmeth

  maskmeth <- control$maskmeth
  # Masks: As at 2016-6-28 do NOT provide for masks in package optimr


  bmtst <- bmchk(par, lower=lower, upper=upper)
  control$have.bounds <- bmtst$bounds # and set a control value
  bdmsk <- bmtst$bdmsk # Only need the masks bit from here on
  # These are set free (1) or set -1 for upper bounds, -3 for lower bounds
  # At this stage should NOT have masks (Or could they be added if upper=lower by bmchk
  control$have.masks <- any(bdmsk == 0)

  if (length(method) == 1 && method == "ALL") control$all.methods <- TRUE
  if (control$all.methods) {
       if (control$have.masks) { method <- maskmeth }
       else { if (control$have.bounds) { method <- bdmeth }
              else { method <- allmeth }
       }
  }
  nmeth <- length(method)

  if (is.null(pstring)) {
      for (j in 1:npar) {  pstring[[j]]<- paste("p",j,sep='')}
  } 
   cnames <- c(pstring, "value", "fevals", "gevals", "convergence", "kkt1", "kkt2", "xtime")
   ans.ret <- matrix(NA, nrow=nmeth, ncol=npar+7)
  if (control$trace > 2) {
      print(ans.ret)
      tmp <- readline("continue after printing ans.ret initial")
  }
  ans.ret <- data.frame(ans.ret)
  colnames(ans.ret)<-cnames
  row.names(ans.ret)<-method
  ans.details <- list()
  if (control$trace > 2) {
     cat("width of ans.ret =", npar+7,"\n")
     print(dim(ans.ret))
  }
  for (i in 1:nmeth) {
    meth <- method[i] # extract the method name
    if (control$trace > 0) cat("Method: ",meth,"\n")
    # Note: not using try() here
    if (is.character(gr) && (control$trace>0)) 
           cat("Using numerical gradient ",gr," for method ", meth,"\n")
    time <- system.time(ans <- optimr(par, fn, gr, hess=hess, method=meth, lower=lower, upper=upper, 
           hessian=hessian, control=control, ...))[1]
    if (control$trace > 2) print(ans)
    # add to list

## --------------------------------------------
## Post-processing -- Kuhn Karush Tucker conditions
#  Ref. pg 77, Gill, Murray and Wright (1981) Practical Optimization, Academic Press
      if (control$trace>0) { cat("Post processing for method ",meth,"\n") }
      if (exists("ans$message")) {
           amsg<-ans$message
           ans$message <- NULL # ?? do we need this
      } else { amsg <- "none" }
      ngatend <- NA
      nhatend <- NA
      hev <- NA
      ans$gevals <- ans$counts[2]
      ans$fevals <- ans$counts[1]
      ans$kkt1<-NA
      ans$kkt2<-NA
      kktres <- list(gmax=NA, evratio = NA, kkt1=NA, kkt2=NA, 
                     hev=rep(NA,npar), ngatend=NA, nhatend=NA)
      if ( control$save.failures || (ans$convergence < 1) ){
           # Save soln if converged or directed to save
          if ((control$trace > 0) && (ans$convergence==0)) cat("Successful convergence! \n") 
# Testing final soln. Use numDeriv for gradient & Hessian; compute Hessian eigenvalues
           if ((control$kkt || hessian) && (ans$convergence < 9900)) { # chg 160917 for no gradient
             wgr <- gr
             if (is.null(wgr)) wgr <- control$defgrapprox
             kktres <- kktchk(ans$par, fn, wgr, hess=NULL, upper=NULL, lower=NULL, 
                    maximize=control$maximize, control=control, ...) 
             ans$kkt1<-as.logical(kktres$kkt1)
             ans$kkt2<-as.logical(kktres$kkt2)
          }
# put together results
          ans$xtimes <- time
          # Do we want more information saved?
          if (control$trace > 1) { 
		cat("Save results from method ",meth,"\n") 
	  	print(ans)
	  }
	  if (control$trace > 2) { 
             cat("Assemble the answers\n") 
             cat("ans.ret now\n")
             print(ans.ret)
          }
          addvec <- c(ans$par, ans$value, ans$fevals, ans$gevals, 
                              ans$convergence, ans$kkt1, ans$kkt2, ans$xtimes)
	  if (control$trace > 2) { 
             cat("length addvec = ", length(addvec),"\n")
             print(addvec)
          }
          ans.ret[meth, ] <- addvec
      }  ## end post-processing of successful solution
      ans.details<-rbind(ans.details, list(method=meth, ngatend=kktres$ngatend, 
             nhatend=kktres$nhatend, hev=kktres$hev, message=amsg))
      # 1303234 try list() not c()
      row.names(ans.details)[[i]]<-meth
    } # End loop  ## end loop over method (index i)
    ansout <- NULL # default if no answers
    if (length(ans$par) > 0) { # cannot save if no answers
	ansout <- ans.ret # Don't seem to need drop=FALSE
        attr(ansout, "details")<-ans.details
        ansout[, "kkt1"] <- as.logical(ansout[, "kkt1"])
        ansout[, "kkt2"] <- as.logical(ansout[, "kkt2"])
    }
    ansout # return(ansout)
    answer <- structure(ansout, details = ans.details, maximize = control$maximize,
            npar = npar, class = c("opm", "data.frame"))

} ## end of opm

