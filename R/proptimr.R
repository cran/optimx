proptimr <- function(opres){
  # to provide compact output of an optimr result
  obname<-deparse(substitute(opres))
  cat("Result ",obname,"  proposes optimum function value =", opres$value," at parameters\n")
  print(opres$par )
  nfn <- opres$counts[[1]]
  ngr <- opres$counts[[2]]
  cat("After ",nfn," fn evals, and ",ngr," gr evals\n")
  cat("Termination code is ", opres$convergence,":", opres$message)
  if (! is.null(opres$grad)) {
    cat("Gradient:")
    print(opres$grad)
  }
  if (! is.null(opres$Hess)) {
    cat("Hessian:\n")
    print(opres$Hess)
  }
  cat("\n---------------------------------------\n")
}
# -------------- end proptimr ----------------- #
#################################################################
