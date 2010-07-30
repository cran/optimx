% File src/library/stats/man/optimx.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{optimx}
\alias{optimx}
\encoding{UTF-8}
\title{General-purpose optimization}
\concept{minimization}
\concept{maximization}
\description{
  General-purpose optimization wrapper function that calls other
  R tools for optimization, including the existing optim() function.
  \code{optimx} also tries to unify the calling sequence to allow
  a number of tools to use the same front-end. These include 
  \code{spg} from the BB package, \code{ucminf}, \code{nlm}, and 
  \code{nlminb}. Note that 
  optim() itself allows Nelder--Mead, quasi-Newton and 
  conjugate-gradient algorithms as well as box-constrained optimization 
  via L-BFGS-B. Because SANN does not return a meaningful convergence code
  (conv), optimx() does not call the SANN method.
}
\usage{
optimx(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=c("Nelder-Mead","BFGS"), itnmax=NULL, hessian=NULL,
            control=list(),
             ...)

}
\arguments{
 \item{par}{A vector of initial values for the parameters for which optimal 
   values are to be found. }
 \item{fn}{A function to be minimized (or maximized), with first
   argument the vector of parameters over which minimization is to take
   place.  It should return a scalar result.}
 \item{gr}{A function to return (as a vector) the gradient for those methods that 
   can use this information.
   This includes the following methods:
		\code{"BFGS"}
		\code{"CG"}
		\code{"L-BFGS-B"}


   If 'gr' is \code{NULL}, a finite-difference approximation will be used.
   An open question concerns whether the SAME approximation code used for all methods, 
   or whether there are differences that could/should be examined? }

 \item{hess}{A function to return (as a symmetric matrix) the Hessian of the objective 
   function for those methods that can use this information.}
 \item{lower, upper}{Bounds on the variables for methods such as \code{"L-BFGS-B"} that can
   handle box (or bounds) constraints.}
 \item{method}{A list of the methods to be used. 
       Note that this is an important change from optim() that allows
       just one method to be specified. See \sQuote{Details}.}
 \item{itnmax}{If provided as a vector of the same length as the list of methods \code{method}, 
	gives the maximum number of iterations or function values for the corresponding 
	method. If a single number is provided, this will be used for all methods. Note that
	there may be control list elements with similar functions, but this should be the
	preferred approach when using \code{optimx}.}
 \item{hessian}{A logical control that if TRUE forces the computation of an approximation 
       to the Hessian at the final set of parameters. If FALSE, ensures it is not calculated.
       This setting is provided for compatibility with optim(). See how to set \code{kkt} in
       \sQuote{Details}.}
 \item{control}{A list of control parameters. See \sQuote{Details}.}
 \item{\dots}{Further arguments to be passed to \code{fn} and \code{gr}.}
}
\details{
  Note that arguments after \code{\dots} must be matched exactly.

  By default this function performs minimization, but it will maximize
  if \code{control$maximize} is TRUE. The original optim() function allows
  \code{control$fnscale} to be set negative to accomplish this. DO NOT
  use both methods. 

  Possible method codes at the time of writing are 'Nelder-Mead', 'BFGS',
  'CG', 'L-BFGS-B', 'nlm', 'nlminb', 'spg', and 'ucminf'.

  The default methods are an implementation of the Nelder and Mead
  (1965) and a Variable Metric method based on the ideas of Fletcher
  (1970) as modified by him in conversation with Nash (1979). Nelder-Mead
  uses only function values and is robust but relatively slow.  It will 
  work reasonably well for non-differentiable functions. The Variable
  Metric method, \code{"BFGS"} updates an approximation to the inverse
  Hessian using the BFGS update formulas, along with an acceptable point
  line search strategy. This method appears to work best with analytic
  gradients. As at June 2009, a box-constrained version is being tested.

  Method \code{"CG"} is a conjugate gradients method based on that by
  Fletcher and Reeves (1964) (but with the option of Polak--Ribiere or
  Beale--Sorenson updates). The particular implementation is now dated,
  and improved yet simpler codes are being implemented (as at June 2009),
  and furthermore a version with box constraints is being tested.
  Conjugate gradient methods will generally be more fragile than the 
  BFGS method, but as they do not store a matrix they may be successful 
  in much larger optimization problems.

  Method \code{"L-BFGS-B"} is that of Byrd \emph{et. al.} (1995) which
  allows \emph{box constraints}, that is each variable can be given a lower
  and/or upper bound. The initial value must satisfy the constraints.
  This uses a limited-memory modification of the BFGS quasi-Newton
  method. If non-trivial bounds are supplied, this method will be
  selected, with a warning.

  Nocedal and Wright (1999) is a comprehensive reference for the
  previous three methods.

  Function \code{fn} can return \code{NA} or \code{Inf} if the function
  cannot be evaluated at the supplied value, but the initial value must
  have a computable finite value of \code{fn}. However, some methods, of
  which \code{"L-BFGS-B"} is known to be a case, require that the values
  returned should always be finite.

  While \code{optim} can be used recursively, and for a single parameter
  as well as many, this may not be true for \code{optimx}. \code{optim}
  also accepts a zero-length \code{par}, and just evaluates the function 
  with that argument.

  Method \code{"nlm"} is from the package of the same name that implements
  ideas of Dennis and Schnabel (1983) and Schnabel et al. (1985). See nlm()
  for more details.

  Method \code{"nlminb"} is the package of the same name that uses the
  minimization tools of the PORT library.  The PORT documentation is at 
  <URL: http://netlib.bell-labs.com/cm/cs/cstr/153.pdf>. See nlminb()
  for details. (Though there is very little information about the methods.)

  Method \code{"spg"} is from package BB implementing a spectral projected 
  gradient method for large-scale optimization with simple constraints due
  R adaptation, with significant modifications, by Ravi Varadhan,
  Johns Hopkins University (March 25, 2008), from the original
  FORTRAN code of Birgin, Martinez, and Raydan (2001). 

  Method \code{"Rcgmin"} is from the package of that name. It implements a
  conjugate gradient algorithm with the Yuan/Dai update (ref??) and also 
  allows bounds constraints on the parameters. (Rcgmin also allows mask 
  constraints -- fixing individual parameters -- but there is no interface
  from \code{"optimx"}.) 

  Methods \code{"bobyqa"}, \code{"uobyqa"} and \code{"newuoa"} are from the 
  package \code{"minqa"} which implement optimization by quadratic approximation
  routines of the similar names due to M J D Powell (2009). See package minqa 
  for details. Note that \code{"uobyqa"} and \code{"newuoa"} are for 
  unconstrained minimization, while \code{"bobyqa"} is for box constrained
  problems.

  The \code{control} argument is a list that can supply any of the
  following components:
  \describe{
    \item{\code{trace}}{Non-negative integer. If positive,
      tracing information on the
      progress of the optimization is produced. Higher values may
      produce more tracing information: for method \code{"L-BFGS-B"}
      there are six levels of tracing. trace = 0 gives no output 
      (To understand exactly what these do see the source code: higher 
      levels give more detail.)}
    \item{\code{follow.on }}{ = TRUE or FALSE. If TRUE, and there are multiple 
      methods, then the last set of 
      parameters from one method is used as the starting set for the next. }
    \item{\code{save.failures}}{ = TRUE if we wish to keep "answers" from runs 
      where the method does not return conv==0. FALSE otherwise (default).}
    \item{\code{maximize}}{ = TRUE if we want to maximize rather than minimize 
      a function. (Default FALSE). Methods nlm, nlminb, ucminf cannot maximize a
      function, so the user must explicitly minimize and carry out the adjustment
      externally. However, there is a check to avoid
      usage of these codes when maximize is TRUE. See \code{fnscale} below for 
      the method used in \code{optim} that we deprecate.}
    \item{\code{all.methods}}{= TRUE if we want to use all available (and suitable)
      methods.}
    \item{\code{kkt}}{=FALSE if we do NOT want to test the Kuhn, Karush, Tucker
      optimality conditions. Default is TRUE. We return logical values kkt1 and 
      kkt2 TRUE if first and second order conditions are satisfied approximately.
      Note, however, that the tests are sensitive to scaling, and users may need
      to perform additional verification. }
    \item{\code{all.methods}}{= TRUE if we want to use all available (and suitable)
      methods.}
    \item{\code{sort.result}}{= TRUE if we want to have the results table sorted 
      into reverse order of final objective function.}
    \item{\code{kkttol}}{= value to use to check for small gradient and negative
      Hessian eigenvalues. Default = .Machine$double.eps^(1/3) }
    \item{\code{kkt2tol}}{= Tolerance for eigenvalue ratio in KKT test of positive 
      definite Hessian. Default same as for kkttol }
    \item{\code{starttests}}{= TRUE if we want to run tests of the function and 
      parameters: feasibility relative to bounds, analytic vs numerical gradient, 
      scaling tests, before we try optimization methods. Default is TRUE.}
    \item{\code{dowarn}}{= TRUE if we want warnings generated by optimx. Default is 
      TRUE.}
    \item{\code{badval}}{= The value to set for the function value when try(fn()) fails.
      Default is (0.5)*.Machine$double.xmax }
  }

  The following \code{control} elements apply only to some of the methods. The list
  may be incomplete. See individual packages for details. 

  \describe{

    \item{\code{fnscale}}{An overall scaling to be applied to the value
      of \code{fn} and \code{gr} during optimization. If negative,
      turns the problem into a maximization problem. Optimization is
      performed on \code{fn(par)/fnscale}. For \code{optim}.}
    \item{\code{parscale}}{A vector of scaling values for the parameters.
	Optimization is performed on \code{par/parscale} and these should be
	comparable in the sense that a unit change in any element produces
	about a unit change in the scaled value.For \code{optim}.}
    \item{\code{ndeps}}{A vector of step sizes for the finite-difference
      approximation to the gradient, on \code{par/parscale}
      scale. Defaults to \code{1e-3}. For \code{optim}.}
    \item{\code{maxit}}{The maximum number of iterations. Defaults to
      \code{100} for the derivative-based methods, and
      \code{500} for \code{"Nelder-Mead"}.}
    \item{\code{abstol}}{The absolute convergence tolerance. Only
      useful for non-negative functions, as a tolerance for reaching zero.}
    \item{\code{reltol}}{Relative convergence tolerance.  The algorithm
      stops if it is unable to reduce the value by a factor of
      \code{reltol * (abs(val) + reltol)} at a step.  Defaults to
      \code{sqrt(.Machine$double.eps)}, typically about \code{1e-8}. For \code{optim}.}
    \item{\code{alpha}, \code{beta}, \code{gamma}}{Scaling parameters
      for the \code{"Nelder-Mead"} method. \code{alpha} is the reflection
      factor (default 1.0), \code{beta} the contraction factor (0.5) and
      \code{gamma} the expansion factor (2.0).}
    \item{\code{REPORT}}{The frequency of reports for the \code{"BFGS"} and
      \code{"L-BFGS-B"} methods if \code{control$trace}
      is positive. Defaults to every 10 iterations for \code{"BFGS"} and
      \code{"L-BFGS-B"}.}
    \item{\code{type}}{for the conjugate-gradients method. Takes value
      \code{1} for the Fletcher--Reeves update, \code{2} for
      Polak--Ribiere and \code{3} for Beale--Sorenson.}
    \item{\code{lmm}}{is an integer giving the number of BFGS updates
      retained in the \code{"L-BFGS-B"} method, It defaults to \code{5}.}
    \item{\code{factr}}{controls the convergence of the \code{"L-BFGS-B"}
      method. Convergence occurs when the reduction in the objective is
      within this factor of the machine tolerance. Default is \code{1e7},
      that is a tolerance of about \code{1e-8}.}
    \item{\code{pgtol}}{helps control the convergence of the \code{"L-BFGS-B"}
      method. It is a tolerance on the projected gradient in the current
      search direction. This defaults to zero, when the check is
      suppressed.}
  }

  Any names given to \code{par} will be copied to the vectors passed to
  \code{fn} and \code{gr}.  Note that no other attributes of \code{par}
  are copied over. (We have not verified this as at 2009-07-29.)

}
\value{
  A list with components:
  \item{par}{The best set of parameters found.}
  \item{value}{The value of \code{fn} corresponding to \code{par}.}
  \item{counts}{A two-element integer vector giving the number of calls
    to \code{fn} and \code{gr} respectively. This excludes those calls needed
    to compute the Hessian, if requested, and any calls to \code{fn} to
    compute a finite-difference approximation to the gradient.}
  \item{convergence}{An integer code. \code{0} indicates successful
    convergence. Error codes are
    \describe{
      \item{\code{1}}{indicates that the iteration limit \code{maxit}
      had been reached.}
      \item{\code{10}}{indicates degeneracy of the Nelder--Mead simplex.}
      \item{\code{51}}{indicates a warning from the \code{"L-BFGS-B"}
      method; see component \code{message} for further details.}
      \item{\code{52}}{indicates an error from the \code{"L-BFGS-B"}
      method; see component \code{message} for further details.}
    }
  }
  \item{message}{A character string giving any additional information
    returned by the optimizer, or \code{NULL}.}
  \item{hessian}{Only if argument \code{hessian} is true. A symmetric
    matrix giving an estimate of the Hessian at the solution found. Note
    that this is the Hessian of the unconstrained problem even if the
    box constraints are active.}

More detail is provided in the attribute "details" to the returned answer object. 
If the returned object from optimx() is \code{ans}, this is accessed 
via the construct
    \code{attr(ans, "details")}
This object is a list of lists, one list for each method that has been successful 
or that has been forced by save.failures=TRUE. To get the details list for the 
third method, we use
    \code{attr(ans, "details")[[3]]]]}

Each such list has possible elements:
  \item{par}{ - the final parameters of the function}
  \item{value}{ - the final value of the objective function}
  \item{convergence}{ - a code indicating whether and/or how the method terminated}
  \item{message}{ - a descriptive message, if one is provided}
  \item{fevals}{ - the number of times the objective function was evaluated during the run;
   this item has an attribute "names" that can take on the value "function"}
  \item{gevals}{ - the number of times the gradient of the objective function was evaluated 
   during the run; this item has an attribute "names" that can take on the value "gradient"}
  \item{kkt1}{ - an indicator if the first order Kuhn, Karush, Tucker condition for a 
  local minimum is satisfied. Note that we must use a tolerance to make this test, and the
  indicator is NOT definitive.}
  \item{kkt2}{  - an indicator if the second order Kuhn, Karush, Tucker condition for a 
  local minimum is satisfied. Note that we must use a tolerance to make this test, and the
  indicator is NOT definitive.}
  \item{ngatend}{ - the gradient vector computed or estimated at termination}
  \item{nhatend}{ - the Hessian matrix computed or estimated at termination}
  \item{evnhatend}{ - a vector of the eigenvalues of the Hessian at termination}
  \item{systime}{ - the system time required for execution of the method in question (named
   in "method" below). This element has the attribute "names".}
  \item{method}{ - The name of the method used. This has the attribute "CPU times (s)" which 
  itself has an attribute "names" that corresponds to the "method".}

}
\note{
  \code{optimx} will work with one-dimensional \code{par}s, but the
  default method does not work well (and will warn).  Use
  \code{\link{optimize}} instead.

  There are a series of demos available. These were set up as tests, but take quite
  a long time to execute. Once the package is loaded (via \code{require(optimx)} or
  \code{library(optimx)}, you may see available demos via 

  demo(package="optimx")

  The demo 'brown_test' may be run with the command
  demo(brown_test, package="optimx")
}
\source{
  The code for methods \code{"Nelder-Mead"}, \code{"BFGS"} and
  \code{"CG"} was based originally on Pascal code in Nash (1990) that was
  translated by \code{p2c} and then hand-optimized.  Dr Nash has agreed
  that the code can be made freely available.

  The code for method \code{"L-BFGS-B"} is based on Fortran code by Zhu,
  Byrd, Lu-Chen and Nocedal obtained from Netlib (file
  \file{opt/lbfgs\_bcm.shar}: another version is in \file{toms/778}).

  See documentation for \code{"nlm"} and \code{"nlminb"} for those
  methods, package \code{"ucminf"} for the funciton of the same
  name, and package BB for method \code{"spg"}.
}
\references{
  Belisle, C. J. P. (1992) Convergence theorems for a class of simulated
  annealing algorithms on \eqn{R^d}{Rd}. \emph{J Applied Probability},
  \bold{29}, 885--895.

   Birgin EG, Martinez JM, and Raydan M (2001): SPG: software for
     convex-constrained optimization, \emph{ACM Transactions on Mathematical
     Software}. ??incomplete reference??


  Byrd, R. H., Lu, P., Nocedal, J. and Zhu, C.  (1995) A limited
  memory algorithm for bound constrained optimization.
  \emph{SIAM J. Scientific Computing}, \bold{16}, 1190--1208.

  Dennis, J. E. and Schnabel, R. B. (1983) _Numerical Methods for
     Unconstrained Optimization and Nonlinear Equations._
     Prentice-Hall, Englewood Cliffs, NJ.

  Fletcher, R. and Reeves, C. M. (1964) Function minimization by
  conjugate gradients. \emph{Computer Journal} \bold{7}, 148--154.

  Fletcher, R (1970) A new approach to variable metric algorithms,
  \emph{Computer Journal}, \bold{13}, 317-322.

  Nash, J. C. (1979, 1990) \emph{Compact Numerical Methods for
    Computers. Linear Algebra and Function Minimisation.} Adam Hilger.

  Nelder, J. A. and Mead, R. (1965) A simplex algorithm for function
  minimization. \emph{Computer Journal} \bold{7}, 308--313.

  Nocedal, J. and Wright, S. J. (1999) \emph{Numerical Optimization}.
  Springer.

  Schnabel, R. B., Koontz, J. E. and Weiss, B. E. (1985) A modular
     system of algorithms for unconstrained minimization. _ACM Trans.
     Math. Software_, *11*, 419-440.

}

\seealso{
  \code{\link{nlm}}, \code{\link{nlminb}}.

  \code{\link{optimize}} for one-dimensional minimization and
  \code{\link{constrOptim}} for constrained optimization.
}

\examples{
require(graphics)

fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
    x1 <- x[1]
    x2 <- x[2]
    c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
       200 *      (x2 - x1 * x1))
}
ans1<-optimx(c(-1.2,1), fr)
print(ans1)
print(attr(ans1,"details"))
cat("\n\n")
ans2<-optimx(c(-1.2,1), fr, grr, method = "BFGS")
print(ans2)
## The next line will fail if executed because 'hessian = TRUE' no longer allowed
# ans3<-optimx(c(-1.2,1), fr, NULL, method = "BFGS", hessian = TRUE)
cat("\n\n")
ans4<-optimx(c(-1.2,1), fr, grr, method = "CG",control=list(trace=TRUE))
print(ans4)
cat("\n\n")
ans5<-optimx(c(-1.2,1), fr, grr, method = "CG", control=list(type=2))
print(ans5)
cat("\n\n")
ans6<-optimx(c(-1.2,1), fr, grr, method = "L-BFGS-B")
print(ans6)
cat("\n\n")

flb <- function(x)
    { p <- length(x); sum(c(1, rep(4, p-1)) * (x - c(1, x[-p])^2)^2) }
## 25-dimensional box constrained
optimx(rep(3, 25), flb, NULL, method = "L-BFGS-B",
      lower=rep(2, 25), upper=rep(4, 25)) # par[24] is *not* at boundary

## "wild" function , global minimum at about -15.81515
fw <- function (x)
    10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80
plot(fw, -50, 50, n=1000, main = "optim() minimising 'wild function'")

## Suppressed for optimx() ans7 <- optimx(50, fw, method="SANN",
##             control=list(maxit=20000, temp=20, parscale=20))
## ans7
## Now improve locally {typically only by a small bit}:
## newpar<-unlist(ans7$par) # NOTE: you need to unlist the parameters as optimx() has multiple outputs
##(r2 <- optimx(newpar, fw, method="BFGS"))
##points(r2$par, r2$value, pch = 8, col = "red", cex = 2)

## Show multiple outputs of optimx using all.methods
# genrose function code
genrose.f<- function(x, gs=NULL){ # objective function
## One generalization of the Rosenbrock banana valley function (n parameters)
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	fval<-1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
        return(fval)
}

genrose.g <- function(x, gs=NULL){
# vectorized gradient for genrose.f
# Ravi Varadhan 2009-04-03
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	gg <- as.vector(rep(0, n))
	tn <- 2:n
	tn1 <- tn - 1
	z1 <- x[tn] - x[tn1]^2
	z2 <- 1 - x[tn]
	gg[tn] <- 2 * (gs * z1 - z2)
	gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1
	return(gg)
}

genrose.h <- function(x, gs=NULL) { ## compute Hessian
   if(is.null(gs)) { gs=100.0 }
	n <- length(x)
	hh<-matrix(rep(0, n*n),n,n)
	for (i in 2:n) {
		z1<-x[i]-x[i-1]*x[i-1]
		z2<-1.0-x[i]
                hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
                hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
                hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
                hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
	}
        return(hh)
}

startx<-4*seq(1:10)/3.
ans8<-optimx(startx,fn=genrose.f,gr=genrose.g, hess=genrose.h, control=list(all.methods=TRUE, save.failures=TRUE), gs=10)
print(ans8)

get.result(ans8, attribute="grs")
get.result(ans8, method="spg")


startx<-4*seq(1:10)/3.
cat("Polyalgorithm with 200 steps NM followed by up to 75 of ucminf\n")
ans9<-optimx(startx,fn=genrose.f,gr=genrose.g, hess=genrose.h, method=c("Nelder-Mead","ucminf"),
             itnmax=c(200,75), control=list(follow.on=TRUE, save.failures=TRUE,trace=TRUE), gs=10)
print(ans9)


}
\keyword{nonlinear}
\keyword{optimize}
