\name{checksolver}
\alias{checksolver}
\encoding{UTF-8}
\title{Test if requested solver is present}
\concept{minimization}
\concept{maximization}
\description{
  Test if requested solver is present.
}
\usage{
   checksolver(method, allmeth, allpkg)
}
\arguments{
 \item{method}{Character string giving the name of the solver requested.}
 \item{allmeth}{Character vector giving the names of the methods optimr can use.}
 \item{allpkg}{Character vector giving the names of the packages where
     solvers are found.}
}
\value{
  \code{checksolver} tests if requested function minimization solver is present.
}
\examples{
   allmeth <- c("Rvmmin", "nlminb","ipopttest")
   allpkg <- c("Rvmmin", "stats","ipoptr")
   
   print(checksolver("nlminb", allmeth, allpkg))
   # If Rvmmin NOT available, get msg that PACKAGE not available.
   print(checksolver("Rvmmin", allmeth, allpkg))
   # Get message that SOLVER not found
   print(checksolver("notasolver", allmeth, allpkg))

}
\keyword{nonlinear}
\keyword{optimize}

