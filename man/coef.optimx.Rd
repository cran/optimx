\name{coef.optimx}
\alias{coef.optimx}
\alias{coef<-.optimx}
\alias{coef<-}
\encoding{UTF-8}
\title{Summarize optimx object}
\concept{minimization}
\concept{maximization}
\description{
  Summarize an \code{"optimx"} object.
}
\usage{
  \method{coef}{optimx}(object, ...) 
  \method{coef}{optimx}(x) <- value 
}
\arguments{
 \item{object}{Object returned by \code{optimx}.}
 \item{\dots}{Further arguments to be passed to the function. Currently not
   used.}
 \item{x}{An \code{optimx} object.}
 \item{value}{Set parameters equal to this value.}
}
\value{
  \code{coef.optimx} returns the best parameters found by each method that
  returned such parameters. The returned coefficients are in the form of a
  matrix with the rows named by the relevant methods and the columns named
  according to parameter names provided by the user in the vector of starting
  values, or else by "p1", "p2", ..., if names are not provided.
}
\examples{
ans <- optimx(fn = function(x) sum(x*x), par = 1:2)
coef(ans)

\dontrun{
proj <- function(x) x/sum(x)
f <- function(x) -prod(proj(x))
ans <- optimx(1:2, f)
ans
coef(ans) <- apply(coef(ans), 1, proj)
ans
}

}
\keyword{nonlinear}
\keyword{optimize}

