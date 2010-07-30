% File src/library/stats/man/getresult.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{get.result}
\alias{get.result}
\encoding{UTF-8}
\title{Result extraction tool for answers from optimx}
\concept{minimization}
\concept{maximization}
\description{
  The optimx general-purpose optimization wrapper function that calls other
  R tools for optimization, including the existing optim() function creates
  a data frame of results. \code{get.results} extracts elements from this 
  data frame.
}
\usage{
get.result(optimx.obj, method, attribute) 
}
\arguments{
 \item{optimx.obj}{object returned by `optimx'}
 \item{method}{Name of the method for which all results are desired (e.g. `spg')}
 \item{attribute}{type of result desired for all methods (e.g. `fvalues').
      There is partial matching on inputs, so `fval' will recover `fvalues'.}
}
\details{
  None
}
\value{
  \code{get.result} returns the items requested or a warning string.
        If both method and attribute are provided, the function returns 
        a warning message as a character string
}
\keyword{nonlinear}
\keyword{optimize}

