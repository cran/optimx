% File src/library/stats/man/optansout.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{optansout}
\alias{optansout}
\encoding{UTF-8}
\title{Output tool for answers from optimx}
\concept{minimization}
\concept{maximization}
\description{
  The optimx general-purpose optimization wrapper function that calls other
  R tools for optimization, including the existing optim() function creates
  a data frame of results. \code{optansout} displays this data frame and 
  optionally saves it in text form to a file. The parameter vectors, which
  may be very long, are summarized by the first 5 parameters if there are 
  more than 5.

  This routine has also been added to other optimization packages for
  convenience in output of optimization results.
}
\usage{
optansout(ansdf, filename)
}
\arguments{
 \item{ansdf}{A dataframe of results from an optimx run}
 \item{filename}{A properly qualified filename where text data is to be saved. If NULL
   output is only to the terminal. }
}
\details{
  None
}
\value{
  \code{optansout} always returns TRUE.
}
\keyword{nonlinear}
\keyword{optimize}
