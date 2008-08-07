
\name{grpl.control-class}
\docType{class}
\alias{grpl.control-class}
\title{Class "grpl.control": Options for the Group Lasso Algorithm}
\description{
  Objects of class "grpl.control" define options such as bounds on the Hessian,
  convergence criteria and output management for the Group Lasso algorithm.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{grpl.control(...)}
}
\section{Slots}{
  \describe{
  \item{\code{save.x}}{a logical indicating whether the design matrix
    should be saved.} 
  \item{\code{save.y}}{a logical indicating whether the response should
    be saved.} 
  \item{\code{update.hess}}{should the hessian be updated in each
    iteration ("always")? update.hess = "lambda" will update
    the Hessian once for each component of the penalty
    parameter "lambda" based on the parameter estimates
    corresponding to the previous value of the penalty
    parameter.}
  \item{\code{update.every}}{Only used if update.hess = "lambda". E.g. set to 3
    if you want to update the Hessian only every third grid point.}
  \item{\code{inner.loops}}{How many loops should be done (at maximum)
    when solving only the active set (without considering the remaining
    predictors). Useful if the number of predictors is large. Set to 0
    if no inner loops should be performed.}
  \item{\code{line.search}}{Should line searches be performed?}
  \item{\code{max.iter}}{Maximal number of loops through all groups}
  \item{\code{tol}}{convergence tolerance; the smaller the more precise.}
  \item{\code{lower}}{lower bound for the diagonal approximation of the
    corresponding block submatrix of the Hessian of the negative
    log-likelihood function.}
  \item{\code{upper}}{upper bound for the diagonal approximation of the
    corresponding block submatrix of the Hessian of the negative
    log-likelihood function.}
  \item{\code{beta}}{scaling factor \eqn{\beta < 1} of the Armijo line search.}
  \item{\code{sigma}}{\eqn{0 < \sigma < 1} used in the Armijo line search.}
  \item{\code{trace}}{integer. \code{1} prints the current lambda value,
    \code{2} prints the improvement in the objective function after each
    sweep through all the parameter groups and additional information.}
  }
}
\details{For the convergence criteria see chapter 8.2.3.2 of Gill et
  al. (1981).}
\references{Philip E. Gill, Walter Murray and Margaret H. Wright (1981)
  \emph{Practical Optimization}, Academic Press.
  
  Dimitri P. Bertsekas (2003) \emph{Nonlinear Programming}, Athena Scientific.}
\keyword{classes}
