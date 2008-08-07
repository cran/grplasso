\name{grpl.control}
\alias{grpl.control}
\title{Options for the Group Lasso Algorithm}
\description{
  Definition of options such as bounds on the Hessian,
  convergence criteria and output management for the Group Lasso algorithm.
}
\usage{
grpl.control(save.x = FALSE, save.y = TRUE,
             update.hess = c("lambda", "always"), update.every = 3,
             inner.loops = 10, line.search = TRUE, max.iter = 500,
             tol = 5 * 10^-8, lower = 10^-2, upper = Inf, beta = 0.5,
             sigma = 0.1, trace = 1)
}
\arguments{
  \item{save.x}{a logical indicating whether the design matrix should be saved.}
  \item{save.y}{a logical indicating whether the response should be saved.}
  \item{update.hess}{should the hessian be updated in each
    iteration ("always")? update.hess = "lambda" will update
    the Hessian once for each component of the penalty
    parameter "lambda" based on the parameter estimates
    corresponding to the previous value of the penalty
    parameter.}
  \item{update.every}{Only used if update.hess = "lambda". E.g. set to 3
    if you want to update the Hessian only every third grid point.}
  \item{inner.loops}{How many loops should be done (at maximum) when solving
    only the active set (without considering the remaining
    predictors). Useful if the number of predictors is large. Set to 0
    if no inner loops should be performed.}
  \item{line.search}{Should line searches be performed?}
  \item{max.iter}{Maximal number of loops through all groups}
  \item{tol}{convergence tolerance; the smaller the more precise, see
    details below.}
  \item{lower}{lower bound for the diagonal approximation of the
    corresponding block submatrix of the Hessian of the negative
    log-likelihood function.}
  \item{upper}{upper bound for the diagonal approximation of the
    corresponding block submatrix of the Hessian of the negative
    log-likelihood function.}
  \item{beta}{scaling factor \eqn{\beta < 1} of the Armijo line search.}
  \item{sigma}{\eqn{0 < \sigma < 1} used in the Armijo line search.}
  \item{trace}{integer. \code{0} omits any output, \code{1} prints the
    current lambda value, \code{2} prints the improvement in the
    objective function after each sweep through all the parameter groups
    and additional information.} 
}
\details{For the convergence criteria see chapter 8.2.3.2 of Gill et
  al. (1981).}
\references{Philip E. Gill, Walter Murray and Margaret H. Wright (1981)
  \emph{Practical Optimization}, Academic Press.
  
  Dimitri P. Bertsekas (2003) \emph{Nonlinear Programming}, Athena Scientific.}
\value{
  An object of class \code{grpl.control}.
}
\keyword{misc}
