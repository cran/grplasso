\name{plot.grplasso}
\alias{plot.grplasso}
\title{Plots the solution path of a grplasso object}
\description{Plots the solution path of a \code{grplasso} object.}
\usage{
\method{plot}{grplasso}(x, type = "coefficients", col = NULL, ...)
}
\arguments{
  \item{x}{a \code{grplasso} object}
  \item{type}{type = "coefficients" plots coefficients with respect to
    the input variables, even if \code{standardize = TRUE} is used in
    \code{grplasso}.}
  \item{col}{a vector indicating the color of the different group
    paths. The length should equal the number of \emph{groups}. The same 
    ordering as in the vector \code{index} is used with the exception
    that the unpenalized coefficients are grouped at the beginning of
    the vector.}
  \item{...}{other parameters to be passed to the plotting functions.}
}
\examples{
data(splice)

contr <- list(Pos.1 = "contr.sum", Pos.2 = "contr.sum")
lambda <- lambdamax(y ~ Pos.1 * Pos.2, data = splice, model = LogReg(),
                    contrasts = contr, standardize = TRUE) * 0.8^(0:8)

fit <- grplasso(y ~ Pos.1 * Pos.2, data = splice, model = LogReg(),
                lambda = lambda, contrasts = contr, standardize = TRUE,
                control = grpl.control(trace = 0, inner.loops = 0,
                update.every = 1, update.hess = "lambda"))
plot(fit, log = "x")
}
\keyword{hplot}
