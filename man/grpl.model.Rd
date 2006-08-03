\name{grpl.model}
\alias{grpl.model}
\alias{LogReg}
\alias{LinReg}
\alias{PoissReg}
\title{Group Lasso models}
\description{Generates models to be used for the Group Lasso algorithm.}
\usage{
grpl.model(invlink, nloglik, ngradient, nhessian, 
           name = "user-specified", comment = "user-specified")
LogReg()
LinReg()
PoissReg()
}
\arguments{
  \item{invlink}{a function with arguments \code{eta}
    implementing the inverse link function.}
  \item{nloglik}{a function with arguments \code{y}, \code{mu} and
    \code{weights} implementing the \emph{negative} log-likelihood function.}
  \item{ngradient}{a function with arguments \code{X}, \code{y},
    \code{mu} and \code{weights} 
    implementing the \emph{negative} gradient of the log-likelihood function.}
  \item{nhessian}{a function with arguments \code{X}, \code{mu} and
    \code{weights} implementing the \emph{negative} hessian of the
    log-likelihood function.}
  \item{name}{a character name}
  \item{comment}{a character comment}
}
\value{
  An object of class \code{grpl.model}.
}
\examples{
    LogReg()
}
\keyword{misc}
