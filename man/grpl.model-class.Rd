\name{grpl.model-class}
\docType{class}
\alias{grpl.model-class}
\alias{show,grpl.model-method}
\title{Class "grpl.model": Group Lasso Models}
\description{Objects of class "grpl.model" define link function, negative
log-likelihood and corresponding gradient and Hessian for the model to
be used in a Group Lasso problem.} 
\section{Objects from the Class}{
Objects can be created by calls of the form \code{grpl.model(...)}}
\section{Slots}{
  \describe{
    \item{\code{invlink}}{a function with arguments \code{eta}
      implementing the inverse link function.}
    \item{\code{link}}{a function with arguments \code{mu}
      implementing the link function.}
    \item{\code{nloglik}}{a function with arguments \code{y}, \code{mu} and
    \code{weights} implementing the \emph{negative}
    log-likelihood function.} 
    \item{\code{ngradient}}{a function with arguments \code{x}, \code{y},
    \code{mu} and \code{weights} 
    implementing the \emph{negative} gradient of the log-likelihood function.} 
    \item{\code{nhessian}}{a function with arguments \code{x}, \code{mu} and
    \code{weights} implementing the \emph{negative} hessian of the
    log-likelihood function.}
    \item{\code{check}}{a function with argument \code{y} to check
      whether the response has the correct format.}
    \item{\code{name}}{a character name}
    \item{\code{comment}}{a character comment}
  }
}
\section{Methods}{
\describe{
\item{show}{object}
}
}
\examples{
  LogReg()
}
\keyword{classes}
