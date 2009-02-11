\name{predict.grplasso}
\alias{predict.grplasso}

\title{Predict method for grplasso objects}
\description{
  Obtains predictions from a \code{grplasso} object. 
}
\usage{
\method{predict}{grplasso}(object, newdata, type = c("link", "response"),
        na.action = na.pass, ...)
}

\arguments{
  \item{object}{a \code{grplasso} object}
  \item{newdata}{\code{data.frame} or design matrix of new observations}
  \item{type}{the type of prediction. \code{type = "link"} is on the
    scale of linear predictors, whereas \code{type = "response"} is on
    the scale of the response variable, i.e. \code{type = "response"}
    applies the inverse link function to the linear predictors.}
  \item{na.action}{function determining what should be done with missing values
    in \code{newdata}. The default is to predict \code{NA}.}
  \item{...}{other options to be passed to the predict function.}
}
\value{
  A matrix whose \emph{columns} correspond to the different values of
  the penalty parameter \code{lambda} of the \code{grplasso} object.
}
\note{
If \code{newdata} is given, offsets specified by \code{offset} in the
fit by \code{grplasso.default} will not be included in predictions,
whereas those specified by an offset term in the formula will be considered.
}
\seealso{\code{\link{grplasso}}}
\examples{
data(splice)

contr <- rep(list("contr.sum"), ncol(splice) - 1)
names(contr) <- names(splice)[-1]

fit <- grplasso(y ~ ., data = splice, model = LogReg(), lambda = 10,
                contrasts = contr, standardize = TRUE)

pred <- predict(fit)
pred.resp <- predict(fit, type = "response")

## The following points should lie on the sigmoid curve
plot(pred, pred.resp)
}
\keyword{methods}
