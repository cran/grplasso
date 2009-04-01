\name{grplasso}
\alias{grplasso}
\alias{grplasso.default}
\alias{grplasso.formula}
\title{Function to fit a solution of a group lasso problem}
\description{Fits the solution of a group lasso problem for a model of type
  \code{grpl.model}.}  
\synopsis{grplasso(x, ...)}
\usage{
grplasso(x, ...)

\method{grplasso}{formula}(formula, nonpen = ~ 1, data, weights = rep(1, length(y)),
         subset, na.action, lambda, coef.init, penscale = sqrt,
         model = LogReg(), center = TRUE, standardize = TRUE,
         control = grpl.control(), contrasts = NULL, ...)

\method{grplasso}{default}(x, y, index, weights = rep(1, length(y)), offset = rep(0,
         length(y)), lambda, coef.init = rep(0, ncol(x)),
         penscale = sqrt, model = LogReg(), center = TRUE,
         standardize = TRUE, control = grpl.control(), ...)
}
\arguments{
  \item{x}{design matrix (including intercept)}
  \item{y}{response vector}
  \item{formula}{\code{formula} of the penalized variables. The response
    has to be on the left hand side of \code{~}.}
  \item{nonpen}{\code{formula} of the nonpenalized variables. This will
    be added to the \code{formula} argument above and doesn't need to have the
    response on the left hand side.}
  \item{data}{\code{data.frame} containing the variables in the model.}
  \item{index}{vector which defines the grouping of the
    variables. Components sharing the same
    number build a group. Non-penalized coefficients are marked with
    \code{NA}.}
  \item{weights}{vector of observation weights.}
  \item{subset}{an optional vector specifying a subset of observations to be
    used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when the data
    contain 'NA's.}
  \item{offset}{vector of offset values; needs to have the same length as the
    response vector.}
  \item{lambda}{vector of penalty parameters. Optimization starts with
    the first component. See details below.}
  \item{coef.init}{initial vector of parameter estimates corresponding
    to the first component in the vector \code{lambda}.} 
  \item{penscale}{rescaling function to adjust the value of the penalty
    parameter to the degrees of freedom of the parameter group. See the
    reference below.}
  \item{model}{an object of class \code{\link{grpl.model}} implementing
    the negative log-likelihood, gradient, hessian etc. See the
    documentation of \code{\link{grpl.model}} for more details.}
  \item{center}{logical. If true, the columns of the design matrix will be
    centered (except a possible intercept column).}
  \item{standardize}{logical. If true, the design matrix will be
    blockwise orthonormalized such that for each block \eqn{X^TX = n 1}
         (*after* possible centering).}
  \item{control}{options for the fitting algorithm, see
    \code{\link{grpl.control}}.}
  \item{contrasts}{an optional list. See the 'contrasts.arg' of
    'model.matrix.default'.}  
  \item{...}{additional arguments to be passed to the functions defined
    in \code{model}.}
}
\details{
  When using \code{grplasso.formula}, the grouping of the variables is
  derived from the type of the variables: The dummy variables of a
  factor will be automatically treated as a group. 

  The optimization process starts using the first component of
  \code{lambda} as penalty parameter \eqn{\lambda} and with starting
  values defined in \code{coef.init} for the parameter vector. Once
  fitted, the next component of \code{lambda} is considered as penalty
  parameter with starting values defined as the (fitted) coefficient
  vector based on the previous component of \code{lambda}.
}

\value{A \code{grplasso} object is returned, for which \code{coef},
  \code{print}, \code{plot} and \code{predict} methods exist.
  \item{coefficients}{coefficients with respect to the \emph{original} input
    variables (even if \code{standardize = TRUE} is used for fitting).}
  \item{lambda}{vector of lambda values where coefficients were calculated.}
  \item{index}{grouping index vector.}
}
\references{Lukas Meier, Sara van de Geer and Peter B\"uhlmann (2008), \emph{The
    Group Lasso for Logistic Regression}, Journal of the Royal
  Statistical Society, 70 (1), 53 - 71, see also
    \url{http://stat.ethz.ch/~meier/logistic-grouplasso.php}} 
\author{Lukas Meier, \email{meier@stat.math.ethz.ch}}
\examples{
## Use the Logistic Group Lasso on the splice data set
data(splice)

## Define a list with the contrasts of the factors
contr <- rep(list("contr.sum"), ncol(splice) - 1)
names(contr) <- names(splice)[-1]

## Fit a logistic model 
fit.splice <- grplasso(y ~ ., data = splice, model = LogReg(), lambda = 20,
                       contrasts = contr, center = TRUE, standardize = TRUE)

## Perform the Logistic Group Lasso on a random dataset
set.seed(79)

n <- 50  ## observations
p <- 4   ## variables

## First variable (intercept) not penalized, two groups having 2 degrees
## of freedom each

index <- c(NA, 2, 2, 3, 3)

## Create a random design matrix, including the intercept (first column)
x <- cbind(1, matrix(rnorm(p * n), nrow = n))
colnames(x) <- c("Intercept", paste("X", 1:4, sep = ""))

par <- c(0, 2.1, -1.8, 0, 0)
prob <- 1 / (1 + exp(-x \%*\% par))
mean(pmin(prob, 1 - prob)) ## Bayes risk
y <- rbinom(n, size = 1, prob = prob) ## binary response vector

## Use a multiplicative grid for the penalty parameter lambda, starting
## at the maximal lambda value
lambda <- lambdamax(x, y = y, index = index, penscale = sqrt,
                    model = LogReg()) * 0.5^(0:5)

## Fit the solution path on the lambda grid
fit <- grplasso(x, y = y, index = index, lambda = lambda, model = LogReg(),
                penscale = sqrt,
                control = grpl.control(update.hess = "lambda", trace = 0))

## Plot coefficient paths
plot(fit)
}
\keyword{models}
\keyword{regression}
