setClass("grpl.control",
         representation = representation(
           save.x       = "logical",
           save.y       = "logical",
           update.hess  = "character",
           update.every = "numeric",
           tol          = "numeric",
           lower        = "numeric",
           upper        = "numeric",
           beta         = "numeric",
           sigma        = "numeric",
           trace        = "numeric"),

         prototype = list(
           save.x       = FALSE,
           save.y       = TRUE,
           update.hess  = "always",
           update.every = 1,
           tol          = 5 * 10^-8,
           lower        = 10^-2,
           upper        = Inf,
           beta         = 0.5,
           sigma        = 0.1,
           trace        = 1),
           
         validity = function(object){
           if(object@beta <= 0 | object@beta >= 1)
             return("beta has to be in (0, 1)")
           
           if(object@sigma <= 0 | object@sigma >= 1)
             return("sigma has to be in (0, 1)")

           if(object@tol <= 0)
             return("tol has to be positive")

           if(object@lower > object@upper)
             return("lower <= upper has to hold")

           return(TRUE)
         }
)

grpl.control <- function(save.x = FALSE, save.y = TRUE,
                         update.hess = c("always", "lambda"),
                         update.every = 1,
                         tol = 5 * 10^-8, lower = 10^-2, upper = Inf,
                         beta = 0.5, sigma = 0.1, trace = 1){
  
  ## Purpose: Options for the Group Lasso Algorithm
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## save.x: a logical indicating whether the design matrix should be saved.
  ## save.y: a logical indicating whether the response should be saved.
  ## update.hess: should the hessian be updated in each
  ##              iteration ("always")? update.hess = "lambda" will update
  ##              the Hessian once for each component of the penalty
  ##              parameter "lambda" based on the parameter estimates
  ##              corresponding to the previous value of the penalty
  ##              parameter. 
  ## tol: convergence tolerance; the smaller the more precise, see
  ##      details below.
               
  ## lower: lower bound for the diagonal approximation of the
  ##        corresponding block submatrix of the Hessian of the negative
  ##        log-likelihood function.
  ## upper: upper bound for the diagonal approximation of the
  ##        corresponding block submatrix of the Hessian of the negative
  ##        log-likelihood function.
  ## beta: scaling factor beta < 1 of the Armijo line search.
  ## sigma: 0 < \sigma < 1 used in the Armijo line search.
  ## trace: integer. "0" omits any output,
  ##        "1" prints the current lambda value,
  ##        "2" prints the improvement in the objective function after each
  ##        sweep through all the parameter groups and additional
  ##        information.
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  1 Jun 2006, 10:02

  
  update.hess <- match.arg(update.hess)

  RET <- new("grpl.control",
             save.x       = save.x,
             save.y       = save.y,
             update.hess  = update.hess,
             update.every = update.every,
             tol          = tol,
             lower        = lower,
             upper        = upper,
             beta         = beta,
             sigma        = sigma,
             trace        = trace)
  RET
}
