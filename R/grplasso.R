grplasso <- function(x, ...)
  UseMethod("grplasso")

grplasso.formula <- function(formula, nonpen = ~ 1, data,
                             weights, subset, na.action,
                             lambda, coef.init,
                             penscale = sqrt, model = LogReg(),
                             standardize = TRUE,
                             control = grpl.control(),
                             contrasts = NULL, ...){
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 27 Jun 2006, 14:52

  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  
  ## Remove not-needed stuff to create the model-frame
  m$nonpen <- m$lambda <- m$coef.init <- m$penscale <- m$model <-
    m$standardize <- m$contrasts <- m$control <- m$... <- NULL

  l <- create.design(m, formula, nonpen, data, weights, subset, na.action,
                     contrasts, parent.frame())
  
  if(missing(coef.init))
    coef.init <- rep(0, ncol(l$x))
  
  fit <- grplasso.default(x = l$x, y = l$y, index = l$index, weights = l$w,
                          offset = l$off, lambda = lambda,
                          coef.init = coef.init,
                          penscale = penscale, model = model,
                          standardize = standardize, control = control,
                          ...)
  fit$terms <- l$Terms
  fit$contrasts <- attr(l$x, "contrasts")
  fit$xlevels <- .getXlevels(l$Terms, l$mf)
  fit$na.action <- attr(l$mf, "na.action")
  fit$call <- match.call() ## Overwrite grplasso.default 
  structure(fit, class = "grplasso")
}

grplasso.default <- function(x, y, index, weights = rep(1, length(y)),
                             offset = rep(0, length(y)), lambda,
                             coef.init = rep(0, ncol(x)),
                             penscale = sqrt, model = LogReg(),
                             standardize = TRUE, control = grpl.control(), ...)
{
  ## Purpose: Function to fit a solution (path) of a group lasso problem
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## x: design matrix (including intercept), already rescaled and
  ##    possibly blockwise orthonormalized.
  ## y: response vector
  ## index: vector which defines the grouping of the variables. Components
  ##        sharing the same number build a group. Non-penalized
  ##        coefficients are marked with "NA".
  ## weights: vector of observation weights.
  ## offset: vector of offset values; needs to have the same length as the
  ##         response vector.
  ## lambda: vector of penalty parameters. Optimization starts with the
  ##         first component. See details below.
  ## coef.init: initial vector of parameter estimates corresponding to the
  ##            first component in the vector "lambda". 
  ## penscale: rescaling function to adjust the value of the penalty
  ##           parameter to the degrees of freedom of the parameter group.
  ##           See the reference below.
  ## model: an object of class "grpl.model" implementing the negative
  ##        log-likelihood, gradient, hessian etc. See the documentation
  ##        of "grpl.model" for more details.
  ## control: options for the fitting algorithm, see "grpl.control".
  ## ...: additional arguments to be passed to the functions defined
  ##      in "model".
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 30 Aug 2005, 09:02

  ## Do some error checking first

  ## Check the design matrix
  if(!is.matrix(x))
    stop("x has to be a matrix")

  if(any(is.na(x)))
    stop("Missing values in x not allowed!")

  ## Check the response
  if(!is.numeric(y))
    stop("y has to be of type 'numeric'")

  if(!model@check(y))
    stop("y has wrong format")

  ## Check the other arguments
  if(length(weights) != length(y))
    stop("length(weights) not equal length(y)")

  if(any(weights < 0))
    stop("Negative weights not allowed")
  
  if(length(offset) != length(y))
    stop("length(offset) not equal length(y)")

  if(length(coef.init) != ncol(x))
    stop("length(coef.init) not equal ncol(x)")

  if(!is.numeric(index))
    stop("index has to be of type 'numeric'!")

  if(is.unsorted(rev(lambda)))
    warning("lambda values should be sorted in decreasing order")

  if(all(is.na(index)))
    stop("None of the predictors are penalized.")
  
  check <- validObject(control) ## will stop the program if error occurs
  
  ## Extract the control information
  update.hess  <- control@update.hess
  update.every <- control@update.every
  inner.loops  <- control@inner.loops
  lower        <- control@lower
  upper        <- control@upper
  save.x       <- control@save.x
  save.y       <- control@save.y
  tol          <- control@tol
  trace        <- control@trace
  beta         <- control@beta
  sigma        <- control@sigma

  nrlambda <- length(lambda)
  ncolx    <- ncol(x)
  nrowx    <- nrow(x)
  
  ## Which are the non-penalized parameters?
  any.notpen    <- any(is.na(index))
  inotpen.which <- which(is.na(index))
  nrnotpen      <- length(inotpen.which)
  
  ## Index vector of the penalized parameter groups
  if(any.notpen){
    ipen <- index[-inotpen.which]
    ipen.which <- split((1:ncolx)[-inotpen.which], ipen)
  }else{
    cat("\n...All groups are penalized. Did you include an intercept in your\n")
    cat("   design matrix and really want to penalize it?\n")
    ipen <- index
    ipen.which <- split((1:ncolx), ipen)
  }

  nrpen      <- length(ipen.which)
  dict.pen   <- sort(unique(ipen))
  
  ## Table of degrees of freedom
  ipen.tab   <- table(ipen)[as.character(dict.pen)]
  
  x.old <- x
  ## Standardize the design matrix -> blockwise orthonormalization
  if(standardize){
    cat("...Using standardized design matrix.\n")
    stand        <- blockstand(x, ipen.which, inotpen.which)
    x            <- stand$x
    scale.pen    <- stand$scale.pen
    scale.notpen <- stand$scale.notpen
  }
  ## From now on x is the *normalized* design matrix!
  
  ## Extract the columns into lists, works faster for large matrices
  if(any.notpen){
    x.notpen <- list(); length(x.notpen) <- nrnotpen
    for(i in 1:length(inotpen.which))
      x.notpen[[i]] <- x[,inotpen.which[[i]], drop = FALSE]
  }
  
  x.pen <- list(); length(x.pen) <- length(nrpen)
  for(i in 1:length(ipen.which))
    x.pen[[i]] <- x[,ipen.which[[i]], drop = FALSE]

  ## Extract the needed functions
  check     <- validObject(model)
  invlink   <- model@invlink
  nloglik   <- model@nloglik
  ngradient <- model@ngradient
  nhessian  <- model@nhessian

  #########################################################################
  ##                                                                     ##
  ## Start the optimization process                                      ##
  ##                                                                     ##
  #########################################################################

  coef      <- coef.init
  coef.pen  <- coef.init
  if(any.notpen)
    coef.pen  <- coef[-inotpen.which]

  norms.pen <- c(sqrt(rowsum(coef.pen^2, group = ipen)))

  norms.pen.m  <- matrix(0, nrow = nrpen, ncol = nrlambda,
                         dimnames = list(NULL, lambda))
  norms.npen.m <- matrix(0, nrow = nrnotpen, ncol = nrlambda,
                         dimnames = list(NULL, lambda))
  nloglik.v <- fn.val.v <- numeric(nrlambda)
  coef.m    <- grad.m <- matrix(0, nrow = ncolx, ncol = nrlambda,
                                dimnames = list(colnames(x), lambda))
  fitted    <- linear.predictors <- matrix(0, nrow = nrowx,
                                           ncol = nrlambda,
                                           dimnames = list(rownames(x), lambda))

  ## *Initial* vector of linear predictors (eta) and transformed to the
  ## scale of the response (mu)
  eta <- offset + c(x %*% coef)
  mu <- invlink(eta)

  if(any.notpen){
    nH.notpen <- numeric(nrnotpen)
  }
  nH.pen <- numeric(nrpen)

  for(pos in 1:nrlambda){
    l <- lambda[pos]

    if(trace >= 2)
      cat("\nLambda:", l, "\n")

    ## Initial Hessian Matrix of the *negative* log-likelihood function
    ## (uses parameter estimates based on the last penalty parameter value)

    ##if(!(update.hess == "always") | pos == 1){
    if(update.hess == "lambda" & pos %% update.every == 0 | pos == 1){
      if(any.notpen){
        for(j in inotpen.which){
          Xj <- x.notpen[[j]] 
          nH.notpen[j] <- min(max(nhessian(Xj, mu, weights, ...), lower), upper)
        }
      }
      
      for(j in 1:nrpen){
        ind <- ipen.which[[j]]
        Xj <- x.pen[[j]] 
        diagH <- numeric(length(ind))
        for(i in 1:length(ind))
          diagH[i] <- nhessian(Xj[, i, drop = FALSE], mu, weights, ...)
        
        nH.pen[j] <- min(max(diagH, lower), upper)
      }
    }
    
    ## Start the optimization process
    fn.val <- nloglik(y, eta, weights, ...) +
      l * sum(penscale(ipen.tab) * norms.pen)

    ## These are needed to get into the while loop the first time
    do.all  <- FALSE
    d.fn <- d.par <- 1

    counter <- 1
    
    while(d.fn > tol | d.par > sqrt(tol) | !do.all){
      fn.val.old <- fn.val
      coef.old   <- coef

      ## Count how many times we are already optimizing
      ## Will be reset (see end of while loop) when subgroup
      ## is finished with optimizing

      ## Check whether we have some useful information from the previous step
       
      if(counter == 0 | counter > inner.loops){
        do.all <- TRUE
        guessed.active <- 1:nrpen
        counter <- 1
        if(trace >= 2)
          cat("...Running through all groups\n")
      }else{
        guessed.active <- which(norms.pen != 0)
        if(length(guessed.active) == 0){ 
          guessed.active <- 1:nrpen
          do.all <- TRUE
          if(trace >= 2)
            cat("...Running through all groups\n")
        }else{
          do.all <- FALSE
          if(counter == 1 & trace >= 2)
            cat("...Starting inner loop\n")
          counter <- counter + 1
        }
      }

      ## These are used for the line search
      start.notpen <- rep(1, nrnotpen)
      start.pen    <- rep(1, nrpen)

      if(any.notpen){
        ## Optimize the *non-penalized* parameters
        for(j in 1:nrnotpen){
          ind <- inotpen.which[j]
          Xj <- x.notpen[[j]]
        
          ## Gradient of the negative log-likelihood function 
          ngrad <- c(ngradient(Xj, y, mu, weights, ...))
          
          if(update.hess == "always"){
            nH <- min(max(nhessian(Xj, mu, weights, ...), lower), upper)
          }else{
            nH <- nH.notpen[j]
          }
          
          d <- -(1 / nH) * ngrad
          d <- zapsmall(c(coef[ind], d))[2]
        
          if(d != 0){
            qh    <- sum(ngrad * d)
            scale <- min(start.notpen[j] / beta, 1) ##1 
            coef.test <- coef
            coef.test[ind] <- coef[ind] + scale * d
            Xjd <- Xj * d
            eta.test     <- eta + Xjd * scale
            fn.val0      <- nloglik(y, eta, weights, ...)
            fn.val.test  <- nloglik(y, eta.test, weights, ...)
            while(fn.val.test - fn.val0 > sigma * scale * qh){
              scale          <- scale * beta
              coef.test[ind] <- coef[ind] + scale * d
              eta.test       <- eta + Xjd * scale
              fn.val.test    <- nloglik(y, eta.test, weights, ...)
            }
            start.notpen[j] <- scale
            coef <- coef.test
            eta  <- eta.test
            mu   <- invlink(eta)
          }
        } 
      }

      ## Optimize the *penalized* parameter groups
      for(j in guessed.active){ ##for(j in 1:nrpen){
        ind  <- ipen.which[[j]]
        npar <- ipen.tab[j]

        coef.ind <- coef[ind]
        cross.coef.ind <- crossprod(coef.ind)

        ## Design matrix of the current group
        Xj <- x.pen[[j]] 


        ngrad <- c(ngradient(Xj, y, mu, weights, ...))
        if(update.hess == "always"){
          diagH <- numeric(length(ind))
          for(i in 1:length(ind)) ## for loop seems to be faster than sapply
            diagH[i] <- nhessian(Xj[,i,drop = FALSE], mu, weights, ...)
          
          nH <- min(max(diagH, lower), upper)
        }else{
          nH <- nH.pen[j]
        }
        cond <- -ngrad + nH * coef.ind 
        cond.norm2 <- crossprod(cond) #sum(cond^2)
        
        ## Check whether the Minimum is at c(0, 0, ...) via the condition on
        ## the subgradient.
        border <- penscale(npar) * l
        if(cond.norm2 > border^2){
          d <- (1 / nH) *
            (-ngrad - l * penscale(npar) * (cond / sqrt(cond.norm2)))
        }else{
          d <- -coef.ind
        }
        if(!all(d == 0)){
          scale <- min(start.pen[j] / beta, 1) ##1 
          qh <- sum(ngrad * d) + 
            l * penscale(npar) * sqrt(crossprod(coef.ind + d)) -
              l * penscale(npar)* sqrt(cross.coef.ind)

          coef.test      <- coef
          coef.test[ind] <- coef.ind + scale * d
          Xjd            <- c(Xj %*% d)
          eta.test       <- eta + Xjd * scale
          fn.val.test    <- nloglik(y, eta.test, weights, ...)
          fn.val0        <- nloglik(y, eta, weights, ...)

          while(fn.val.test - fn.val0 +
                l  * penscale(npar) * sqrt(crossprod(coef.test[ind])) -
                l  * penscale(npar) * sqrt(cross.coef.ind) >
                sigma * scale * qh){
            scale          <- scale * beta
            coef.test[ind] <- coef.ind + scale * d
            eta.test       <- eta + Xjd * scale
            fn.val.test    <- nloglik(y, eta.test, weights, ...)
          }
          coef <- coef.test
          eta  <- eta.test
          mu   <- invlink(eta)
          start.pen[j] <- scale
        }
        norms.pen[j] <- sqrt(crossprod(coef[ind]))
      }
      
      fn.val <- nloglik(y, eta, weights, ...) +
        l * sum(penscale(ipen.tab) * norms.pen)
      
      ## Relative convergence criteria
      d.par <- sqrt(crossprod(coef - coef.old)) / (1 + sqrt(crossprod(coef)))

      d.fn <- (fn.val.old - fn.val) / (1 + abs(fn.val))
      if(trace >= 2){
        cat("d.fn:", d.fn, " d.par:", d.par, " nr.var:",
            sum(coef != 0), "\n")
      }

      ## Check whether the sub-problem is already finished
      if(d.fn <= tol & d.par <= sqrt(tol)){
        counter <- 0 ## will force a run through all groups
        if(trace >= 2 & !do.all)
          cat("...Subproblem (active set) solved\n")
      } 
        
    } ## end of while loop

    if(trace == 1)
      cat("Lambda:", l, " nr.var:", sum(coef != 0), "\n")
    
    coef.m[,pos]            <- coef
    fn.val.v[pos]           <- fn.val
    norms.pen.m[,pos]       <- norms.pen
    nloglik.v[pos]          <- nloglik(y, eta, weights, ...)
    grad.m[,pos]            <- ngradient(x, y, mu, weights, ...)
    linear.predictors[,pos] <- eta
    fitted[,pos]            <- invlink(eta)
  }

  ## Transform the coefficients back to the original scale
  if(standardize){
    if(any.notpen)
      coef.m[inotpen.which,] <- scale.notpen * coef.m[inotpen.which,] 
    for(j in 1:length(ipen.which)){
      ind <- ipen.which[[j]]
      coef.m[ind,] <- solve(scale.pen[[j]], coef.m[ind,,drop = FALSE])
    }
  }
  
  if(!save.x)
    x.old <- NULL
  if(!save.y)
    y <- NULL

  out <- list(x = x.old, ## use untransformed values
              y = y, 
              coefficients = coef.m,
              norms.pen    = norms.pen.m,
              lambda       = lambda,
              index        = index,
              penscale     = penscale,
              model        = model,
              ngradient    = grad.m,
              nloglik      = nloglik.v,
              fitted       = fitted,
              linear.predictors = linear.predictors,
              fn.val   = fn.val.v,
              weights  = weights,
              offset   = offset,
              control  = control,
              call     = match.call())
  structure(out, class = "grplasso")
}


