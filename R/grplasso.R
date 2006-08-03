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
  if(!is.matrix(x))
    stop("x has to be a matrix")

  if(!is.numeric(y))
    stop("y has to be of type 'numeric'")
  
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
  update.hess <- control@update.hess
  lower       <- control@lower
  upper       <- control@upper
  save.x      <- control@save.x
  save.y      <- control@save.y
  tol         <- control@tol
  trace       <- control@trace
  beta        <- control@beta
  sigma       <- control@sigma
  
  ## Which are the non-penalized parameters?
  any.notpen <- any(is.na(index))

  inotpen.which <- which(is.na(index))
  
  ## Index vector of the penalized parameter groups
  ipen <- index[!is.na(index)]

  dict.pen <- sort(unique(ipen))
  ipen.tab <- table(ipen)[as.character(dict.pen)] ## Table of degrees of freedom
  ## Indices of parameter groups
  ipen.which <- list(); length(ipen.which) <- length(dict.pen)
  for(j in 1:length(dict.pen))
    ipen.which[[j]] <- which(index == dict.pen[j])

  check     <- validObject(model)
  invlink   <- model@invlink
  nloglik   <- model@nloglik
  ngradient <- model@ngradient
  nhessian  <- model@nhessian

  x.old <- x
  ## Standardize the design matrix -> blockwise orthonormalization
  if(standardize){
    stand        <- blockstand(x, ipen.which, inotpen.which)
    x            <- stand$x
    scale.pen    <- stand$scale.pen
    scale.notpen <- stand$scale.notpen
  }
  
  #########################################################################
  ##                                                                     ##
  ## Start the optimization process                                      ##
  ##                                                                     ##
  #########################################################################

  coef <- coef.init
  
  norms.pen    <- sqrt(sapply(1:length(ipen.which),
                              function(j) crossprod(coef[ipen.which[[j]]])))
  norms.pen.m  <- matrix(0, nrow = length(norms.pen), ncol = length(lambda),
                         dimnames = list(NULL, lambda))
  norms.npen.m <- matrix(0, nrow = sum(is.na(index)), ncol = length(lambda),
                         dimnames = list(NULL, lambda))
  nloglik.v <- fn.val.v <- numeric(length(lambda))
  coef.m    <- grad.m <- matrix(0, nrow = ncol(x), ncol = length(lambda),
                                dimnames = list(colnames(x), lambda))
  fitted    <- linear.predictors <- matrix(0, nrow = nrow(x),
                                           ncol = length(lambda),
                                           dimnames = list(rownames(x), lambda))

  ## *Initial* vector of linear predictors (eta) and transformed to the
  ## scale of the response (mu)
  eta <- offset + c(x %*% coef)
  mu <- invlink(eta)

  for(pos in 1:length(lambda)){
    l <- lambda[pos]
    if(trace >= 1)
      cat("Lambda: ", l, "\n")

    d.fn <- d.par <- 1

    ## Initial Hessian Matrix of the *negative* log-likelihood function
    ## (uses parameter estimates based on the last penalty parameter value)
    ## If we update the Hessian we don't
    if(!(update.hess == "always") | pos == 1){
      if(any.notpen){
        nH.notpen <- numeric()
        for(j in inotpen.which){
          Xj <- x[,j, drop = FALSE]
          nH.notpen[j] <- min(max(nhessian(Xj, mu, weights, ...), lower), upper)
        }
      }
      
      nH.pen <- numeric()
      for(j in 1:length(ipen.which)){
        ind <- ipen.which[[j]]
        Xj <- x[,ind, drop = FALSE] ##ipen.which[[j]]]
        diagH <- numeric(length(ind))
        for(i in 1:length(ind))
          diagH[i] <- nhessian(Xj[, i, drop = FALSE], mu, weights, ...)
        
        nH.pen[j] <- min(max(diagH, lower), upper)
      }
    }
    
    ## Start the optimization process
    fn.val <- nloglik(y, eta, weights, ...) +
      l * sum(penscale(ipen.tab) * norms.pen)
    
    while(d.fn > tol | d.par > sqrt(tol)){
      fn.val.old <- fn.val
      coef.old   <- coef

      if(any.notpen){
        ## Optimize the *non-penalized* parameters
        for(j in 1:length(inotpen.which)){
          ind <- inotpen.which[j]
          ##mu <- invlink(eta) ## to be removed
          Xj <- x[,ind,drop = FALSE]
        
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
            scale <- 1 
            coef.test <- coef
            coef.test[ind] <- coef[ind] + scale * d
            Xjd <- Xj * d
            eta.test     <- eta + Xjd * scale
            fn.val0      <- nloglik(y, eta, weights, ...)
            fn.val.test  <- nloglik(y, eta.test, weights, ...)
            while(fn.val.test - fn.val0 > sigma * scale * qh){
              scale <- scale * beta
              coef.test[ind] <- coef[ind] + scale * d
              eta.test <- eta + Xjd * scale
              fn.val.test <- nloglik(y, eta.test, weights, ...)
            }
            coef <- coef.test
            eta  <- eta.test
            mu <- invlink(eta)
          }
        }
      }
      
      ## Optimize the *penalized* parameter groups
      for(j in 1:length(ipen.which)){
        ind  <- ipen.which[[j]]
        npar <- ipen.tab[j]

        ## Design matrix of the current group
        Xj <- x[,ind, drop = FALSE]

        ngrad <- c(ngradient(Xj, y, mu, weights, ...))
        if(update.hess == "always"){
          diagH <- numeric(length(ind))
          for(i in 1:length(ind)) ## for loop seems to be faster than sapply
            diagH[i] <- nhessian(Xj[,i,drop = FALSE], mu, weights, ...)
          
          nH <- min(max(diagH, lower), upper)
        }else{
          nH <- nH.pen[j]
        }
        cond <- -ngrad + nH * coef[ind] 
        cond.norm2 <- crossprod(cond) #sum(cond^2)
        
        ## Check whether the Minimum is at c(0, 0, ...) via the condition on
        ## the subgradient.
        border <- penscale(npar) * l
        if(cond.norm2 > border^2){
          d <- (1 / nH) *
            (-ngrad - l * penscale(npar) * (cond / sqrt(cond.norm2)))
        }else{
          d <- -coef[ind]
        }
        if(!all(d == 0)){
          scale <- 1 
          qh <- sum(ngrad * d) + 
            l * penscale(npar) * sqrt(crossprod(coef[ind] + d)) -
              l * penscale(npar)* sqrt(crossprod(coef[ind]))

          coef.test      <- coef
          coef.test[ind] <- coef[ind] + scale * d
          Xjd            <- c(Xj %*% d)
          eta.test       <- eta + Xjd * scale
          fn.val.test    <- nloglik(y, eta.test, weights, ...)
          fn.val0        <- nloglik(y, eta, weights, ...)

          while(fn.val.test - fn.val0 +
                l  * penscale(npar) * sqrt(crossprod(coef.test[ind])) -
                l  * penscale(npar) * sqrt(crossprod(coef[ind])) >
                sigma * scale * qh){
            scale <- scale * beta
            coef.test[ind] <- coef[ind] + scale * d
            eta.test <- eta + Xjd * scale
            fn.val.test <- nloglik(y, eta.test, weights, ...)
          }
          coef <- coef.test
          eta <- eta.test
          mu <- invlink(eta)
        }
        norms.pen[j] <- sqrt(crossprod(coef[ind]))
      }

      fn.val <- nloglik(y, eta, weights, ...) +
        l * sum(penscale(ipen.tab) * norms.pen)
      
      ## Relative convergence criteria
      d.par <- sqrt(crossprod(coef - coef.old)) / (1 + sqrt(crossprod(coef)))

      d.fn <- (fn.val.old - fn.val) / (1 + abs(fn.val))
      if(trace >= 2){
        cat("d.fn: ", d.fn, "d.par: ", d.par, "nr.var: ",
            sum(coef != 0), "\n")
      }
    }
    coef.m[,pos]            <- coef
    fn.val.v[pos]           <- fn.val
    #norms.npen.m[,pos]      <- abs(coef[inotpen.which])
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
              norms.pen = norms.pen.m,
              lambda = lambda,
              index = index,
              penscale = penscale,
              model = model,
              ngradient = grad.m,
              nloglik = nloglik.v,
              fitted = fitted,
              linear.predictors = linear.predictors,
              fn.val = fn.val.v,
              weights = weights,
              offset = offset,
              control = control,
              call = match.call())
  structure(out, class = "grplasso")
}


