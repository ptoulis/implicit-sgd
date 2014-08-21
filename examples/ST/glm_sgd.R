# Copyright (c) 2014
# Panos Toulis, ptoulis@fas.harvard.edu
#
# Explicit and implicit SGD implementation for fitting GLM.
# Reference paper:
#  "Statistical analysis of stochastic gradient methods for generalized linear models"
#   Panagiotis Toulis, Edoardo Airoldi, Jason Rennie ; JMLR W&CP 32 (1) : 667–675, 2014
# 
rm(list=ls())
library(mvtnorm)

simulate.data <- function(p=2, nsamples=1000) {
  Sigma = diag(p)
  x = rmvnorm(nsamples, mean=rep(0, p), sigma = Sigma)
  beta = matrix(1:p, nrow=p, ncol=1)
  y = x %*% beta + rnorm(nsamples)
  
  d = data.frame(x=x, y = y)
}

family_h <- function(family) {
  fh = list(h=list(), h.prime=list())
  if(family=="gaussian") {
    fh = list(h=function(x) x, h.prime=function(x) 1)
  } else if(family=="binomial") {
    fh = list(h=function(x) exp(x) / (1 + exp(x)), h.prime=function(x) h(x) * (1-h(x)))
  } else if (family=="poisson") {
    fh = list(h=function(x) exp(x), h.prime=function(x) exp(x))
  } else if(family=="quasi") {
    # ..
  } else {
    stop(sprintf("Family %s is not supported. Please send a report to ptoulis@fas.harvard.edu", family))
  }
  return(fh)
}

optimal.alpha.lr <- function(fisher.info) {
  lam = eigen(fisher.info)$values
  lam = Re(sqrt(lam * Conj(lam)))
  optim(par=c(0), fn=function(a) sum(a^2 * lam / (2 * a * lam - 1)), 
        lower = 1 / (2 * min(lam)) + 1e-8, upper = 1e6 / min(lam),
        method="L-BFGS-B")$par
}

glm_sgd <- function(formula, family = "gaussian", data, 
                    sgd.params = list(alpha=1),
                    method = "implicit") {
  # Implementation of the SGD methods for GLM models.
  # Default is implicit stochastic gradient descent.
  #
  call <- match.call()
  glm.model <- list() # Transfer function and first derivative.
  
  if (is.character(family)) {
    glm.model = family_h(family)
  } else if(is.list(family) && identical(names(family), c("h", "h.prime"))) {
    glm.model = family
  } else {
    print(family)
    stop("Argument 'family' is not valid.")
  }

  if (missing(data)) {
    data <- environment(formula)
  }
  if (!is.character(method) || (!is.element(method, c("explicit", "implicit")))) 
    stop("invalid SGD 'method' argument")

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "control"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  
  # Read model.
  Y <- model.response(mf, "any")
  X <- if (!is.empty.model(mt)) 
    model.matrix(mt, mf, contrasts)
  else matrix(, nrow(Y), 0L)
  
   # Helper functions for implicit SGD only.
  implicit.coeff <- function(ksi, yn, eta_n, norm_xn) {
    yn - glm.model$h(eta_n + norm_xn * ksi)
  }
  
  solve.implicit <- function(yn, eta_n, norm_xn, an) {
    rn = an * implicit.coeff(0, yn, eta_n, norm_xn)
    Bn = c(0, 0) # search bounds.
    if(rn > 0) {
      Bn = c(0, rn)
    } else {
      Bn <- c(rn, 0)
    }
    ksi_n = NA
    if(Bn[2] != Bn[1]) {
      ksi_n = uniroot(function(ksi) ksi - an * implicit.coeff(ksi, yn, eta_n, norm_xn), 
                      interval=Bn)$root
    } else {
      ksi_n = Bt[1]
    }
    return(ksi_n)
  }
  
  # Control parameters
  all.iters = 1:length(Y)
  alpha.lr = sgd.params$alpha
  an.list = alpha.lr / (alpha.lr + all.iters) # TODO(ptoulis): User-defined learning rates?
  approx.fisher = ncol(X) * diag(ncol(X))  # TODO(ptoulis): Estimating the Fisher matrix is costly atm.
  approx.fisher.updates = 1
  
  # Save iterations. TODO(ptoulis): User-defined?
  nmax.iters.store = max(10, length(Y) / 1000)
  store.points = as.integer(seq(10, length(Y), length.out=nmax.iters.store))
  theta.iters = matrix(NA, nrow=nmax.iters.store, ncol=ncol(X))
  
  # Current parameter estimate.
  theta = rep(0, ncol(X))  # TODO(ptoulis): Initial value should be user-defined?
    
  for(iter in all.iters) {
    # 1. Get n-th data point
    an = an.list[iter]
    yn = Y[iter]
    xn = X[iter, ]
    eta_n = sum(xn * theta)
    yn_pred = glm.model$h(eta_n)
    # 2. Perform the update.
    if(method=="explicit") {
      theta = theta + an * (yn - yn_pred) * xn # TODO(ptoulis): Keep track of estimates?
    } else if(method=="implicit") {
      norm_xn = sum(xn^2)
      ksi_n = solve.implicit(yn, eta_n, norm_xn, an)
      theta = theta + ksi_n * xn
    } else {
      stop(sprintf("Method %s is not recognized.", method))
    }
    # 3. Save the estimate in a matrix, and update estimated Fisher information.
    #    TODO(ptoulis): Should be user-defined.
    index = match(iter, store.points)
    if(!is.na(index)) {
      approx.fisher.updates <- approx.fisher.updates + 1
      # TODO(ptoulis): Should be user-controlled.
      approx.fisher = approx.fisher + glm.model$h.prime(eta_n) * xn %*% t(xn)
      
      # TODO(ptoulis): Output should be user-defined.
      # print(theta)
      print(sprintf("(%d/%d) Fisher updates=%d, Best learning rate %.3f. Current rate is %.3f", 
                    iter, length(all.iters), approx.fisher.updates,
                    optimal.alpha.lr(fisher.info = (1/approx.fisher.updates) * approx.fisher),
                    alpha.lr))
      theta.iters[index, ] <- theta
    }
  }
  # TODO(ptoulis): Output should have estimated covariance matrix,
  #   and thus include standard errors.
  # Output should be similar to output of glm().
  return(theta.iters)  
}
