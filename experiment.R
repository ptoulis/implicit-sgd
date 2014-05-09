# Copyright (c) 2013
# Panos Toulis, ptoulis@fas.harvard.edu
#
library(mvtnorm)  # recall rmvnorm(n,...) returns n x p matrix.

base.learning.rate <- function(t, gamma0, alpha, c) {
  # Computes a learning rate of the form g * (1 + a * t * g)^-c
  #
  # Typically a, g have to be set according to the curvature 
  # of the loss function (log-likelihood)
  # c = is usually problem-independent but may get different values 
  # according to convexity of loss.
  CHECK_TRUE(all(c(gamma0, alpha, c) >= 0), msg="Positive params in learning rate.")
  CHECK_INTERVAL(c, 0, 1, msg="c in [0,1]")
  x = exp(log(gamma0) - c * log(1 + alpha * gamma0 * t))
  y = gamma0 * (1 + alpha * gamma0 * t)^-c
  CHECK_NEAR(x, y, tol=1e-4)
  return(y)
}

glm.score.function <- function(h.transfer, theta, datapoint) {
  # Computes  (yt - h(theta' xt)) * xt = score function
  # for a GLM model with transfer function "h.transfer"
  # 
  # Examples:
  #   normal model : h(x) = x  identity function
  #   poisson model : h(x) = e^x
  #   logistic regression : h(x) = logit(x)
  #
  # Returns: px1 vector of the score (gradient of log-likelihood) 
  xt = datapoint$xt
  yt = datapoint$yt
  # Check dimensions
  CHECK_EQ(length(xt), length(theta))
  CHECK_EQ(length(yt), 1, msg="need one-dimensional outcome")
  yt.hat = h.transfer(sum(xt * theta))
  with(datapoint, matrix((yt - yt.hat) * xt, ncol=1))
}

copy.experiment <- function(experiment) {
  e = empty.experiment(experiment$niters)
  for(i in names(experiment)) {
    e[[i]] = experiment[[i]]
  }
  return(e)
}

empty.experiment <- function(niters) {
  # Returns an empty EXPERIMENT object.
  # Useful for initialization and inspection.
  return(list(name="default",
              theta.star=matrix(0, nrow=1, ncol=1),
              niters=niters,
              score.function=function(theta, datapoint) {},
              sample.dataset=function() {}))
}

get.experiment.description <- function(experiment) {
  # Returns a description of an experiment as a string.
  return(sprintf(" Experiment %s: iters=%d p=%d limit.lr=%.2f",
                 experiment$name, 
                 experiment$niters, experiment$p,
                 experiment$learning.rate(10^9) * 10^9))  
}

get.experiment <- function(name="normal",
                           niters=1000) {
  # Creates an EXPERIMENT object (see terminology) from those 
  # that are available. Use the mnemonic name and the specified niters variable
  function.name = sprintf("%s.experiment", name)
  e = do.call(function.name, args=list(niters=niters))
  e$name = name
  return(e)
}

limit.variance <- function(experiment) {
  # Computes the asymptotic variance from the Theorem of (Toulis et al, 2014)
  J = experiment$J
  limit.a = experiment$learning.rate(10^9) * 10^9
  if(limit.a > 10^3) stop("Error. Learning rate grows indefinitely.")
  I = diag(experiment$p)
  return(limit.a * solve(2 * limit.a * J - I) %*% J)
}

best.alpha <- function(max.alpha=2, nalphas=10^3, p=10) {
  # TODO(ptoulis): Not sure what this is doing atm.
  coeffs = seq(0, max.alpha, length.out=nalphas)
  base.experiment = normal.experiment(p=p, niters=100)
  print(limit.variance(base.experiment))
  print(eigen(limit.variance(base.experiment))$values)
  traces = sapply(coeffs, function(a) {
    experiment = normal.experiment(niters=100, p=p)
    experiment$learning.rate = function(t) a * base.experiment$learning.rate(t)
    sum(diag(limit.variance(experiment)))
  })
  plot(coeffs, traces, type="l")       
  l0 = min(eigen(base.experiment$J)$values)
  x.critical = 0.5 * l0 / eigen(base.experiment$J)$values
  abline(v=x.critical, col="red")
}

sample.covariance.matrix <- function(p) {
  # Samples a low-rank covariance matrix.
  #
  u1 = 0.5 * seq(-1, 1, length.out=p)
  u2 = seq(0.2, 1, length.out=p)
  C = matrix(0, nrow=p, ncol=p)
  diag(C) <- u2
  V =  (C + u1 %*% t(u1))
  CHECK_TRUE(all(eigen(V)$values > 0))
  V
}

normal.experiment <- function(niters, p=100, lr.scale=1.0) {
  # Normal experiment (linear regression)
  #
  # Assume xt ~ N(0, A)  where A has eigenvalues from 0.01 to 1
  #        yt | xt = xt'θ* + ε ,  where ε ~ N(0, 1) ind.
  #
  # Thus the score function is equal to
  #     (yt - xt'θ) * xt  since h(.) transfer = identity
  #
  # Args:
  #   niters = number of samples (also #iterations for online algorithms)
  #   p = #parameters (dimension of the problem)
  #   lr.scale = scale of learning rate
  # 1. Define θ*
  experiment = empty.experiment(niters)
  experiment$name = "normal"
  # experiment$theta.star = matrix(runif(p, min=0, max=5), ncol=1) 
  experiment$theta.star = matrix(rep(1, p), ncol=1)  # all 1's
  experiment$p = p
  A = sample.covariance.matrix(p)
  experiment$Vx = A
  # 2. Define the sample dataset function.
  experiment$sample.dataset = function() {
    epsilon = matrix(rnorm(niters), ncol=1)
    X = rmvnorm(niters, mean=rep(0, p), sigma=A)
    Y = X %*% experiment$theta.star + epsilon
    if(niters > 1000) {
      # CHECK_MU0(as.vector(Y), 0)
    }
    CHECK_TRUE(nrow(X) == niters)
    return(list(X=X, Y=Y))
  }
  
  id.fn = function(x) x
  gamma0 = 1 / sum(diag(A))
  lambda0 = min(eigen(A)$values)
  
  # 3. Define the score function
  experiment$h.transfer <- function(u) u
  experiment$score.function = function(theta, datapoint) {
    glm.score.function(h.transfer=id.fn, theta, datapoint)
  }
  
  # 4. Define the learning rate
  experiment$learning.rate <- function(t) {
    # stop("Need to define learning rate per-application.")
    lr.scale * base.learning.rate(t, gamma0=gamma0, alpha=lambda0, c=1)
  }
  
  # 4b. Fisher information
  experiment$J = A
  
  # 5. Define the risk . This is usually the negative log-likelihood
  truth = experiment$theta.star
  experiment$risk <- function(theta) {
    CHECK_EQ(length(theta), length(truth))
    tmp = 0.5 * t(theta - truth) %*% A %*% (theta - truth)
    CHECK_EQ(nrow(tmp), 1)
    CHECK_EQ(ncol(tmp), 1)
    CHECK_TRUE(all(tmp >= 0))
    return(as.numeric(tmp))
  }

  return(experiment)
}

sparse.sample <- function(nsamples, nsize) {
  ## Will sample n elements from 0 to (n-1) such that 
  ## P(X=i) ~ 1 / (1+i)   marginally
  # Creates some pseudo-sparsity.
  s = seq(0, nsize-1)
  return(sample(s, size=nsamples, replace=T, prob=(1+s)^-1))
}

poisson.experiment <- function(niters, p=100, max.rate=5, lr.scale=1.0) {
  # Poisson regression
  #
  # Assume xt ~ N(0, A)  where A has eigenvalues from 0.01 to 1
  #        yt | xt = Pois(λ) , log(λ) = xt'θ*
  #
  # Thus the score function is equal to
  #     (yt - exp(xt'θ)) * xt
  #
  # Args:
  #   niters = number of samples (also #iterations for online algorithms)
  #   p = #parameters (dimension of the problem)
  #   lr.scale = scale of learning rate
  # 1. Define θ*
  experiment = empty.experiment(niters)
  experiment$name = "poisson"
  experiment$theta.star = matrix(log(c(2, 12)), ncol=1)  # just a bivariate experiment
  experiment$p = 2
  
  sample.X <- function(n) {
    # code=0 then x=(0, 0), code=1 x=(1,0) etc.
    code = sample(0:2, size=n, replace=T, prob=c(8, 1, 1))
    X = matrix(0, nrow=n, ncol=2)
    X[,1] <- as.numeric(code==1)
    X[,2] <- as.numeric(code==2)
    return(X)
  }
  
  ## 1b Set Covariance of X
  empirical.X = sample.X(20000) ## used for some empirical methods.
  CHECK_EQ(sum(apply(empirical.X, 1, prod)), 0, msg="No (1,1) vectors")
  CHECK_MU0(apply(empirical.X, 1, sum)==0, 0.8)
  experiment$Vx = cov(empirical.X)
  
  # 2. Define the sample dataset function.
  experiment$sample.dataset = function() {
    epsilon = matrix(rnorm(niters), ncol=1)
    X = sample.X(niters)
    log.lambdas = X %*% experiment$theta.star
    CHECK_EQ(length(log.lambdas), niters)
    Y = matrix(rpois(niters, lambda=exp(log.lambdas)), ncol=1)
    CHECK_TRUE(nrow(X) == niters)
    return(list(X=X, Y=Y))
  }

  experiment$h.transfer <- function(u) {
    exp(u)
  }
  # 3. Define the score function
  experiment$score.function = function(theta, datapoint) {
    glm.score.function(h.transfer=exp, theta, datapoint)
  }
  
  
  A = experiment$Vx
  gamma0 = 1 / sum(diag(A))
  lambda0 = min(eigen(A)$values)
  # 4. Define the learning rate
  experiment$learning.rate <- function(t) {
    # stop("Need to define learning rate per-application.")
    lr.scale * base.learning.rate(t, gamma0=gamma0, alpha=lambda0, c=1)
  }
  
  # 4b. Fisher information
  # Recall J= E(h'(theta.star' x) x x')
  experiment$J = matrix(0, nrow=2, ncol=2)
  N = nrow(empirical.X)
  for(i in 1:N) {
    x = empirical.X[i, ]
    h.prime = exp(sum(experiment$theta.star * x))
    experiment$J <- experiment$J + h.prime * x %*% t(x)
  }
  experiment$J = experiment$J / N
  CHECK_NEAR(diag(experiment$J), 0.1 * exp(experiment$theta.star), tol=0.05)
  # 5. Define the risk . This is usually the negative log-likelihood
  truth = experiment$theta.star
  experiment$risk <- function(theta) {
    vector.dist(theta, experiment$theta.star)
  }
  
  return(experiment)
}

