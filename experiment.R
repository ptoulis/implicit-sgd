# Copyright (c) 2013
# Panos Toulis, ptoulis@fas.harvard.edu
#
library(mvtnorm)  # recall rmvnorm(n,...) returns n x p matrix.

base.learning.rate <- function(t, gamma0, alpha, c) {
  # Computes a learning rate of the form g * (1 + a * t * g)^-c
  #
  # Typically a, g have to be set according to the curvature of the loss function (log-likelihood)
  # c = is usually problem-independent but may get different values according to convexity of loss
  #
  CHECK_TRUE(all(c(gamma0, alpha, c) >= 0), msg="Positive params in learning rate.")
  CHECK_INTERVAL(c, 0, 1, msg="c in [0,1]")
  x = exp(log(gamma0) - c * log(1 + alpha * gamma0 * t))
  y = gamma0 * (1 + alpha * gamma0 * t)^-c
  CHECK_NEAR(x, y, tol=1e-4)
  return(y)
}

glm.score.function <- function(h.transfer, theta, datapoint) {
  # Computes  (yt - h(theta' xt) * xt) = score function
  # for a GLM model with transfer function "h.transfer"
  # 
  # Examples:
  #   normal model : h(x) = x  identity function
  #   poisson model : h(x) = e^x
  #   logistic regression : h(x) = logit(x)
  # 
  xt = datapoint$xt
  yt = datapoint$yt
  CHECK_EQ(length(xt), length(theta))
  CHECK_EQ(length(yt), 1, msg="one-dimensional")
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
  # returns an empty EXPERIMENT object.
  # Useful for initialization and inspection.
  return(list(name="default",
              theta.star=matrix(0, nrow=1, ncol=1),
              niters=niters,
              score.function=function(theta, datapoint) {},
              sample.dataset=function() {}))
}

get.experiment.description <- function(experiment) {
  return(sprintf(" Experiment %s: iters=%d p=%d limit.lr=%.2f",
                 experiment$name, 
                 experiment$niters, experiment$p,
                 experiment$learning.rate(10^9) * 10^9))  
}

get.experiment <- function(name="normal",
                           niters=1000) {
  # Creates an EXPERIMENT object (see terminology)
  # from the mnemonic name and the specified niters variable
  function.name = sprintf("%s.experiment", name)
  e = do.call(function.name, args=list(niters=niters))
  e$name = name
  return(e)
}

limit.variance <- function(experiment) {
  J = experiment$J
  limit.a = experiment$learning.rate(10^9) * 10^9
  if(limit.a > 10^6) stop("There is a problem most likely")
  I = diag(experiment$p)
  return(limit.a^2 * solve(2 * limit.a * J - I) %*% J)
}

best.alpha <- function(max.alpha=2, nalphas=10^3, p=10) {
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

covariance.matrix <- function(p) {
  u1 = seq(0.3, 1.2, length.out=p)
  u2 = seq(0.2, 1, length.out=p)
  V =  (diag(u1) + u2 %*% t(u2))
  CHECK_TRUE(all(eigen(V)$values > 0))
  V
}

normal.experiment <- function(niters, p=100, lr.scale=1.0) {
  # Normal experiment (linear regression)
  # Defined in Xu (2011), Section 6.2, p.8
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
  A = covariance.matrix(p)
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
