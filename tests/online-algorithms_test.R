## Unit-tests for the online algorithms module.
source("online-algorithms.R")

test.sgd <- function() {
  p = 10
  n = 100
  e = normal.experiment(niters=n, p=p)
  e$learning.rate <- function(t) {
    1 / (t+1)
  }
  # Create dataset such that SGD (for the specified learning rate)
  # will give estimates that will satisfy:
  #  theta_t = x1 + x2 + x3 + ....x_t
  d = list()
  d$X = matrix(rpois(n * p, lambda=10), nrow=n, ncol=p, byrow=T)
  xsums = t(apply(d$X, 2, function(s) cumsum(s)))
  y = sapply(1:nrow(d$X), function(i) {
    if (i > 1)  {
      return(i+1 + sum(d$X[i,] * colSums(matrix(d$X[1:(i-1), ], ncol=e$p))))
    } else {
      return(2)
    }})
  d$Y = matrix(y, ncol=1)
  
  out = run.onlineAlgorithm(d, e, algorithm=sgd.onlineAlgorithm)
  CHECK_TRUE(all(out$estimates == xsums))
}

test.asgd <- function() {
  p = 4
  n = 20
  e = normal.experiment(niters=n, p=p)
  e$learning.rate <- function(t) {
    1 / (t+1)
  }
  # Create dataset such that SGD (for the specified learning rate)
  # will give estimates that will satisfy:
  #  theta_t = x1 + x2 + x3 + ....x_t
  d = list()
  d$X = matrix(rpois(n * p, lambda=10), nrow=n, ncol=p, byrow=T)
  xsums = t(apply(d$X, 2, function(s) cumsum(s)))
  y = sapply(1:nrow(d$X), function(i) {
    if (i > 1)  {
      return(i+1 + sum(d$X[i,] * colSums(matrix(d$X[1:(i-1), ], ncol=e$p))))
    } else {
      return(2)
    }})
  d$Y = matrix(y, ncol=1)
  out = run.onlineAlgorithm(d, e, algorithm=asgd.onlineAlgorithm)
  out.sgd = run.onlineAlgorithm(d, e, algorithm=sgd.onlineAlgorithm)
  rand.t = sample(1:n, 1)
  theta.t = onlineOutput.estimate(out, rand.t)
  should.be = rowMeans(out.sgd$estimates[, 1:rand.t])
  CHECK_TRUE(all(abs(theta.t - should.be) < 1e-3))
}

test.implicit <- function() {
  p = 3
  n = 10
  e = normal.experiment(niters=n, p=p)
  alpha = 10 * runif(1, min=0, max=1)
  e$learning.rate <- function(t) {
    alpha
  }
  # Create dataset such that SGD (for the specified learning rate)
  # will give estimates that will satisfy:
  #  theta_t = x1 + x2 + x3 + ....x_t
  d = list()
  d$X = matrix(1, nrow=n, ncol=p, byrow=T)
  y = rep(1/alpha, n)
  d$Y = matrix(y, ncol=1)
  
  out = run.onlineAlgorithm(d, e, algorithm=implicit.onlineAlgorithm)
  
  U = matrix(1, nrow=p, ncol=p)
  B = (diag(p) - (alpha / (1+alpha * p)) * U)  # inverse of I + aU
  # print(B %*% (diag(p) + alpha* U))
  matrix.pow <- function(k) {
    # calculates B + B^2 + B^3 + ...B^k
    # Observe  Sk = B * (I + Sk-1)
    if(k == 0) return(matrix(0, nrow=p, ncol=p))
    return(B %*% (diag(p) + matrix.pow(k - 1)))
  }
  
  rand.t = sample(1:n, 1)
  # recall for the normal model
  #   θ_t = (I + a_t * Xt)^-1  * (θ_t-1 + at yt xt)
  # Here we assume at = a and yt = 1/a so at yt = 1
  # Therefore
  #   θ_t = (I + a Xt)^-1  * (θ_t-1 + xt)
  # Also xt = (1 1 1 ..) = u and so Xt = u u' = U = all ones
  # Also (I + aU)^-1 = I - a/(1 + ap) * U = B (by definition)
  # Therefore:
  #   θ_t = B * (θ_t-1 + u)
  # The solution is finally:
  #   θ_t = (B + B^2 + ...B^t) u = S(t) * u
  # Then
  # B * (θ_t-1 + u) = B * (S(t-1) * u + u) = (B + B S(t-1)) * u = S(t) * u
  theta.t = onlineOutput.estimate(out, rand.t)
  should.be = rowSums(matrix.pow(rand.t))
  CHECK_TRUE(all(abs(theta.t - should.be) < 1e-2))
}



