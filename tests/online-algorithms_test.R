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

test.implicit <- function() {
  p = 3
  n = 10
  e = normal.experiment(niters=n, p=p)
  alpha = 0.01
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
  matrix.pow <- function(arg.M, arg.n) {
    if(arg.n == 0) return(diag(nrow(arg.M)))
    if(arg.n == 1) return(arg.M)
    return(arg.M %*% matrix.pow(arg.M, arg.n-1))
  }
  
  rand.t = sample(1:n, 1)
  theta.t = onlineOutput.estimate(out, rand.t)
  print(theta.t)
  # print(B)
  print(rand.t)
  print(dim(xsums))
  print(matrix.pow(B, rand.t-1) %*% xsums[, rand.t-1])
}



