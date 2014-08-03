## Online EM algorithm
## Example 1. Poisson mixture (this is probably a bad example)
source("../../../r-toolkit/checks.R")

true.w = c(0.35, 0.65)
true.lambda = c(1, 5.5)
K = length(true.w) # no. of mixtures.

sample.y <- function(n) {
  I = sample(1:K, size=n, replace=T, prob=true.w)
  y = rep(0, n)
  for(i in 1:K) {
    y = y + (I==i) * rpois(n, lambda=true.lambda[i])
  }
  # poisson mixture.
  return(y)
}

w.ave <- function(y, w, lambdas) {
  w.conditional = w * lambdas^y * exp(-lambdas)
  N = sum(w.conditional)
  if(N == 0)
    return(rep(1, K)/K)
  
  w.bar = w.conditional / N
}


log.likelihood <- function(Y, w, lambdas) {
  el = 0
  for(y in Y) {
    w.bar = w.ave(y, w, lambdas)
    el = el + sum(w.bar * (log(w) - lambdas)) + sum(log(lambdas) * w.bar) * y
  }
  return(el)
}

run.online.em <- function(n) {
  Y = sample.y(n) 
  s.hat <- c(0, 2 * K)
  w.hat <- rep(0, K)
  lambda.hat <- rep(0, K)
  
  for(iter in 1:n) {
    ai = 1 / (1 + iter)
    yn = Y[iter]
    sn.expected = c(w.ave(yn, w.hat, lambda.hat), 
                    yn * w.ave(yn, w.hat, lambda.hat))
    #print(sn.expected)
    #print(sprintf("y = %.3f", yn))
    s.hat = s.hat + ai * (sn.expected - s.hat)
    w.hat = head(s.hat, K)
    lambda.hat = tail(s.hat, K) / w.hat
    
    # print(w.hat)
    # print(lambda.hat)
  }
  print(sprintf("w.hat = %s", paste(round(w.hat, 4), collapse=", ")))
  print(sprintf("lambda.hat = %s", paste(round(lambda.hat, 2), collapse=", ")))
  print(sprintf("Log-likelihood = %.3f", log.likelihood(Y, w.hat, lambdas = lambda.hat)))
  print(sprintf("Log-likelihood = %.3f", log.likelihood(Y, true.w, true.lambda)))
  
}


# Tests and checks.
test.w.ave <- function() {
  CHECK_NEAR(w.ave(1, rep(1, K), rep(1, K)), rep(1, K) / K, msg="by symmetry")
}
