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

rao.example <- function(niters) {
  y = c(125, 18, 20, 34)
  
  typical.em <- function(niters) {
    pi.list <- c(0)
    for(i in 1:niters) {
      # impute missing sufficient statistic (E step)
      pi.t <- tail(pi.list, 1)
      x.imp = c(125 * 2  / (2 +  pi.t), 
                125 * pi.t / (2 + pi.t),
                18, 20, 34)
      # M (step)
      pi.list <- c(pi.list, (x.imp[2] + 34) / (x.imp[2] + 34 + 18 + 20))
    }
    print("Typical EM")
    print(tail(pi.list, 10))
    plot(pi.list, type="l", ylim=c(0, 1))
  }
  
  C = matrix(c(1, -1, -1, 1,
               -1, 1, 1, -1,
               -1, 1, 1, -1,
               1, -1, -1, 1), nrow=4, ncol=4, byrow=T)
  # eigen=4 and 0
  pi.true = 0.62682
  st.em <- function(niters, alpha) {
    get.prob <- function(pi) {
      c(0.5 + 0.25 * pi, 
        0.25 * (1-pi), 
        0.25 * (1-pi), 
        0.25  * pi)
    }
    prob.true = get.prob(pi.true)
    N = sum(y)
    pi.est <- c(0)
    pi.multi = rep(1, 4) / 4
    get.single.pi <- function(pi.mul) {
      x = 4 * pi.mul - c(2, 1, 1, 0)
      mean(x * c(1, -1, -1, 1))
    }
      
    for(i in 1:niters) {
      pi.t = tail(pi.est, 1)
      ai = 1 / (1 + (1/a) * i)
      y.rep = rmultinom(1, size = length(y), prob = get.prob(pi.true))
      y.exp <- rmultinom(1, size = length(y), prob = pi.multi)
      
      pi.multi.next  <- pi.multi + ai * (y.rep - y.exp)
      if(all(pi.multi.next >= 0) & all(pi.multi.next <= 1))
        pi.multi = pi.multi.next
      
      pi.est <- c(pi.est, get.single.pi(pi.multi))
    }
    return(pi.est)
  }
  
  # 1. Run typical EM

  # typical.em(niters)
  alpha.values = seq(0.05, 0.5, length.out=10)
  for(a in alpha.values) {
    a.vals = c()
    for(j in 1:500) {
      a.vals = c(a.vals, tail(st.em(niters, a)))
    }
    print(sprintf("a=%.3f, Var=%.3f bias=%.3f MSE=%.3f", 
                  a, 1000 * var(a.vals), 
                  1000 * mean(a.vals - pi.true)^2,
                  1000 * (var(a.vals) + mean(a.vals-pi.true)^2)))
  }
}




# Tests and checks.
test.w.ave <- function() {
  CHECK_NEAR(w.ave(1, rep(1, K), rep(1, K)), rep(1, K) / K, msg="by symmetry")
}
