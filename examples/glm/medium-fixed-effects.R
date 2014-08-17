# Copyright (c) 2013
# Panos Toulis, ptoulis@fas.harvard.edu
#
## This is a self-contained file as it needs to work without reference to the 
## sgd R library here.
## 
## Goal: Create 10 million observations (Y) and 10k binary covariates (X)
##       and fit a linear normal model Y = X b + e
##
## Example usage:
#   create.dataset()  # creates a dataset with 100 covariates, 100,000 obs.
#   analyze.dataset(method="sgd") # uses sgd (implicit) to fit
#   analyze.dataset(method="lm") # uses R's lm() method.
#
#  Note: create.dataset() saves the dataset in a file (current wd)
#
rm(list=ls())
library(Matrix)
source("datasets.R")



analyze.dataset <- function(dim.p, dim.n, method="glm", verbose=F) {
  # Fits the model using the dataset and prints performance metrics.
  # method = { lm, sgd }
  if(dim.p > 1e3 || dim.n > 1e5 || dim.p * dim.n > 1e8)
    stop("This is a Big dataset. Please use the large-fixed-effects.R file to analyze it.")
  # 1. Load dataset
  df = load.dataset(dim.p, dim.n)
  CHECK(nrow(df)==dim.n, msg="Correct #rows")
  CHECK(ncol(df)==dim.p + 1, msg="Correct #cols")
  # 2. Load truth
  true.beta = load.beta(dim.p, dim.n)
  print(sprintf("> Fitting using %s()..", method))
  ## Time performance of method
  perf = 0
  ## beta estimates from the method
  beta.hat = rep(0, dim.p)
  if(method=="glm") {
    perf = system.time({ beta.hat = analyze.dataset.glm(df) })
  } else if(method=="sgd"){
    perf = system.time({ beta.hat = analyze.dataset.sgd(df) })
  } else {
    stop(sprintf("Method  %s not supported..", method))
  }
  if(verbose) {
    print(sprintf("> %s() took %.2f secs. MSE=%.3f", 
                  method,
                  perf[["elapsed"]],
                  vector.dist(true.beta, beta.hat)))
    plot(true.beta, beta.hat, ylim=c(min(true.beta), max(true.beta)), pch=20)
    points(true.beta, true.beta, pch="x", col="magenta")
    for(b in true.beta) {
      abline(h=b, lty=3, lwd=0.3)
    }
  }
  return(list(time=perf[["elapsed"]], mse=vector.dist(beta.hat, true.beta)))
}

analyze.dataset.glm <- function(dataset) {
  fit = lm(Y ~ 0 + ., data=dataset)
  return(as.numeric(fit$coefficients))
}


solve.best.alpha <- function(q.hat, p) {
  # Solves for the best possile alpha rate based on (Toulis et al, 2014)
  # Recall that the variance is proportional to (2aJ-I)^-1 J (see below)
  # Therefore, imposing a structure on the Fisher information matrix
  # (e.g. consider the normal linear model), then we can figure out 
  # the learning rate that minimizes some metric.
  # (in this case we minimize the trace).
  # 
  # The computation assumes that all columns but the first of the design X matrix
  # has independent probability of being 1, equal to q.
  # q.hat = estimate of this probability based on "1" in X.
  #
  # q.hat = as.numeric((rowSums(S.hat)[1] + colSums(S.hat)[1] + sum(diag(S.hat)) - 3) / (3 * (p-1)))
  # print(q.hat)
  # This solves the two eigenvalues that are not equal to q (1-q) with mul (p-2)
  # (recall p =#covariates)#
  #
  lambdas = as.numeric(polyroot(c(q.hat * (1-q.hat), -(1+q.hat + q.hat^2* (p-2)), 1)))
  other.lambdas = rep(q.hat * (1-q.hat), p-2)
  lambdas = c(lambdas, other.lambdas)
  out = optim(par=0, 
              fn=function(x) sum(x^2 * lambdas / (2 * x * lambdas-1)), 
              lower=1 / (2 * min(lambdas)) + 1e-4,
              method="L-BFGS-B")
  return(out$par)
}

best.alpha <- function(J) {
  lambdas = eigen(J)$values
  f = function(x) sum(x^2 * lambdas / (2 * x * lambdas - 1))
  m = 1 / (2 * min(lambdas))
  m = ifelse(is.finite(m), m, 0)
  optim(par=c(0), fn=f, lower= m + 1e-4, upper=100, method="L-BFGS-B")$par
}

analyze.dataset.sgd <- function(dataset) {
  p = ncol(dataset) - 1
  n = nrow(dataset)
  Y = as.numeric(dataset[, p+1])
  X = as.matrix(dataset[, 1:p])
  # X = dataset[, 1:p]  # this will be very slow. Random access in data-frames is very slow.

  beta.old = rep(0, p) ## initial estimate
  beta.new = beta.old
  # print(sprintf("> Setting learning rates."))
  # 1. Pick the optimal learning rate:
  # By (Toulis et.al., 2014) the variance will be
  # V = a^2 f^2 (2af J - I)^-1 * J
  #  where a=learning rate, f=dispersion param=Var(y), J=fisher information, I=identity matrix.
  #  But in this model J = (1/f) E(xx') = (1/f) q I  where q=P(xi=1)
  #  and so, f * J = q I.  Thus the variance is V = f * q * (a^2) / (2aq-1) I
  # To minimize variance we thus need to minimize a^2 / (2aq-1)
  # which leads to a.opt = 1/q. If we assume a_n = 1 / (1 + h * n) then
  # clearly a_n * n -> a = 1/h and thus we need to set h=q=1/a for the optimal value.
  #
  #a.optimal = 1 / (sum(X) / (n * p))
  #learning.rates = 1 / (1 + (1/a.optimal) * seq(1, n))
  
  # 2. Set the SGD method (either explicit or implicit)
  use.explicit = F # what method to use
  print(sprintf("> Using %s updates in SGD.",
                ifelse(use.explicit, "explicit", "implicit")))
  pb = txtProgressBar(style=3)
  # J.est <- diag(p)
  q.hat = (sum(X) - n) / (n * (p-1))
  alpha.opt = solve.best.alpha(q.hat, p)
  print(sprintf("q.hat=%.3f -- Learning rate limit a = %.3f", q.hat, alpha.opt))
  for(i in 1:n) {
    xi = X[i, ] # current covariate vector
    yi = Y[i]
    ai = 1 / (1 + (1/alpha.opt) * i)
    Yi.pred = sum(beta.old * xi)
    Si = sum(xi^2)
    if(use.explicit) {
      beta.new = beta.old + ai * (yi - Yi.pred) * xi      
    } else {
      ksi = ai * (yi - Yi.pred) / (1 + Si * ai) 
      beta.new = beta.old + ksi * xi # implicit sgd  
    }
    beta.old = beta.new
    setTxtProgressBar(pb, value=i/n)
  }
  return(beta.new)
}


## Main function to create the benchmark
run.small.experiment <- function(dim.p.vector=seq(10, 500),
                          dim.n.vector=seq(100, 50000),
                          methods=c("glm", "sgd"),
                          nsamples=5,
                          results.file="datasets/results/results.small.experiment.rda") {
  require(stringr)
  results = matrix(0, nrow=0, ncol=6)
  metrics = c("time.", "mse.")
  colnames(results) = c("p", "n", as.character(sapply(methods, function(m) str_c(metrics, m))))
  
  for(i in 1:nsamples) {
    n = sample(dim.n.vector, size=1)
    p = sample(dim.p.vector, size=1)
    print(sprintf(">> Iteration %d / %d --> Experiment (%d, %d)", 
                  i, nsamples, p, n))
    # filename = dataset.filename(dim.p=p, dim.n=n)
    create.dataset(dim.p=p, dim.n=n)
    out.row = c(p, n)
    for(m in 1:length(methods)) {
      method = methods[m]
      out = analyze.dataset(p, n, method=method)
      out.row = c(out.row, out$time, out$mse)
    }
    print(round(out.row, 3))
    results = rbind(results, out.row)
    rownames(results) <- NULL
    remove.dataset(p, n)
    ## Save results
    res.small.exp = as.data.frame(results)
    save(res.small.exp, file=results.file)
  }
  return(results)
}

results.small.experiment <- function() {
  load("datasets/results/results.small.experiment.rda")
  df = res.small.exp
  log.p = log(df$p)
  log.n = log(df$n)
  # 1. time GLM
  fit = lm(log(df$time.glm) ~ log.p + log.n)
  print("======   Time GLM() parameters (log p, log n)")
  print(summary(fit))
  
  # 3. mse GLM
  fit = lm(log(df$mse.glm) ~ log.p + log.n)
  print("======   MSE GLM() parameters (log p, log n)")
  print(summary(fit))
  
  # 2. time SGD
  fit = lm(log(df$time.sgd) ~ log.p + log.n)
  print("======   Time SGD() parameters (log p, log n)")
  print(summary(fit))
  
  
  # 4. mse SGD
  fit = lm(log(df$mse.sgd) ~ log.p + log.n)
  print("======   MSE SGD() parameters (log p, log n)")
  print(summary(fit))
  
  ## average MSE difference
  se = (df$mse.sgd - df$mse.glm) / df$mse.glm
  avg.se = mean(se)
  hist(se)
  se = sd(replicate(1000, { x = sample(se, size=length(se), replace=T); mean(x)}))
  print(sprintf("Difference MSE = %.2f (%.2f)", avg.se, 2 * se))
  
}






