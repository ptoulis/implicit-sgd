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
  if(dim.p >= 1e3 || dim.n >= 1e6 || dim.p * dim.n > 1e8)
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

analyze.dataset.sgd <- function(dataset, tol=1e-5) {
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
  a.optimal = 1 / (sum(X) / (n * p))
  learning.rates = 1 / (1 + (1/a.optimal) * seq(1, n))
  
  # 2. Set the SGD method (either explicit or implicit)
  use.explicit = F # what method to use
  print(sprintf("> Using %s updates in SGD. Optimal a=%.3f",
                ifelse(use.explicit, "explicit", "implicit"), a.optimal))
  pb = txtProgressBar(style=3)
  
  for(i in 1:n) {
    ai = learning.rates[i]
    xi = X[i, ] # current covariate vector
    yi = Y[i]
    Yi.pred = sum(beta.old * xi)
    Si = sum(xi)
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
run.benchmark <- function(dim.p.vector=c(10, 100),
                          dim.n.vector=c(100, 1000, 10000),
                          methods=c("glm", "biglm", "sgd"),
                          ntrials=5) {
  avoid.method <- function(p, n, method) {
    c1 = (method=="glm" && n >= 10**5)
    c2 = (method=="biglm" && p >= 1000)
    return(c1 || c2)
  }
  benchmark.results = list(method=)
}







