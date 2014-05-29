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
library(mvtnorm)
library(biglm)
rm(list=ls())

create.dataset <- function(dim.p=10**2, dim.n=10**5) {
  print(sprintf("> Creating dataset with #obs=%d and #covariates=%d", 
                dim.n, dim.p))
  sample.beta <- function() {
    # Samples a beta vector (px1) of model parameters.
    # The model parameters are either "high", "medium" or "low" (levels)
    beta.levels = c(-1, -0.35, -0.001, 0.001, 0.35, 1)
    # How frequent are the levels
    beta.level.freq = c(0.05, 0.1, 0.25, 0.3, 0.2, 0.1)    
    beta = c()
    for(i in 1:length(beta.levels)) {
      beta <- c(beta, 
                rep(beta.levels[i], as.integer(dim.p * beta.level.freq[i])))
    }
    if(length(beta) != dim.p)
      stop("Beta vector of parameters has wrong size.")
    return(sample(beta))
  }
  
  # Samples the covariates
  sample.X <- function() {
    X = matrix(rbinom(dim.n * dim.p, size=1, prob=0.1), nrow=dim.n)
  }
  
  X = sample.X()
  beta = sample.beta()
  Y = X %*% beta + rnorm(dim.n, sd=1.0)
  dataset = list(X=X, beta=beta, Y=Y)
  save(dataset, file="dataset.Rdata")
  print("> File saved as dataset.Rdata")
}

vector.dist <- function(x, y) {
  if(length(x) != length(y))
    stop("Vectors should have equal length to calculate distance.")
  # Compute MSE.
  sqrt(sum((x-y)^2))
}
CHECK <- function(claim, msg) {
  if(!claim) 
    stop(msg)
}

analyze.dataset <- function(method="lm") {
  # Fits the model using the dataset and prints performance metrics.
  # method = { lm, sgd }
  #
  print("> Loading dataset..")
  load(file="dataset.Rdata")
  print(sprintf("> Dataset loaded. no.obs=%d no.covariates=%d",
                nrow(dataset$X), ncol(dataset$X)))
  true.beta = dataset$beta
  print(sprintf("> Fitting using %s()..", method))
  ## Time performance of method
  perf = 0
  ## beta estimates from the method
  beta.hat = rep(0, ncol(dataset$X))
  dataset$beta = NULL # remove ground truth
  if(method=="lm") {
    perf = system.time({ beta.hat = analyze.dataset.lm(dataset) })
  } else if(method=="sgd"){
    perf = system.time({ beta.hat = analyze.dataset.sgd(dataset) })
  } else if(method=="biglm") {
    perf = system.time({ beta.hat = analyze.dataset.biglm(dataset) })
  } else {
    stop(sprintf("Method  %s not supported..", method))
  }
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

analyze.dataset.lm <- function(dataset) {
  fit = lm(Y ~ 0 + X, data=dataset)
  return(as.numeric(fit$coefficients))
}

analyze.dataset.biglm <- function(dataset) {
  dat = as.data.frame(dataset)
  data(dat)
  formula = Y ~ 0 + .
  print(names(dat))
  fit = bigglm(formula, data=dat)
  return(as.numeric(fit$coefficients))
}

analyze.dataset.sgd <- function(dataset) {
  p = ncol(dataset$X)
  n = nrow(dataset$X)
  beta.old = rep(0, p) ## initial estimate
  beta.new = beta.old
  print(sprintf("> Setting learning rates."))
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
  a.optimal = 1 / (sum(dataset$X) / (p * n))
  print(sprintf("Optimal a = %.3f", a.optimal))
  learning.rates = sapply(1:n, function(i) 1 / (1 + (1/a.optimal) * i))
  
  # 2. Set the SGD method (either explicit or implicit)
  Y = dataset$Y
  X = dataset$X
  use.explicit = F # what method to use
  print(sprintf("> Using %s updates in SGD",
                ifelse(use.explicit, "explicit", "implicit")))
  pb = txtProgressBar(style=3)
  
  for(i in 1:n) {
    ai = learning.rates[i]
    xi = X[i, ] # current covariate vector
    Yi.pred = sum(beta.old * xi)
    if(use.explicit) {
      beta.new = beta.old + ai * (Y[i] - Yi.pred) * xi      
    } else {
      ksi = ai * (Y[i] - Yi.pred) / (1 + sum(xi^2) * ai) 
      beta.new = beta.old + ksi * xi # implicit sgd  
    }
    beta.old = beta.new
    setTxtProgressBar(pb, value=i/n)
  }
  return(beta.new)
}
