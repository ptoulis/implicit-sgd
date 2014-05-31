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
library(Matrix)
library(mvtnorm)
library(biglm)
rm(list=ls())

dataset.filename <- function(dim.p, dim.n, is.sparse) {
  return(sprintf("datasets/%sDataset-p%d-n%d.Rdata", 
                 ifelse(is.sparse,"sparse", ""),
                 log(dim.p, 10), log(dim.n, 10)))
}

create.dataset <- function(dim.p=10**2, dim.n=10**5, sparse=F) {
  print(sprintf("> Creating dataset with #obs=%d and #covariates=%d", 
                dim.n, dim.p))
  sample.beta <- function() {
    # 1. Samples a beta vector (px1) of model parameters.
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
    if(sparse) {
      return(as(sample(beta), "sparseMatrix"))
    } else {
      return(sample(beta))
    }
  }
  
  # 2. Sample the covariates (nxp) sparse matrix.
  sample.X <- function() {
    # Create a sparse binary feature matrix
    L = dim.p * dim.n 
    if(sparse) {
      # 1. Pick non-zeros (and positions in the matrix)
      nnzeros = min(0.05 * L, L**0.7)
      print(sprintf("> Sampling %e non-zeros out of %e observations.", nnzeros, L))
      nnzero.pos = sample(L, size=nnzeros, replace=F)
      print("> Sampling finished. Sampling dataset.")
      nnzero.i = 1 + as.integer(nnzero.pos / dim.p)
      nnzero.j = nnzero.pos - (nnzero.i-1) * dim.p
      k =  which(nnzero.j==0)
      nnzero.i[k] <- nnzero.i[k] - 1
      nnzero.j[k] <- dim.p
      return(sparseMatrix(i=nnzero.i, j=nnzero.j, dims=c(dim.n, dim.p)))
    } else {
      # 1. Pick non-zeros (and positions in the matrix)
      nnzeros = min(0.05 * L)
      print(sprintf("> Sampling %e non-zeros out of %e observations.", nnzeros, L))
      nnzero.pos = sample(L, size=nnzeros, replace=F)
      print("> Sampling finished. Sampling dataset.")
      
      X = matrix(0, nrow=dim.n, ncol=dim.p)
      X[nnzero.pos] <- 1
      return(X)
    }
  }
  
  X = sample.X()
  beta = sample.beta()
  # 3. Sample the outcomes (nx1)
  Y = X %*% beta + rnorm(dim.n, sd=1.0)
  print("> Dataset size...")
  print(object.size(X), units="Mb")
  filename = dataset.filename(dim.p, dim.n, sparse)
  dataset = list(X=X, beta=beta, Y=Y)
  print(sprintf("> Saving file to %s", filename))
  save(dataset, file=filename)
  rm(list=c("beta","X", "Y")) # free up some memory
  gc()
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

analyze.dataset <- function(dim.p, dim.n, method="lm", sparse=F) {
  # Fits the model using the dataset and prints performance metrics.
  # method = { lm, sgd }
  #
  filename = dataset.filename(dim.p, dim.n, is.sparse=sparse)
  print(sprintf("> Loading dataset..%s. Is it sparse? %s", filename, sparse))
  load(file=filename)
  print(sprintf("> Dataset loaded. no.obs=%d no.covariates=%d",
                nrow(dataset$X), ncol(dataset$X)))
  true.beta = as.numeric(dataset$beta)
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

analyze.dataset.sparse.lm <- function(dataset) {
  
}

analyze.dataset.biglm <- function(dataset) {
  df = as.data.frame(dataset)
  fit = bigglm(terms(Y ~ 0 + ., data=df), data=df)
  return(as.numeric(coef(fit)))
}

analyze.dataset.sgd <- function(dataset) {
  p = ncol(dataset$X)
  n = nrow(dataset$X)
  Y = dataset$Y
  X = dataset$X
  
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
  a.optimal = 1 / (sum(X) / (n * p))
  print(sprintf("Optimal a = %.3f", a.optimal))
  learning.rates = sapply(1:n, function(i) 1 / (1 + (1/a.optimal) * i))
  
  # 2. Set the SGD method (either explicit or implicit)
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

fit.implicit <- function(x, y) {
  
}