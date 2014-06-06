# Copyright (c) 2013
# Panos Toulis, ptoulis@fas.harvard.edu
#
# Using the simulation setup as in glmnet JoSS paper(Friedman, Hastie, Tibshirani)
#
rm(list=ls())
library(mvtnorm)
library(glmnet)

genjerry = function(x,snr){
  # generate data according to Friedman's setup
  n=nrow(x)
  p=ncol(x)
  b=((-1)^(1:p))*exp(-2*( (1:p)-1)/20)
  f=x%*%b
  e=rnorm(n)
  k=sqrt(var(f)/(snr*var(e)))
  y=f+k*e
  return(list(y=y, beta=b))
}

genx2 = function(n,p,rho){
  #    generate x's multivariate normal with equal corr rho
  if(abs(rho)<1){
    beta=sqrt(rho/(1-rho))
    x0=matrix(rnorm(n*p),ncol=p)
    z=rnorm(n)
    A = matrix(z, nrow=n, ncol=p, byrow=F)
    x= beta * A + x0
  }
  if(abs(rho)==1){ x=matrix(z,nrow=n,ncol=p,byrow=F)}
  
  return(x)
}

sample.data <- function(dim.n, dim.p, rho, 
                        model="gaussian") {
  # Samples the covariates as normal with the specific correlation
  X = genx2(dim.n, dim.p, rho)
  d = genjerry(X, 3)
  beta = d$beta
  Y = d$y
  if(model=="logit") {
    Y = as.vector(1* (runif(dim.n) < 1/(1 + exp(-Y))))
  }
  return(list(Y=Y, X=X, beta=beta))
}

optimal.alpha <- function(p, rho) {
  lambdas = c(1 + (p-1) * rho, rep(1-rho, p-1))
  f = function(x) sum(x^2 * lambdas / (2 * x * lambdas - 1))
  optim(par=c(0), fn=f, method="L-BFGS-B", lower=1/(2 * min(lambdas)) + 1e-3)$par
}

implicit.sgd <- function(x, y, rho, model="gaussian") {
  n = nrow(x)
  p = ncol(x)
  beta.hat = rep(0, p)
  a.optimal = optimal.alpha(p, rho)  # not sensitive to that
  # a.optimal = 1
  learning.rate = 1 / (1 + (1/a.optimal) * seq(1, n))
  for(i in 1:n) {
    xi = x[i, ]
    yi = y[i]
    yi.pred = sum(beta.hat * xi)
    ai = learning.rate[i]
    ksi = ai * (yi - yi.pred) / (1 + sum(xi^2) * ai)
    beta.hat = beta.hat + ksi * xi
  }
  return(beta.hat)
}

mse <- function(x, y) {
  if(length(x) != length(y))
    stop("MSE should compare vectors of same length")
  sqrt(mean((x-y)^2))
}
time.gaussian <- function(dim.n, dim.p,
                          rho.values=c(0.0, 0.1, 0.2),
                          nreps=3, 
                          methods=c("glmnet")) {
  niters = 0
  cols = c("rho", "rep", "time","mse", "method")
  timings = matrix(nrow=0, ncol=length(cols))
  colnames(timings) <- cols
  rownames(timings) = NULL
  total.iters = nreps * length(rho.values) * length(methods)
  pb = txtProgressBar(style=3)
  seeds=sample(1:1000000, size=total.iters)
  
  for(i in 1:nreps) {
    for(rho in rho.values) {
      for(method in methods) {
        niters = niters + 1
        set.seed(seeds[niters])
        dataset = sample.data(dim.n, dim.p, rho)
        true.beta = dataset$beta
        x = dataset$X
        y = dataset$Y
        new.dt = 0
        beta.est = NULL
        if(method=="glmnet") {
          new.dt = unix.time({ fit = glmnet(x,y, alpha=1, standardize=FALSE, type.gaussian="naive")})["user.self"]
          B = fit$beta
          best.j = which.min(sapply(1:ncol(B), function(j) {
            x1 = as.numeric(B[, j])
            mse(x1, true.beta) }))
          beta.est = B[, best.j]
        } else {
          new.dt = unix.time({ beta.est= implicit.sgd(x, y, rho)})["user.self"]
        }
        timings = rbind(timings, c(rho, i, round(new.dt, 5), 
                                   round(mse(beta.est, true.beta), 3), 
                                   method))
        setTxtProgressBar(pb, niters/total.iters)
      }
    }
    
  }
  return(as.data.frame(timings))
}

glmnet.wrapper=function(sx,ly,folds,lambda=lambda,standardize=FALSE,family="binomial"){
  nfolds=max(folds)
  for(i in 1:nfolds){
    fit2<-glmnet(sx[folds!=i,],ly[folds!=i],lambda=lambda,standardize=standardize,family=family)
  }
  return()
}

run.logit <- function(dim.n, dim.p,
                      rho.values=c(0.0, 0.1, 0.2),
                      nreps=3) {
  timings = matrix(nrow=0, ncol=2)
  colnames(timings) <- c("rho", "time")
  nlam=100
  nfolds=10
  pb = txtProgressBar(style=3)
  niters = 0
  for(rho in rho.values) {
    set.seed(22)
    niters = niters + 1
    dataset = sample.data(dim.n, dim.p, rho, model="logit")
    x = dataset$X
    y = dataset$Y
    # make sure n is a multiple of nfolds=10
    folds=sample(rep(1:nfolds, dim.n/nfolds))
    fit2=glmnet(x,y,standardize=FALSE,family="binomial")
    tim.fht=system.time(glmnet.wrapper(x,y,folds=folds,lam=fit2$lambda,standardize=FALSE,family="binomial"))[1]
    timings = rbind(timings, c(rho, as.numeric(tim.fht)))
    setTxtProgressBar(pb, niters/length(rho.values))
  }
  return(as.data.frame(timings))
}
