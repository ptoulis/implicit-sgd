# Copyright (c) 2013
# Panos Toulis, ptoulis@fas.harvard.edu
#
# Using the simulation setup as in glmnet JoSS paper(Friedman, Hastie, Tibshirani)
# http://www.jstatsoft.org/v33/i01/
rm(list=ls())
source("datasets.R")
library(mvtnorm)
library(glmnet)
library(plyr)

# genjerry, genx2 are functions taken from the above paper.
#
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
  # Xi = b Z + Wi, and Z, Wi are independent normal.
  # Then Var(Xi) = b^2 + 1
  #  Cov(Xi, Xj) = b^2  and so cor(Xi, Xj) = b^2 / (1+b^2) = rho
  z=rnorm(n)
  if(abs(rho)<1){
    beta=sqrt(rho/(1-rho))
    x0=matrix(rnorm(n*p),ncol=p)
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
  if(model=="logistic") {
    Y = as.vector(1* (runif(dim.n) < 1/(1 + exp(-Y))))
  }
  return(list(Y=Y, X=X, beta=beta))
}

optimal.alpha <- function(p, rho) {
  b = sqrt(rho / (1-rho))
  lambdas = c(p * b^2 + 1, rep(1, p-1))
  f = function(x) sum(x^2 * lambdas / (2 * x * lambdas - 1))
  optim(par=c(0), fn=f, method="L-BFGS-B", lower=1/(2 * min(lambdas)) + 1e-4)$par
}

implicit.sgd <- function(x, y, rho, model="gaussian") {
  n = nrow(x)
  p = ncol(x)
  beta.hat = rep(0, p)
  a.optimal = 1 #optimal.alpha(p, rho)  # not sensitive to that
  # print(sprintf("Optimal a=%.3f  r=%.3f", a.optimal, rho))
  # Define transfer function.
  h <- function(eta) {
    if(model=="gaussian") {
      return(eta)
    } else if(model=="logistic") {
      if(eta > 20) {
        return(1)
      } else if (eta < -20) {
        return(0)
      } else {
        return(exp(eta) / (1 + exp(eta)))
      }
    } else {
      stop(print(sprintf("Model %s is not supported..", model)))
    }
  }
  
  learning.rate = 1 / (1 + (1/a.optimal) * seq(1, n))
  solve.implicit <- function(ai, yi, eta.i, xi.norm, Bi) {
    if(Bi[1]==Bi[2]) return(Bi[1])
    f = function(w) w - ai * (yi - h(eta.i + xi.norm * w))
    uniroot(f, interval=Bi, tol=1e-4)$root
  }
  
  for(i in 1:n) {
    xi = x[i, ]
    xi.norm = sum(xi^2)
    yi = y[i]
    ai = learning.rate[i]
    eta.i = sum(beta.hat * xi)
    yi.pred = h(eta.i)
    ri = ai * (yi - yi.pred)
    Bi = c(0, ri)
    if(ri <= 0)
      Bi = c(ri, 0)
    
    # Solve implicit equation
    ksi = NA
    if(model=="gaussian") {
      # For the Gaussia model this is easy.
      ksi = ri / (1 + ai * xi.norm)
    } else {
      # Call implicit solver
      ksi = solve.implicit(ai, yi, eta.i, xi.norm, Bi)
    }
    # Update beta vector.
    beta.hat = beta.hat + ksi * xi
  }
  return(beta.hat)
}

mse <- function(x, y) {
  if(length(x) != length(y))
    stop("MSE should compare vectors of same length")
  sqrt(mean((x-y)^2))
}

mse.glmnet <- function(fit, true.beta) {
  B = fit$beta
  mse.lambda = sapply(1:ncol(B), function(j) {
    x1 = as.numeric(B[, j])
    mse(x1, true.beta) })
  sort.mse.index = rev(order(mse.lambda))
  k = length(sort.mse.index)
  q1 = mse.lambda[sort.mse.index[as.integer(0.25 * k)]]
  q2 = mse.lambda[sort.mse.index[as.integer(0.5 * k)]]
  q3 = mse.lambda[sort.mse.index[as.integer(0.75 * k)]]
  return(c(q1, q2, q3))
}

analyze.timings <- function(timings, rho.values, methods) {
  times = as.data.frame(timings)
  for(m in 1:length(methods)) {
    print(sprintf("----> Results for method %s", methods[m]))
    df = ddply(subset(times, method==m), .(rho), function(z) { c(mean(z$time), 
                                                                 mean(z$mse1),
                                                                 mean(z$mse2),
                                                                 mean(z$mse3)) })
    print(paste(rho.values, collapse= " & "))
    print(sprintf("Times= %s", paste(round(df[,2], 3), collapse=" & ")))
    print(sprintf("MSE (Q1) = %s", paste(round(df[,3], 3), collapse=" & ")))
    print(sprintf("MSE (Q2) = %s", paste(round(df[,4], 3), collapse=" & ")))
    print(sprintf("MSE (Q3) = %s", paste(round(df[,5], 3), collapse=" & ")))
  }
}

run.experiment.glmnet <- function(dim.n, dim.p,
                                  model="gaussian",
                                  rho.values=c(0.0, 0.1, 0.2, 0.5, 0.9, 0.95),
                                  nreps=3, 
                                  methods=c("glmnet", "sgd")) {
  print(sprintf(">>> Running model %s...Methods= %s N=%d, p=%d", 
                model, paste(methods, collapse=","),
                dim.n, dim.p))
  niters = 0
  cols = c("rho", "rep", "time","mse1", "mse2", "mse3", "method")
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
        dataset = sample.data(dim.n=dim.n, dim.p=dim.p, rho=rho, model=model)
        true.beta = dataset$beta
        x = dataset$X
        y = dataset$Y
        true.beta = dataset$beta
        CHECK(nrow(x) == dim.n)
        CHECK(ncol(x) == dim.p)
        new.dt = 0
        new.mse = NA
        
        if(method=="glmnet") {
          if(model=="gaussian") {
            new.dt = system.time({ fit = glmnet(x, y, alpha=1, standardize=FALSE, type.gaussian="naive")})[1]
          } else if(model=="logistic") {
            new.dt = system.time({ fit = glmnet(x,y,standardize=FALSE,family="binomial") })[1]
          }
          new.mse = mse.glmnet(fit, true.beta)
        } else if (method=="sgd") {
          ## method==sgd
          new.dt = system.time({ beta.est= implicit.sgd(x, y, rho, model=model) })[1]
          mse.im = mse(beta.est, true.beta)
          new.mse = rep(mse.im, 3)
        } else {
          stop(sprintf(">> Method %s is not supported ...", method))
        }
        m = which(methods==method)
        timings = rbind(timings, c(rho, i, 
                                   new.dt, 
                                   new.mse, m))
        setTxtProgressBar(pb, niters/total.iters)
      }
    }
    
  }
  analyze.timings(timings, rho.values, methods)
  return(timings)
}



glmnet.wrapper=function(sx,ly,folds,lambda=lambda,standardize=FALSE,family="binomial"){
  nfolds=max(folds)
  for(i in 1:nfolds){
    fit2<-glmnet(sx[folds!=i,],ly[folds!=i],lambda=lambda,standardize=standardize,family=family)
  }
  return()
}