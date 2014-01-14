# Copyright (c) 2013
# Panos Toulis, ptoulis@fas.harvard.edu
#
# Benchmarks. Contains code to make benchmark comparisons
# across algorithms. Defines evaluation metrics and creates plots.
#
source("online-algorithms.R")

experiment.benchmark <- function(experiment.name, niters=10000) {
  e = get.experiment(name="normal", niters=niters)
  d = e$sample.dataset()
  algos =  c("sgd", "asgd", "implicit")
  cols = c("red", "black", "green")
  
  for(i in 1:length(algos)) {
    out = NA
    if(i==1) {
      out = run.onlineAlgorithm(dataset=d, experiment=e, algorithm=sgd.onlineAlgorithm)
    } else if(i==2) {
      out = run.onlineAlgorithm(dataset=d, experiment=e, algorithm=asgd.onlineAlgorithm)
    } else {
      out = run.onlineAlgorithm(dataset=d, experiment=e, algorithm=implicit.onlineAlgorithm)
    }
    x = log(1:ncol(out$estimates), base=10)
    y = log(apply(out$estimates, 2, e$risk), base=10)
    if(i==1) {
      plot(x, y, type="l", col=cols[i])
      legend(0.5, -1, legend=c("sgd", "asgd", "implicit"), col=cols, lty=rep(1, 3))
    } else {
      lines(x, y, col=cols[i])
    }
  }
}

variance.benchmark <- function(nsamples, nT) {
  niters = nT
  p = 4
  
  nrepeats = nsamples
  theta.T = matrix(0, nrow=p, ncol=0)
  
  e = normal.experiment(niters=niters, p=p)
  alpha = e$learning.rate(100000) * 100000
  print(sprintf("Alpha=%.2f", alpha))
  I = diag(p)
  A = e$Vx
  print("A=")
  print(A)

  print("S theoretical=")
  C.theoretical = alpha^2 * solve(2 * alpha * A - I) %*% A
  print(C.theoretical)
  print(eigen(C.theoretical)$values)
  CHECK_TRUE(all(eigen(C.theoretical)$values > 0), msg="Should be > 0")
  for(i in 1:nrepeats) {
    d = e$sample.dataset()
    out = run.onlineAlgorithm(d, e, sgd.onlineAlgorithm)
    theta.T <- cbind(theta.T, out$last)
    print("Last estimate=")
    print(out$last)
    C.empirical = nT * cov(t(theta.T))
    if(i %% 5==0) {
      print(sprintf("Out of %d samples : Empirical covariance=", i))
      print(C.empirical)
      print("theoretical=")
      print(C.theoretical)
      print(mean((C.empirical - C.theoretical)^2))
    }
  }
}




