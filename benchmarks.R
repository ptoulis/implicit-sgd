# Copyright (c) 2013
# Panos Toulis, ptoulis@fas.harvard.edu
#
# Benchmarks. Contains code to make benchmark comparisons
# across algorithms. Defines evaluation metrics and creates plots.
#
source("online-algorithms.R")

bias.benchmark <- function() {
  algos = c("sgd.onlineAlgorithm", "implicit.onlineAlgorithm")
  e = normal.experiment(niters=5000, p=5)
  nsamples = 20
  # of the form LIST[[algorithm]][[no.sample]]
  # e.g. LIST[[sgd.onlineAlgorithm]][[3]] is the OnlineOutput object
  # of the 3rd run
  out.algo.list = run.online.algorithm.many(e, algos, nsamples=nsamples)
  bias=list()
  # Update the bias data.
  for(algoName in algos) {
    bias[[algoName]] = list(min=rep(0, e$niters), max=rep(0, e$niters))
    out.list = out.algo.list[[algoName]]
    get.bias.samples <- function() {
      bias.matrix = matrix(0, nrow=nsamples, ncol=e$niters)
      for(i in 1:nsamples) {
        out = out.list[[i]]
        bias.matrix[i,] <- onlineOutput.bias(out, experiment=e)
      }
      return(bias.matrix)
    }
    
    B = get.bias.samples()
    bias[[algoName]]$min = apply(B, 2, function(x) quantile(x, probs=c(0.1)))
    bias[[algoName]]$max = apply(B, 2, function(x) quantile(x, probs=c(0.9)))
  }
  
  cols = topo.colors(length(algos))
  ## Plotting.
  for(i in 1:length(algos)) {
    algoName = algos[i]
    x = 1:e$niters
    ymax = log(bias[[algoName]]$max)
    ymin = log(bias[[algoName]]$min)
    m = min(ymin)
    M = max(ymax)
    if(i==1) {
      plot(x, ymax, main="Bias", xlab="iterations", type="l", col="white",
           ylim=c(m,M))
      legend(0.5 * e$niters, max(ymax), col=cols, legend=algos, lty=c(1, 2))
    }
    polygon(c(x, rev(x)), c(ymin, rev(ymax)), col=cols[i], lty=i)
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




