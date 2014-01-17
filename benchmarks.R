# Copyright (c) 2013
# Panos Toulis, ptoulis@fas.harvard.edu
#
# Benchmarks. Contains code to make benchmark comparisons
# across algorithms. Defines evaluation metrics and creates plots.
#
source("online-algorithms.R")
library(scales)

vector.dist <- function(x1, x2) {
  m1 = matrix(x1, ncol=1)
  m2 = matrix(x2, ncol=1)
  return(matrix.dist(m1, m2))
}

matrix.dist <- function(m1, m2) {
  CHECK_EQ(nrow(m1), nrow(m1))
  CHECK_EQ(ncol(m1), ncol(m1))
  norm(m1-m2, "F") / sqrt(nrow(m1) * ncol(m1))
}

plot.low.high <- function(data, experiment) {
  # Given a LIST{algoName} it will assume fields (low, high)
  # and then plot the error lines.
  #
  algos = setdiff(names(data), "x")
  niters = experiment$niters
  cols = topo.colors(length(algos))
  ## Plotting.
  for(i in 1:length(algos)) {
    algoName = algos[i]
    x = data$x
    ymax = log(data[[algoName]]$low)
    ymin = log(data[[algoName]]$high)
    m = min(ymin)
    M = max(ymax)
    if(i==1) {
      plot(x, ymax, main="Bias", xlab="iterations", type="l", col="white",
           ylim=c(min(-3, m), max(3, M)))
      legend(0.5 * niters, max(ymax), col=cols, legend=algos, lty=c(1, 2))
    }
    polygon(c(x, rev(x)), c(ymin, rev(ymax)), col=alpha(cols[i], 0.4), lty=i)
  }
}

bias.benchmark.asymptotics <- function() {
  algos = c("sgd.onlineAlgorithm", "implicit.onlineAlgorithm")
  e = normal.experiment(niters=1000, p=100)
  niters = e$niters
  e$learning.rate <- function(t) base.learning.rate(t, gamma0=1/sum(diag(e$Vx)),
                                                    alpha=0.1, c=1)
  nsamples = 10
  # 1. Run the algorithms. Get a MultipleOnlineOutput object.
  mul.out = run.online.algorithm.many(e, algos, nsamples=nsamples)
  # 2. Process to get the bias
  bias=list(x=1:niters)
  for(algoName in algos) {
    theta.t.fn <- function(theta) vector.dist(theta, e$theta.star)
    summary.min = function(x) quantile(x, 0.05)
    summary.max = function(x) quantile(x, 0.95)
    bias[[algoName]]$low = mul.OnlineOutput.vapply(e, mul.out, algoName,
                                                   theta.t.fn, summary.min)
    bias[[algoName]]$high = mul.OnlineOutput.vapply(e, mul.out, algoName,
                                                   theta.t.fn, summary.max)
  }
  # 3. Plot low/high curves.
  plot.low.high(bias, e)
}

bias.benchmark.learningRate <- function() {
  algos = c("sgd.onlineAlgorithm", "implicit.onlineAlgorithm")
  e = normal.experiment(niters=2000, p=5)
  nsamples = 10
  alpha.values = seq(0.01, 5.0, length.out=20)
  bias=list(x=alpha.values)
  
  for(alpha.index in 1:length(alpha.values)) {
    alpha = alpha.values[alpha.index]
    print(sprintf("Iteration %d/%d: alpha=%.2f", alpha.index, length(alpha.values), 
                  alpha))
    e$learning.rate <- function(t) alpha / (t+1)
    # 1. Run the algorithms. Get a MultipleOnlineOutput object.
    mul.out = run.online.algorithm.many(e, algos, nsamples=nsamples)
    # 2. Process to get the bias
    for(algoName in algos) {
      theta.t.fn <- function(theta) vector.dist(theta, e$theta.star)
      summary.min = function(x) quantile(x, 0.05)
      summary.max = function(x) quantile(x, 0.95)
      bias[[algoName]]$low[alpha.index] = min(mul.OnlineOutput.vapply(e, mul.out, algoName,
                                                                      theta.t.fn, summary.min))
      
      bias[[algoName]]$high[alpha.index] = max(mul.OnlineOutput.vapply(e, mul.out, algoName,
                                                      theta.t.fn, summary.max))
    }
    bias$x = alpha.values[1:alpha.index]
    plot.low.high(bias, e)
  }
  # 3. Plot low/high curves.
  plot.low.high(bias, e)
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




