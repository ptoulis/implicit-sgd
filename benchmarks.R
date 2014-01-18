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

plot.low.high <- function(data, experiment, draw) {
  # Given a LIST{algoName} it will assume fields (low, high)
  # and then plot the error lines.
  #
  algos = names(data)  # algorithm names (character)
  niters = experiment$niters  # no. of iterations
  cols = topo.colors(length(algos))  # colors.
  x = draw$x  # x-axis
  logY = draw$logY  # T or F, whether to get the log() of the outcome.
  logX = draw$logX
  if(logX) {
    x = log(x)
  }
  # Draw parameters.
  title = draw$main
  xlab = draw$xlab
  ylab = draw$ylab
  ## Plotting.
  for(i in 1:length(algos)) {
    algoName = algos[i]
    ymin = data[[algoName]]$low
    ymax = data[[algoName]]$high
    if(logY) {
      ymin = log(ymin)
      ymax = log(ymax)
    }
    
    defaultYlimMin = ifelse(logY, -3, 10^-3)
    defaultYlimMax = ifelse(logY, 3, 10^3)
    ylims = c(min(defaultYlimMin, min(ymin)), min(defaultYlimMax, max(ymax)))
    
    if(i==1) {
      plot(x, ymax, main=title, 
           xlab=xlab,
           ylab=ylab,
           type="l",
           col="white",
           ylim=ylims)
      legend(0.6 * niters, 0.8* max(ylims), col=cols, legend=algos, lty=1:length(algos))
    }
    polygon(c(x, rev(x)), c(ymin, rev(ymax)), col=alpha(cols[i], 0.4), lty=i)
  }
}

generic.benchmark <- function(algos, experiment, nsamples,
                              process.params) {
  # Runs a generic benchmark. The idea is the following:
  #   1) Define the algorithms to be tested
  #   2) Define the experiment (sample.dataset, score function, learning rate)
  #   3) Get #nsamples for each experiment
  #   4) For each sample use process.params to process the output.
  #
  # Returns a DATA object (not defined in terminology)
  # A DATA object is DATA{algorithmName}{low/high} = (niters x 1) vector.
  # process.params has two fields (vapply, theta.fn):
  # 
  # Define theta_tj = j-sample of vector θt. Then:
  #  If vapply=TRUE then DATA{algo}{low} = V, 
  #     where V(t) = 5% quantile (fn(theta_t1), fn(theta_t2)...)
  #
  # If vapply=FALSE then DATA{algo}{low} = V, 
  #     where V(t) = 5% quantile of fn(theta_tj)  i.e. it is applied in the entire matrix.
  #
  # Returns an object, LIST{algoName}{low/high} = []  of size #niters.
  #
  niters = experiment$niters
  # 1. Run the algorithms. Get a MultipleOnlineOutput object.
  mul.out = run.online.algorithm.many(experiment, algos, nsamples=nsamples)
  # 2. Process to get the data
  data = list()
  summary.min = function(x) quantile(x, 0.05)
  summary.max = function(x) quantile(x, 0.95)
  for(algoName in algos) {
    if(process.params$vapply) {
      theta.t.fn <- process.params$theta.fn
      data[[algoName]]$low = mul.OnlineOutput.vapply(experiment, mul.out, algoName,
                                                     theta.t.fn, summary.min)
      data[[algoName]]$high = mul.OnlineOutput.vapply(experiment, mul.out, algoName,
                                                      theta.t.fn, summary.max)
    } else {
      theta.fn <- process.params$theta.fn
      data[[algoName]]$low = summary.min(mul.OnlineOutput.mapply(experiment, mul.out, algoName, theta.fn))
      data[[algoName]]$high = summary.max(mul.OnlineOutput.vapply(experiment, mul.out, algoName, theta.fn))
    }
  }
  return(data)
}

bias.benchmark.asymptotics <- function() {
  algos = c("sgd.onlineAlgorithm", "implicit.onlineAlgorithm")
  e = normal.experiment(niters=300, p=100)
  e$learning.rate <- function(t) 2 / (1+t)
  #base.learning.rate(t, gamma0=1/sum(diag(e$Vx)),      alpha=0.1, c=1)
  nsamples = 10
  dist = function(theta) vector.dist(theta, e$theta.star)
  process.params = list(vapply=T, theta.fn=dist)
  # 1. Run the algorithms. Get a MultipleOnlineOutput object.
  data = generic.benchmark(algos=algos, 
                           experiment=e,
                           nsamples=nsamples,
                           process.params=process.params)
  # 3. Plot low/high curves.
  draw = list(x=1:e$niters, logY=F, logX=F,
              main="Bias asymptotics", xlab="Iterations", ylab="|| bias ||")
  plot.low.high(data, e, draw)
}

bias.benchmark.learningRate <- function() {
  algos = c("sgd.onlineAlgorithm", "implicit.onlineAlgorithm")
  e = normal.experiment(niters=100, p=5)
  nsamples = 5
  alpha.values = seq(0.01, 5.0, length.out=5)
  dist = function(theta) vector.dist(theta, e$theta.star)
  process.params = list(vapply=T, theta.fn=dist)
  
  alpha.data = list()
   
  for(alpha.index in 1:length(alpha.values)) {
    alpha = alpha.values[alpha.index]
    print(sprintf("Iteration %d/%d: alpha=%.2f", alpha.index, length(alpha.values), 
                  alpha))
    e$learning.rate <- function(t) alpha / (t+1)
    # 1. Run the algorithms. Get a MultipleOnlineOutput object.
    data = generic.benchmark(algos=algos, 
                             experiment=e,
                             nsamples=nsamples,
                             process.params=process.params)
    for(algoName in algos) {
      alpha.data[[algoName]]$low[alpha.index] = min(data[[algoName]]$low)
      alpha.data[[algoName]]$high[alpha.index] = max(data[[algoName]]$high)
    }
  }
  draw = list(x=alpha.values, logY=T, logX=F,
              main="Bias learning rate", xlab="alpha", ylab="|| bias ||")
  # 3. Plot low/high curves.
  plot.low.high(alpha.data, e, draw)
}


variance.benchmark.asymptotics <- function() {
  algos = c("sgd.onlineAlgorithm", "implicit.onlineAlgorithm")
  e = normal.experiment(niters=300, p=100)
  e$learning.rate <- function(t) 0.05 / t
  nsamples = 10
}




