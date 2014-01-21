# Copyright (c) 2013
# Panos Toulis, ptoulis@fas.harvard.edu
#
# Benchmarks. Contains code to make benchmark comparisons
# across algorithms. Defines evaluation metrics and creates plots.
#
source("online-algorithms.R")
library(scales)

get.benchmark.filename <- function(prefix, experiment, nsamples, ext) {
  folder = ifelse(length(grep("/n/home", getwd()))==1, "out/odyssey", "out")
  
  filename = sprintf("%s/%s-%s-p%d-t%d-s%d.%s",
                     folder,
                     prefix, experiment$name, 
                     experiment$p, 
                     experiment$niters,
                     nsamples,
                     ext)
  return(filename)
}

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

plot.benchmarkFile <- function(benchmarkFile) {
  # Will plot the low-high values of the benchmark
  # according to the drawing parameters in "draw"
  #
  # Args:
  #   A BenchmarkFile = LIST{benchmark, experiment, draw}
  #   Recall that:
  #     A BENCHMARK is {algo}{low/high} = [] vector of values
  #     A DRAW has information about the drawing (e.g. x-axis etc)
  #
  # Does not return a value, but plots the low/high polygons
  #
  benchmark = NA
  experiment = NA
  draw = NA
  
  if(is.character(benchmarkFile)) {
    if(!file.exists(benchmarkFile)) {
      print(sprintf("File %s does not exist..", benchmarkFile))
      return()
    }
    load(benchmarkFile)
    benchmark = benchmarkFile$benchmark
    experiment = benchmarkFile$experiment
    draw = benchmarkFile$draw
  } else {
    benchmark = benchmarkFile$benchmark
    experiment = benchmarkFile$experiment
    draw = benchmarkFile$draw
  }
  CHECK_benchmark(benchmark)
  # Done loading
  algos = names(benchmark)  # algorithm names (character)
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
    ymin = benchmark[[algoName]]$low
    ymax = benchmark[[algoName]]$high
    if(logY) {
      ymin = log(ymin)
      ymax = log(ymax)
    }
    # Limits in the y-axis
    defaultYlimMin = ifelse(logY, -3, 10^-3)
    defaultYlimMax = ifelse(logY, 3, 10^3)
    ylims = c(min(defaultYlimMin, min(ymin)), min(defaultYlimMax, max(ymax)))
    if(is.element("ylims", names(draw)))
      ylims = draw$ylims
    
    if(i==1) {
      plot(x, ymax, main=title, 
           xlab=xlab,
           ylab=ylab,
           type="l",
           col="white",
           ylim=ylims)
      legend(0.6 * niters, 0.8 * max(ylims), col=cols, legend=algos, lty=1:length(algos))
    }
    polygon(c(x, rev(x)), c(ymin, rev(ymax)), col=alpha(cols[i], 0.2), lty=i)
  }
}

save.benchmarkFile <- function(description,
                               benchmark, experiment,
                               nsamples, draw) {
  benchmarkFile = list(benchmark=benchmark,
                       experiment=experiment,
                       draw=draw)
  filename = get.benchmark.filename(prefix=description,
                                    experiment,
                                    nsamples,
                                    ext="Rdata")
  save(benchmarkFile, file=filename)
}

generic.benchmark <- function(algos, experiment, nsamples,
                              process.params) {
  # Runs a generic benchmark. The idea is the following:
  #   1) Define the algorithms to be tested
  #   2) Define the experiment (sample.dataset, score function, learning rate)
  #   3) Get #nsamples for each experiment
  #   4) For each sample use process.params to process the output.
  #
  # Returns a BENCHMARK object (not defined in terminology)
  # A BENCHMARK object is BENCHMARk{algorithmName}{low/high} = (niters x 1) vector.
  #
  # process.params has two fields (vapply, theta.fn):
  # 
  # Define theta_tj = j-sample of vector θt. Then:
  #  If vapply=TRUE then BENCHMARK{algo}{low} = V, 
  #     where V(t) = 5% quantile (fn(theta_t1), fn(theta_t2)...)
  #
  # If vapply=FALSE then BENCHMARK{algo}{low} = V, 
  #     where V(t) = 5% quantile of fn(theta_tj)  i.e. it is applied in the entire matrix.
  #
  # Args:
  #   algos = vector of algorithms (character)
  #   experiment = EXPERIMENT object.
  #   nsamples = #samples.
  #   process.params = LiST(vapply, fn) - process parameters
  #
  # Returns:
  #   A BENCHMARK object.
  #
  niters = experiment$niters
  # 1. Run the algorithms. Get a MultipleOnlineOutput object.
  mul.out = run.online.algorithm.many(experiment, algos, nsamples=nsamples)
  # 2. Process to get the data
  benchmark = list()
  summary.min = function(x) quantile(x, 0.05)
  summary.max = function(x) quantile(x, 0.95)
  for(algoName in algos) {
    if(process.params$vapply) {
      theta.t.fn <- process.params$theta.fn
      benchmark[[algoName]]$low = mul.OnlineOutput.vapply(experiment, mul.out, algoName,
                                                          theta.t.fn, summary.min)
      benchmark[[algoName]]$high = mul.OnlineOutput.vapply(experiment, mul.out, algoName,
                                                           theta.t.fn, summary.max)
    } else {
      theta.fn <- process.params$theta.fn
      benchmark[[algoName]]$low = mul.OnlineOutput.mapply(experiment, mul.out, algoName, theta.fn)
      benchmark[[algoName]]$high = mul.OnlineOutput.mapply(experiment, mul.out, algoName, theta.fn)
    }
  }
  CHECK_benchmark(benchmark)
  return(benchmark)
}

generic.benchmark.many <- function(algos, experiment.list, nsamples,
                                   process.params) {
  # Runs a generic.benchmark() for every experiment.
  #
  # Args:
  #   algos = vector of algorithms (character)
  #   experiment.list = LIST of experiments
  #   nsamples = #samples.
  #   process.params = LiST(vapply, fn) - process parameters
  #
  # Returns:
  #   A LIST of BENCHMARK objects, one for each EXPERIMENT
  #
  benchmark.list = list()
  for(i in 1:length(experiment.list)) {
    experiment = experiment.list[[i]]
    CHECK_experiment(experiment)
    print(sprintf("Running benchmark for experiment #%d/%d", i, length(experiment.list)))
    benchmark = generic.benchmark(algos, experiment, nsamples, process.params)
    benchmark.list[[i]] = benchmark
  }
  return(benchmark.list)
}

bias.benchmark.asymptotics <- function(p=10, niters=300, nsamples=10) {
  # 0. Define algorithms, basic setup.
  algos = c("sgd.onlineAlgorithm", "implicit.onlineAlgorithm")
  e = normal.experiment(niters=niters, p=p)
  nsamples = nsamples
  
  # 1. Post-processing functions (aggregation)
  dist = function(theta, t) vector.dist(theta, e$theta.star)
  process.params = list(vapply=T, theta.fn=dist)
  
  # 2. Run the algorithms. Get a benchmark object.
  benchmark = generic.benchmark(algos=algos, 
                                experiment=e,
                                nsamples=nsamples,
                                process.params=process.params)
  # 3. Plot low/high curves.
  draw = list(x=1:e$niters, logY=T, logX=F,
              main="Bias asymptotics", xlab="Iterations", ylab="|| bias ||")
  save.benchmarkFile(description="bias-asymp",
                     benchmark, e,
                     nsamples, draw)
  return(benchmark)
}

bias.benchmark.learningRate <- function(p=10, niters=300, nsamples=10) {
  algos = c("sgd.onlineAlgorithm", "implicit.onlineAlgorithm")
  experiment.list = list()
  base.experiment = normal.experiment(niters=niters, p=p)
  nsamples = nsamples
  dist = function(theta, t) vector.dist(theta, base.experiment$theta.star)
  process.params = list(vapply=T, theta.fn=dist)
  
  # 1. Different learning rates to check.
  alpha.values = seq(0.01, 2, length.out=10)
  for(i in 1:length(alpha.values)) {
    experiment.list[[i]] = copy.experiment(base.experiment)
    experiment.list[[i]]$learning.rate <- function(t) alpha.values[i] * base.experiment$learning.rate(t)
  }
  
  for(i in 1:length(alpha.values)) {
    kCurrentLogLevel <<- 0
    t = sample(1:10000, 1)
    CHECK_NEAR(experiment.list[[i]]$learning.rate(t),
               alpha.values[i] * base.experiment$learning.rate(t),
               msg="Check if learning rates are set correctly")
  }

  # 2. Run all benchmarks. Get a LIST of benchmarks.
  benchmark.list = generic.benchmark.many(algos, experiment.list, nsamples, process.params)
  CHECK_EQ(length(experiment.list), length(alpha.values))
 
  ## 3. Populate return object.
  benchmark = list()
  
  for(algoName in algos) {
    benchmark[[algoName]]$low = sapply(1:length(experiment.list),
                                       function(i) {
                                         min(benchmark.list[[i]][[algoName]]$low)
                                       })
    benchmark[[algoName]]$high = sapply(1:length(experiment.list),
                                        function(i) {
                                          max(benchmark.list[[i]][[algoName]]$high)
                                        })
  }
  # 4. Define draw parameters.
  draw = list(x=alpha.values, logY=T, logX=F,
              main="Bias learning rate", xlab="alpha", ylab="|| bias ||")
  
  # 5. Save the benchmark file.
  save.benchmarkFile(description="bias-lr",
                     benchmark, base.experiment,
                     nsamples, draw)
  return(benchmark)
}

variance.benchmark.asymptotics <- function(p=10, niters=300, nsamples=10) {
  # 0. Define algorithms, basic setup.
  kCurrentLogLevel <<- 0
  algos = c("sgd.onlineAlgorithm", "implicit.onlineAlgorithm")
  experiment = normal.experiment(niters=niters, p=p)
  # TODO(ptoulis)
  # Make this into a function.
  J = experiment$Vx
  alpha = experiment$learning.rate(10^9) * 10^9
  loginfo(sprintf("alpha = %.3f", alpha))
  I = diag(experiment$p)
  Sigma.theoretical = alpha^2 * solve(2 * alpha * J - I) %*% J
  CHECK_NEAR(sum(diag(Sigma.theoretical)), sum(diag(experiment$Sigma)))
  CHECK_TRUE(all(eigen(Sigma.theoretical)$values >= 0))
  # 1. Post-processing functions (aggregation)
  
  dist = function(theta.matrix, t) {
    maxT = experiment$niters
    CHECK_EQ(maxT, niters)
    CHECK_EQ(nrow(theta.matrix), experiment$p)
    CHECK_EQ(ncol(theta.matrix), nsamples)
    C = cov(t(theta.matrix))
    matrix.dist(t * C, Sigma.theoretical)
  }
  process.params = list(vapply=F, theta.fn=dist)
  
  # 2. Run the algorithms. Get a MultipleOnlineOutput object.
  benchmark = generic.benchmark(algos=algos, 
                                experiment=experiment,
                                nsamples=nsamples,
                                process.params=process.params)
  # 3. Plot low/high curves.
  draw = list(x=1:experiment$niters, logY=F, logX=F,
              main="Variance asymptotics", xlab="Iterations", ylab="|| Covariance ||")
  
  # 5. Save the benchmark file.
  save.benchmarkFile(description="variance-asymp",
                     benchmark, experiment,
                     nsamples, draw)
  
  return(benchmark)
}

variance.benchmark.learningRate <- function() {
  algos = c("sgd.onlineAlgorithm", "implicit.onlineAlgorithm")
  base.experiment = normal.experiment(niters=400, p=5)
  nsamples = 50
 
  dist = function(theta.matrix, shouldBe.matrix, t) {
    maxT = e$niters
    CHECK_EQ(nrow(theta.matrix), e$p)
    CHECK_EQ(ncol(theta.matrix), nsamples)
    C = cov(t(theta.matrix))
    matrix.dist(t * C, shouldBe.matrix)
  }
  
  # 1. Different learning rates to check.
  alpha.values = seq(0.01, 10, length.out=20)
  benchmark = list()
  
  for(i in 1:length(alpha.values)) {
    e = base.experiment
    alpha = alpha.values[i]
    e$learning.rate <- function(t) alpha * base.experiment$learning.rate(t)
    I = diag(e$p)
    A = alpha^2 * solve(2 * alpha * e$Vx - I) %*% e$Vx
    
    fn = function(theta.matrix, t) dist(theta.matrix, A, t)
    process.params = list(vapply=T, theta.fn=fn)
    
    tmp.bench = generic.benchmark(algos=algos, 
                                  experiment=e,
                                  nsamples=nsamples,
                                  process.params=process.params)
    
    
  }
  
  # 4. Define draw parameters.
  draw = list(x=alpha.values, logY=T, logX=F,
              main="Bias learning rate", xlab="alpha", ylab="|| bias ||")
  
  # 5. Plot low/high curves.
  plot.low.high(benchmark, base.experiment, draw)
}




