# Copyright (c) 2013
# Panos Toulis, ptoulis@fas.harvard.edu
#
# Benchmarks. Contains code to make benchmark comparisons
# across algorithms. Defines evaluation metrics and creates plots.
#
source("online-algorithms.R")
library(scales)

limit.variance <- function(experiment) {
  J = experiment$J
  a = experiment$learning.rate(10^9) * 10^9
  if(a > 10^6) stop("There is a problem most likely")
  I = diag(experiment$p)
  return(a^2 * solve(2 * a * J - I) %*% J)
}

get.out.folder <- function() {
  # The "out" folder is where we save Rdata (benchmarks raw data)
  # and plots (png files)
  #
  folder = ifelse(length(grep("/n/home", getwd()))==1, "out/odyssey", "out")
  return(folder)
}

get.benchmark.filename <- function(prefix, benchmark, ext) {
  folder = get.out.folder()  
  filename = sprintf("%s/%s-%s-p%d-t%d-s%d.%s",
                     folder,
                     prefix,
                     benchmark$experiment$name, 
                     benchmark$experiment$p, 
                     benchmark$experiment$niters,
                     benchmark.nsamples(benchmark),
                     ext)
  return(filename)
}

plot.all.benchmarks <- function() {
  # Goes through all the Rdata files, and then plots the BENCHMARK object
  # for each one of them
  #
  all.files = list.files(get.out.folder(), full.names=T)
  rdata.files = all.files[grep("Rdata", all.files)]
  kCurrentLogLevel <<- 0
  for(file in rdata.files) {
    loginfo(sprintf("Loading benchmark %s for plotting. Please wait..", file))
    plot.benchmark(file, T)
    loginfo("Done.")
  }
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

plot.benchmark <- function(benchmarkObjectORFile, toPng=F) {
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
  
  if(is.character(benchmarkObjectORFile)) {
    if(!file.exists(benchmarkObjectORFile)) {
      print(sprintf("File %s does not exist..", benchmarkObjectORFile))
      return()
    }
    load(benchmarkObjectORFile)
    experiment = benchmark$experiment
    draw = benchmark$draw
  } else {
    benchmark = benchmarkObjectORFile
    experiment = benchmark$experiment
    draw = benchmark$draw
  }
  CHECK_benchmark(benchmark)
  # Done loading
  algos = benchmark.algos(benchmark)  # algorithm names (character)
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
  if(toPng) {
    png(file=get.benchmark.filename(sprintf("plot-%s", benchmark$name),
                                    benchmark=benchmark, ext="png"))
  }
  
  for(i in 1:length(algos)) {
    algoName = algos[i]
    ymin = benchmark.algo.low(benchmark, algoName)
    ymax = benchmark.algo.high(benchmark, algoName)
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
  if(toPng)
    graphics.off()
}

save.benchmark <- function(description, benchmark) {
  # Saves the BenchmarkFile object.
  # Gets the name from get.benchmark.filename()
  #
  # Args:
  #   description: A character vector as a filename prefix
  #   multipleOut: A MultipleOnlineOutput object
  #   benchmark: A BENCHMARK object.
  #   experiment : EXPERIMENT (see terminology for all)
  # 
  CHECK_benchmark(benchmark)
  filename = get.benchmark.filename(prefix=description,
                                    benchmark,
                                    ext="Rdata")
  save(benchmark, file=filename)
}

summarize.benchmark.list <- function(benchmark.list) {
  if(length(benchmark.list) == 1)
    return(benchmark.list[[1]])
  
  all.index = 1:length(benchmark.list)
  b1 = benchmark.list[[1]] 
  data.lowHigh = list()
  algos = benchmark.algos(b1)
  
  for(algoName in algos) {
    data.lowHigh[[algoName]]$low = 
      sapply(all.index, function(i) {
        b = benchmark.list[[i]]
        min(benchmark.algo.low(b, algoName))
      })
    data.lowHigh[[algoName]]$high = 
      sapply(all.index, function(i) {
        b = benchmark.list[[i]]
        max(benchmark.algo.high(b, algoName))
      })
  }
  
  return(list(name=b1$name, mulOut=b1$mulOut, lowHigh=data.lowHigh, 
              experiment=b1$experiment))
}

best.scale <- function() {
  alpha.values = seq(0.001, 2.5, length.out=100)
  vars = sapply(alpha.values, function(a) {
    e = normal.experiment(niters=100, p=10, lr.scale=a)
    return(sum(diag(limit.variance(e))))
  })
  # print(vars)
  neg = which(vars < 0)
  vars = vars[-neg]
  alpha.values = alpha.values[-neg]
  plot(alpha.values, vars, type="l")
  i = order(vars)
  abline(v=alpha.values[i[1]], col="red")
  return(alpha.values[i])
}

# CORE functionality of benchmarks.R
execute.benchmarks <- function(mulOutParams.list, processParams.list) {
  # Returns a LIST of BENCHMARK objects 
  #
  # Args:
  #   mulOutParams.list = LIST of MultipleOnlineParams object
  #     This defines the different experiments to run (e.g. possibly diferent learning rates)
  #   processParams.list = LIST of processParams
  #     Each one essentially defines a bias or variance distance.
  #     A bit unecessary because for each run.benchmark.* function
  #     we are using the same list of processParams. However, the names of 
  #     the benchmarks are obtained from the "name" field in that list.
  #
  # Returns a LIST of BENCHMARK objects, one for each processParams.list argument.
  #
  benchmark.list.out = list()
  
  isLearningRateBenchmark = length(mulOutParams.list) > 1
  
  # Object to hold the multiple online output
  for(j in 1:length(mulOutParams.list)) {
    params = mulOutParams.list[[j]]
    experiment = params$experiment
    algos = params$algos
    nsamples = params$nsamples
    print(sprintf("LR-benchmarks=%s -- Running new experiment (%d/%d): %s",
                  isLearningRateBenchmark,
                  j, length(mulOutParams.list),
                  get.experiment.description(experiment)))
    # 1. Run j-th experiment. Obtain multiple online output.
    mulOut = run.online.algorithm.many(experiment, algos, nsamples)
    
    # 2. For every benchmark, create the BENCHMARK object
    #   for this particular MultipleOnlineOutput object.
    for(m in 1:length(processParams.list)) {
      benchName = processParams.list[[m]]$name
      processParams = processParams.list[[m]]
      # 2. Process the raw data using the processParams specification.
      data.lohi = list()
      summary.min = function(x) quantile(x, 0.05)
      summary.max = function(x) quantile(x, 0.95)
      for(algoName in algos) {
        if(processParams$vapply) {
          # Hard-coding the transformation function
          theta.t.fn <- default.bias.dist(experiment)
          # Computing low/high values.
          data.lohi[[algoName]]$low = mul.OnlineOutput.vapply(mulOut, experiment, algoName,
                                                              theta.t.fn, summary.min)
          data.lohi[[algoName]]$high = mul.OnlineOutput.vapply(mulOut, experiment, algoName,
                                                               theta.t.fn, summary.max)
        } else {
          theta.fn <- default.var.dist(experiment, nsamples)
          if(isLearningRateBenchmark) {
            print(sprintf("Variance LR benchmark. dist() will be average trace."))
            theta.fn <- function(theta.matrix, t) {
              C = t * cov(t(theta.matrix))
              sum(diag(C)) / nrow(C) # average trace.
            }
          }
          data.lohi[[algoName]]$low = mul.OnlineOutput.mapply(mulOut, experiment, algoName, theta.fn)
          data.lohi[[algoName]]$high = mul.OnlineOutput.mapply(mulOut, experiment, algoName, theta.fn)
        }
      }
      
      benchmark = list(name=benchName, mulOut=mulOut, lowHigh=data.lohi, experiment=experiment)
      CHECK_benchmark(benchmark)
      L = length(benchmark.list.out[[benchName]])
      # 3. Save the BENCHMARK object for each one.
      benchmark.list.out[[benchName]][[L + 1]] <- benchmark
    }
  }
  # Final step. Now the benchmark list is {benchmarkName} = LIST{b1, b2, b...}
  # where b_i = benchmark
  # If we have multiple benchmarks that means there were multiple experiments.
  # so, we need to summarize.
  for(benchName in names(benchmark.list.out)) {
    benchmark.list.out[[benchName]] = 
      summarize.benchmark.list(benchmark.list.out[[benchName]])
  }
  return(benchmark.list.out)
}

run.benchmark.learningRate <- function(niters=100, p=100, nsamples=10,
                                       max.lr.scale=2.0, nlr.scales=2) {
  #  1. Define processParams.
  var.processParams = list(name="variance-LR", vapply=F)
  bias.processParams = list(name="bias-LR", vapply=T)
  
  # 2. Define experiments
  base.experiment = normal.experiment(niters=niters, p=p)
  mulOutParams.list = list()
  
  # 2. Different learning rates to check.
  alpha.values = seq(0.01, max.lr.scale, length.out=nlr.scales)
  for(i in 1:length(alpha.values)) {
    kCurrentLogLevel <<- 0
    experiment = normal.experiment(niters=niters, p=p,
                                   lr.scale=alpha.values[i])
    
    loginfo(sprintf("Learning rate at 10 for alpha=%.2f is %.4f",
                    alpha.values[i], experiment$learning.rate(10)))
    mulOutParams.list[[i]] = list(experiment=experiment,
                                  nsamples=nsamples,
                                  algos=c(kSGD, kIMPLICIT))
  }
  # 2b. Checking to see everything is fine w/learning rates.
  for(alpha.i in 1:length(alpha.values)) {
    kCurrentLogLevel <<- 0
    sample.t = sample(1:niters, 1)
    lr = mulOutParams.list[[alpha.i]]$experiment$learning.rate(sample.t)  # learning-rate
    lr.shouldBe = base.experiment$learning.rate(sample.t) * alpha.values[alpha.i]
    
    CHECK_NEAR(lr, lr.shouldBe, msg="Check if learning rates are set correctly")
  }
  
  # 3. Define input for run.generic.benchmarks()
  # mulOutParams.list already defined.
  processParams.list = list(bias.processParams, var.processParams)
  
  # 3. Run the algorithms. Get a MultipleOnlineOutput object.
  benchmark.list = execute.benchmarks(mulOutParams.list, processParams.list)
  
  # 4. (OPTIONAL) Define draw params -- can be changed later.
  for(benchmarkName in names(benchmark.list)) {
    draw = list(x=alpha.values, logY=F, logX=F,
                main="Variance/Learning rate (lr)", xlab="alpha", ylab="|| Covariance ||")
    benchmark = benchmark.list[[benchmarkName]]
    # Draw params for bias.
    if(length(grep("bias", benchmarkName)) > 0) {
      draw$logY = T
      draw$ylab = "|| bias ||"
      draw$main = "Bias/Learning rate (lr)"
    }
    benchmark$draw = draw
    # 5. Save the benchmark file.
    save.benchmark(description=benchmarkName, benchmark) 
  }
}

run.benchmark.asymptotics <- function(niters=100, p=100, nsamples=10) {
  
  mulOutParams = list(algos = c(kSGD, kIMPLICIT),
                      experiment = normal.experiment(niters=niters, p=p),
                      nsamples = nsamples)
  
  # 0. Define algorithms, basic setup.
  CHECK_mulOutParams(mulOutParams)
  experiment = mulOutParams$experiment
  Sigma.theoretical = limit.variance(experiment)
  
  CHECK_TRUE(all(eigen(Sigma.theoretical)$values >= 0))
  
  # Define processParams.
  var.processParams = list(name="variance-asymp", vapply=F)
  bias.processParams = list(name="bias-asymp", vapply=T)
  
  # 2. Define input for run.generic.benchmarks()
  mulOutParams.list = list(mulOutParams)  # the experiments to run
  processParams.list = list(bias.processParams, var.processParams)
  
  # 3. Run the algorithms. Get a MultipleOnlineOutput object.
  benchmark.list = execute.benchmarks(mulOutParams.list, processParams.list)
  
  # 4. (OPTIONAL) Define draw params -- can be changed later.
  # 4. Add draw parameters and save
  for(benchmarkName in names(benchmark.list)) {
    draw = list(x=1:experiment$niters, logY=F, logX=F,
                main="Variance asymptotics", xlab="Iterations", ylab="|| Covariance ||")
    
    benchmark = benchmark.list[[benchmarkName]]
    # Draw params for bias.
    if(length(grep("bias", benchmarkName)) > 0) {
      draw$logY = T
      draw$ylab = "|| bias ||"
    }
    benchmark$draw = draw
    # 5. Save the benchmark file.
    save.benchmark(description=benchmarkName, benchmark) 
  }
}

