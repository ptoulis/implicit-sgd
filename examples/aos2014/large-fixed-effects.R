# Copyright (c) 2013
# Panos Toulis, ptoulis@fas.harvard.edu
#
# Experiments with big datasets (roughly N * p > 1e9)
# Dataset filenames are defined by dataset.filename() function of datasets.R
# These are stored by default in ./datasets/big 
# 
# Synopsis:
#     run.big.experiment(1e2, 1e5)
#
#  This runs both biglm() and implicit SGD on the dataset with p=1e2 params
#   and N = 1e5 observations.
#
rm(list=ls())
library(ff)
require(biglm)
source("medium-fixed-effects.R")

run.biglm <- function(dim.p, dim.n) {
  # Runs biglm() on the dataset(dim.p, dim.n) using ff package.
  #
  # Returns (time, mse, beta.hat) as a list.
  #
  filename = dataset.filename(dim.p, dim.n, is.big=T)
  print(sprintf("> Loading data from %s", filename))
  # TODO(ptoulis): Haven't found some values that are consistently good across datasets.
  df = read.csv.ffdf(file=filename, header=T, first.rows=10000, next.rows=10000)
  print("> Data loaded..")
  # 1. Fit biglm() using the ffdf frame (loads data successively from the disk)
  f = system.time({ mymodel <- biglm(terms(Y ~ 0 + ., data=df), data = df) })
  # 2. Get the estimates.
  beta.hat = as.numeric(coef(mymodel))
  # 3. Elapsed time.
  dt = f[["elapsed"]]
  # Get ground-truth
  true.beta = load.beta(dim.p, dim.n, is.big=T)
  # 4.. and compute metrics.
  mse = vector.dist(true.beta, as.numeric(coef(mymodel)))
  print(sprintf("Big GLM time=%.2f secs and MSE=%.4f", dt, mse))
  return(list(time=dt, beta.hat=beta.hat, 
              mse=vector.dist(beta.hat, load.beta(dim.p, dim.n, is.big=T))))
}

run.implicit.sgd <- function(dim.p, dim.n) {
  # Runs implicit SGD on the dataset(p, n)
  # The dataset is read is chunks. The # and size of chunks is 
  # determined heuristically in a very crude way (roughly chunk size ~ 1Gb)
  # TODO(ptoulis): Better strategy is possible here.
  # The chunks are fed to the implicit SGD algorithm one-by-one.
  # In small p, we check for sparsity and remove rows that include 
  # only the intercept term, and make the update at once.
  # 
  # The optimal learning rate is determined by sub-sampling the data.
  # This assumes that the no. of covariates p is small-ish (e.g. 1e3)
  # 
  # Rprof(filename="sgd.profile.txt") -- if we need profiling
  filename = dataset.filename(dim.p, dim.n, is.big=T)
  print(sprintf("> Loading data from %s", filename))
  # 1. Get the no. of rows (N) in the dataset
  if(Sys.info()["sysname"] != "Linux") {
    stop("This is supposed to run on a Linux machine. Quitting..")
  }
  cmd.out = system(sprintf("wc -l %s", filename), intern=T)
  nrows = as.numeric(unlist(regmatches(cmd.out, regexpr("\\d+", cmd.out)))) - 1

  read.matrix <- function(start, nlines) {
    # Function to read the dataset in chunks.
    # (start, nlines) determine where to start reading (#line) and 
    # how many lines to read respectively.
    print(sprintf("> Reading file=%s from=%d total=%d lines",
                  filename, start, nlines))
    return(matrix(scan(filename, skip=start-1, sep=",", nlines=nlines, quiet=T), 
                  nrow=nlines, byrow=T))
  }
  
  # 2. Determined how many chunks we need to read.
  nReadRows = 50 # sample that many rows.
  f = system.time({x = read.matrix(2, nReadRows)}) # time access to sample rows
  rowSize = object.size(x) / nReadRows
  rowSpeed = max(0.1, 1000 * f[["elapsed"]] / nReadRows)
  availableRAM = 1024* 1024**2
  chunkSize = floor(availableRAM / as.numeric(rowSize))  # chunk size
  print(sprintf("> File has %d rows with %.2f kB/row read at %f sec/1k row. Chunk size can be %d. Chunks=", 
                nrows, rowSize/1024,  rowSpeed, floor(chunkSize))) # chunks
  chunks = list(start=seq(2, nrows, by=chunkSize))
  nChunks = length(chunks$start)  # total number of chunks
  chunks$len = rep(chunkSize, nChunks)
  chunks$len[nChunks] = nrows-chunks$start[nChunks] + 1
  print(sprintf("> Total read chunks for SGD = %d", nChunks))
  
  # 3. Prepare the variables for the main loop = (p, beta.old, beta,new)
  p = ncol(x) - 1
  beta.old = rep(0, p)
  beta.new = rep(0, p)
  true.beta = load.beta(dim.p, dim.n, is.big=T)
  update.alpha.once = T
  alpha.optimal = 1

  # 4. Define timing and other iteration counters.
  chunk.times = c()  # times to analyze each chunk
  learning.rate.times = c() # times to compute optimal rate for each chunk
  nDatapoints = 0  # current data points fitted.
  
  # 5. Main loop : Load chunk h and then iterate the data
  for(h in 1:length(chunks$start)) {
    print("")
    print(">> Reading new chunk..Please wait...")
    D = read.matrix(chunks$start[h], chunks$len[h])
    print(sprintf("> Data chunk %d/%d read. Calculating sums/means ...etc",
                  h, nChunks))
    x.vars = 1:p
    Sx = rowSums(D[, x.vars]) # take only the sums of the X's
    
    t0 = proc.time()[["elapsed"]]
    # 5a. Fast update for the intercept i.e., xi = (1 0 0 0 0...)
    only.intercept.i = which(Sx==1)
    if(length(only.intercept.i) > 0) {
      w1 = length(only.intercept.i)
      w2 = nDatapoints
      beta.old[1] = (w1 * mean(D[only.intercept.i, p+1]) + w2 * beta.old[1]) / (w1 + w2)
      print(sprintf("Setting b1=%.3f", beta.old[1]))
    }
    
    # 5b. "good examples" = those who do not only "inform" on the intercept.
    good.examples = setdiff(seq(1, nrow(D)), only.intercept.i)
    print(sprintf("> Removed %d lines..", length(only.intercept.i)))
    
    # 5c. Determine optimal learning rate for this chunk.
    #     Sample 10000 rows and then calculate E(x x') and the eigenvalues.
    if(h==1 || !update.alpha.once) {
      print("> Calculating optimal learning rate.")
      lr.t0 = proc.time()[["elapsed"]]
      J = matrix(0, nrow=p, ncol=p)
      ntrials = min(as.integer(0.1 * nrow(D)), 5000)
      D.samples = D[sample(1:nrow(D), size=ntrials, replace=F), ]
      for(i in 1:ntrials) {
        x = c(1, sample(c(1, rbinom(p-2, size=1, prob=0.08))))
        J = (1/i) * ((i-1) * J + x %*% t(x))
      }
      alpha.optimal = best.alpha(J)
      if(alpha.optimal < 1e-3)
        warning("Learning rate is very small. Implicit SGD is potentially unstable.")
      learning.rate.times = c(learning.rate.times, proc.time()[["elapsed"]] - lr.t0)
      print(sprintf("Best alpha = %.1f. Fitting dataset...", alpha.optimal))
    } else {
      print("> Skipping update of learning rate.")
    }
    # 5d. Main loop on "good examples" (see above)
    for(i in good.examples) {
      xi = D[i, x.vars] # current covariate vector
      yi = D[i, p+1]
      Yi.pred = sum(beta.old * xi)
      # learning rate -- need to shift by #datapoints already fit (because of chunks)
      ai = 1 / (1 + (1/alpha.optimal) * (i + nDatapoints))
    
      ksi = ai * (yi - Yi.pred) / (1 + sum(xi^2) * ai)
      beta.new = beta.old + ksi * xi # implicit sgd update
      beta.old = beta.new
    }
    # 5e. Update the count (nDatapoints = #all examples considered so far)
    nDatapoints = nDatapoints + nrow(D)
    
    # 5f. Store times: (t1-t0) = time (secs) to process the current chunk
    t1 = proc.time()[["elapsed"]]
    chunk.times = c(chunk.times, t1-t0)
    print(sprintf("MSE so far %.4f. Total Time =%.2f Learning-rate-Time=%.2f ETA=%.2f (secs)", 
                  vector.dist(true.beta, beta.new), 
                  sum(chunk.times),
                  sum(learning.rate.times),
                  sum(chunk.times) * nChunks / h))
  }
  # 6. Compute metrics and report.
  mse = vector.dist(beta.new, true.beta)
  print(sprintf("Implicit SGD with chunked access %.2f secs, MSE=%.4f", sum(chunk.times), mse))
  return(list(time=sum(chunk.times),
              mse=mse,
              beta.hat=beta.new))
}

run.big.experiment <- function(p=1e2, n=1e6) {
  # Wrapper function to run both biglm() and implicit SGD
  #
  print(">> Running biglm()")
  y = run.biglm(p, n)
  print("")
  print(">> Running Implicit...")
  x = run.implicit.sgd(p, dim.n=n) # about 5mins
}
