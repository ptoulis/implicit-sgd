## big glm example
rm(list=ls())
library(ff)
require(biglm)
source("datasets.R")

run.biglm <- function(dim.p, dim.n) {
  # Use biglm
  filename = dataset.filename(dim.p, dim.n)
  print(sprintf("> Loading data from %s", filename))
  df = read.csv.ffdf(file=filename, header=T, first.rows=50000, next.rows=100000)
  print("Data loaded..")
  f = system.time({ mymodel <- biglm(terms(Y ~ 0 + ., data=df), data = df) })
  beta.hat = as.numeric(coef(mymodel))
  dt = f[["elapsed"]]
  #  Does this overflow the RAM?
  # summary(glm(payment ~ sex + age + place.served, data = x[,c("payment","sex","age","place.served")]))
  return(list(time=dt, beta.hat=beta.hat, 
              mse=vector.dist(beta.hat, load.beta(dim.p, dim.n))))
}

run.sgd <- function(dim.p, dim.n, tol=1e-5) {
  # 1. Find total number of rows.
  filename = dataset.filename(dim.p, dim.n)
  print(sprintf("> Loading data from %s", filename))
  cmd.out = system(sprintf("wc -l %s", filename), intern=T)
  nrows = as.numeric(unlist(regmatches(cmd.out, regexpr("\\d+", cmd.out)))) - 1

  read.matrix <- function(start, nlines) {
    # print(sprintf("> Reading from=%d total=%d lines", start, nlines))
    return(matrix(scan(filename, skip=start-1, sep=",", nlines=nlines, quiet=T), 
                  nrow=nlines, byrow=T))
  }
  
  # 2. Data file storage/access statistics.
  nReadRows = 50
  f = system.time({x = read.matrix(2, nReadRows)})
  rowSize = object.size(x) / nReadRows
  rowSpeed = max(0.1, 1000 * f[["elapsed"]] / nReadRows)
  availableRAM = 512 * 1024**2
  chunkSize = floor(availableRAM / as.numeric(rowSize))
  print(sprintf("> File has %d rows with %.2f kB/row read at %f sec/1k row. Chunk size can be %d. Chunks=", 
                nrows, rowSize/1024,  rowSpeed, floor(chunkSize)))
  chunks = list(start=seq(2, nrows, by=chunkSize))
  m = length(chunks$start)
  chunks$len = rep(chunkSize, m)
  chunks$len[m] = nrows-chunks$start[m] + 1
  print(chunks)
  n = nrows
  p = ncol(x) - 1
  beta.old = rep(0, p) ## initial estimate
  beta.new = beta.old
  # print(sprintf("> Setting learning rates."))
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
  S.optimal = 0
  iter = 0
  pb = txtProgressBar(style=3)
  ## Main SGD loop
  t0 = proc.time()[["elapsed"]]
  use.explicit = F # what method to use
  print(sprintf("> Using %s updates in SGD",
                ifelse(use.explicit, "explicit", "implicit")))
  for(h in 1:length(chunks$start)) {
    X = read.matrix(chunks$start[h], chunks$len[h])
    Y = as.vector(X[, p+1])
    X = X[, 1:p]
    S.optimal = (1/h) * ((h-1) * S.optimal + sum(X) / length(X))
    a.optimal = ifelse(is.infinite(1 / S.optimal), 1.0, 1/S.optimal)
    print(sprintf("Optimal a = %.3f", a.optimal))
    # 2. Set the SGD method (either explicit or implicit)
    for(i in 1:nrow(X)) {
      iter = iter + 1
      ai = 1 / (1 + (1/a.optimal) * iter)
      xi = X[i, ] # current covariate vector
      yi = Y[i]
      xi = xi[1:p]
      Yi.pred = sum(beta.old * xi)
      if(use.explicit) {
        beta.new = beta.old + ai * (yi - Yi.pred) * xi      
      } else {
        ksi = ai * (yi - Yi.pred) / (1 + sum(xi^2) * ai) 
        beta.new = beta.old + ksi * xi # implicit sgd  
      }
      if(iter > 20 && vector.dist(beta.new, beta.old) < tol)
        break;
      beta.old = beta.new
      setTxtProgressBar(pb, value=iter/nrows)
    }
  }
  t1 = proc.time()[["elapsed"]]
  # print(sprintf("SGD with chunked access %.1f secs", t1-t0))
  return(list(time=t1-t0, mse=vector.dist(beta.new, load.beta(p, n)), beta.hat=beta.new))
}
