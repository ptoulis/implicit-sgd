## big glm example
rm(list=ls())
library(ff)
require(biglm)
source("medium-fixed-effects.R")

run.biglm <- function(dim.p, dim.n) {
  # Use biglm
  filename = dataset.filename(dim.p, dim.n, is.big=T)
  print(sprintf("> Loading data from %s", filename))
  df = read.csv.ffdf(file=filename, header=T, first.rows=50000, next.rows=100000)
  print("Data loaded..")
  f = system.time({ mymodel <- biglm(terms(Y ~ 0 + ., data=df), data = df) })
  beta.hat = as.numeric(coef(mymodel))
  dt = f[["elapsed"]]
  #  Does this overflow the RAM?
  # summary(glm(payment ~ sex + age + place.served, data = x[,c("payment","sex","age","place.served")]))
  true.beta = load.beta(dim.p, dim.n, is.big=T)
  mse = vector.dist(true.beta, as.numeric(coef(mymodel)))
  print(sprintf("Big GLM time=%.2f secs and MSE=%.2f", dt, mse))
  return(list(time=dt, beta.hat=beta.hat, 
              mse=vector.dist(beta.hat, load.beta(dim.p, dim.n, is.big=T))))
}

run.sgd <- function(dim.p, dim.n) {
  # 1. Find total number of rows.
  filename = dataset.filename(dim.p, dim.n, is.big=T)
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
  availableRAM = 1024 * 1024**2
  chunkSize = floor(availableRAM / as.numeric(rowSize))
  print(sprintf("> File has %d rows with %.2f kB/row read at %f sec/1k row. Chunk size can be %d. Chunks=", 
                nrows, rowSize/1024,  rowSpeed, floor(chunkSize)))
  chunks = list(start=seq(2, nrows, by=chunkSize))
  m = length(chunks$start)
  chunks$len = rep(chunkSize, m)
  chunks$len[m] = nrows-chunks$start[m] + 1
  print(sprintf("> Total read chunks for SGD = %d", m))
  # 3. Prepare the variables for the main loop.
  n.all = nrows
  p.all = ncol(x) - 1
  beta.old = rep(0, p.all)
  beta.new = rep(0, p.all)
  true.beta = load.beta(dim.p, dim.n, is.big=T)
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
  q.hat=0.0
  iter = 0
  pb = txtProgressBar(style=3)
  ## Main SGD loop
  t0 = proc.time()[["elapsed"]]
  use.explicit = F # what method to use
  print(sprintf("> Using %s updates in SGD. ",
                ifelse(use.explicit, "explicit", "implicit")))
  nonzero.iter = 0
  for(h in 1:length(chunks$start)) {
    print(">> Reading data chunk. Please wait...")
    X = read.matrix(chunks$start[h], chunks$len[h])
    Y = as.vector(X[, p.all+1])
    X = X[, 1:p.all]
    Sx = rowSums(X)
    n = nrow(X) # re-define n=#obs and p=#covariates in this chunk.
    p = ncol(X)
    # 5. Fast update for the intercept
    only.intercept.i = which(Sx==1)
    # print(sprintf("Total %d / %d lines with intercept only.", length(only.intercept.i), n))
    if(length(only.intercept.i) > 0) {
      iter = iter + length(only.intercept.i)
      beta.old[1] = mean(Y[only.intercept.i])
      print(sprintf("> Start with Intercept=%.3f", beta.new[1]))
      # 4. Find optimal a
      X = X[-only.intercept.i, ]
      Y = Y[-only.intercept.i]
      n = nrow(X)
    }

    # q.hat = (1/h) * ((h-1) * q.hat + (sum(X) - 2 * n) / (n * (p-2)))
    # a.optimal = solve.best.alpha(q.hat, p=p)
    # print(sprintf("n=%d q.hat = %.3f -- Optimal a = %.3f", n, q.hat, a.optimal))
    XXt = 0
    a.optimal = 10
    update.alpha = T
    get.alpha.optimal <- function(J) {
      lambdas = eigen(J)$values
      lambdas = lambdas[lambdas > 0]
      if(min(lambdas) < 1e-10) return(9.99)  # hope this is because sample is small.
      optim(par=0, fn=function(x) sum(x^2 * lambdas / (2 * x * lambdas-1)),
            method="L-BFGS-B", lower=1/(2*min(lambdas)) + 1e-5)$par
    }
    
    for(i in 1:n) {
      nonzero.iter = nonzero.iter + 1
      iter = iter + 1
      xi = X[i, ] # current covariate vector
      yi = Y[i]
      Yi.pred = sum(beta.old * xi)
      if(update.alpha) {
        XXt = (1/nonzero.iter) * ((nonzero.iter-1) * XXt + xi %*% t(xi))
        a.new = get.alpha.optimal(XXt)
        # print(sprintf("New optimal a=%.3f", a.new))
        if(abs(a.new - a.optimal) < 5e-2 && a.new != 9.99) {
          update.alpha = F
          print(sprintf("> Best alpha reached at iteration %d...%.3f", 
                        nonzero.iter, a.optimal))
        }
        a.optimal = a.new
      }
      ai = 1 / (1 + (1/a.optimal) * nonzero.iter)
      if(use.explicit) {
        beta.new = beta.old + ai * (yi - Yi.pred) * xi      
      } else {
        ksi = ai * (yi - Yi.pred) / (1 + sum(xi^2) * ai) 
        beta.new = beta.old + ksi * xi # implicit sgd  
      }
      beta.old = beta.new
      setTxtProgressBar(pb, value=iter/n.all)
    }
    print(sprintf("MSE so far %.2f", vector.dist(true.beta, beta.new)))
  }
  t1 = proc.time()[["elapsed"]]
  mse = vector.dist(beta.new, true.beta)
  print(sprintf("Implicit SGD with chunked access %.1f secs, MSE=%.2f", t1-t0, mse))
  return(list(time=t1-t0,
              mse=mse,
              beta.hat=beta.new))
}

run.big.experiment <- function(p=1e2, n=1e6) {
  y = run.biglm(p, n)  # about 1min
  print("")
  print(">> Running Implicit...")
  x = run.sgd(p, dim.n=n) # about 5mins
}

run.experiment.two <- function() {
  # run.biglm(1e4, 1e4) # hangs.
  run.sgd(dim.p=1e4, dim.n=1e4) # about
}

