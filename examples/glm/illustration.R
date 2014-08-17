## Illustrative example on Section 2.
rm(list=ls())

example.params = list(a=1.5, sigma = 1.5, theta.true = 5.67)
random.params <- function() list(a = runif(1, min=1, 5),
                                 sigma = 5 * runif(1),
                                 theta.true = 20 * runif(1))

theoretical.variance = function(par) with(par, a^2 * sigma^2 / (2 * a -1))
theoretical.variance.sample.based <- function(par, L) {
  with(par, (1 + 1/L) * a^2 * sigma^2 / (2 * a -1))
}

theoretical.variance.plot <- function(npoints, sample.based=NA) {
  var.emp = c()
  var.the = c()
  kIters = 1000
  kSamples = 500
  for(k in 1:npoints) {
    par = random.params()
    m = NA
    V.theor = NA
    if(!is.na(sample.based)) {
      m = run.sgd(par, kIters, kSamples, sample.based.sgd = list(L=2))
      V.theor = theoretical.variance.sample.based(par, L=2)
    } else {
      m = run.sgd(par, kIters, kSamples) 
      V.theor = theoretical.variance(par)
    }
    if(length(m) != kSamples)
      stop("Not correct #samples")
    var.emp = c(var.emp, var(m) * kIters)
    var.the = c(var.the, V.theor)
  }
  plot(var.emp, var.the, pch=20, ylab="Theoretical variance", xlab="Empirical variance")
  lines(var.emp, var.emp, col="red")
}

run.sgd <- function(par, niters, nsamples, sample.based.sgd=NA) {
  M = c()
  for(j in 1:nsamples) {
    y = rnorm(niters, mean=par$theta.true, sd=par$sigma)
    theta.est = 0
    for(i in 1:niters) {
      ai = par$a / i
      if(is.na(sample.based.sgd)) {
        theta.est = theta.est + ai * (y[i] - theta.est) 
      } else {
        L = sample.based.sgd$L
        y.hat = rnorm(L, mean=theta.est, sd=par$sigma)
        theta.est = theta.est + ai * (y[i] - mean(y.hat))
      }
    }
    M[j] = theta.est
  }
  return(M)
}