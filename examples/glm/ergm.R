## ERGM code
rm(list=ls())
library(ergm)
data("sampson")
K = 18
all.edges = 2 * choose(K, 2)
theta.true = -0.9
p.true = exp(theta.true) / (1 + exp(theta.true))
V = p.true * (1-p.true) / all.edges


theoretical.var <- function(L, a) (1 + 1/L) * a^2 * V  / (2 * a * V * all.edges - 1)
alpha.opt =  1 / (V * all.edges)
alpha.min = 0.6  / (V * all.edges)


run.many.sgd <- function(nsamples, alpha.inv=0.01, L=1, nruns=1000, verbose=F) {
  all.est <- c()
  pb = txtProgressBar(style=3)
  for(j in 1:nsamples) {
    x = run.sgd(nruns, alpha.inv, L=L)
    all.est <- c(all.est, tail(x, 1))
    setTxtProgressBar(pb, j/nsamples)
  }

  abline(v=theta.true, col="red")
  a = 1 / alpha.inv
  theoretical.variance <- theoretical.var(L, a)
  if(verbose) {
    hist(all.est)
      print(sprintf("empirical variance=%.3f -- theoretical =%.3f", 
                    var(all.est) * nruns,
                    theoretical.variance))
  }
  return(var(all.est) * nruns)
}

logistic <- function(x) exp(x) / (1+exp(x))

plot.sgd <- function(alpha.inv.list, var.samples=100) {
  nruns = 1000
  m = length(alpha.inv.list)
  cols = rainbow(m)
  labels = c() #sapply(alpha.inv.list, function(i) sprintf("a = %5.2f", 1/i))
  nplotPoints = 30
  plot(1:nplotPoints, rep(theta.true, nplotPoints), type="l", col="black",
       xlab=sprintf("iterations (1 to %d)", nruns), ylab="fitted parameter (theta)",
       xaxt="n", lty=2,
       main="simple bernoulli example:  p(G) ~ exp(theta * #edges)",
       ylim=c(theta.true-0.1, -0.6))
  
  for(i in 1:m) {
    a.inv = alpha.inv.list[i]
    theta.sgd = run.sgd(nruns=nruns, alpha.inv = a.inv, L = 1)
    lines(theta.sgd[seq(1, nruns, length.out=nplotPoints)], col=cols[i], type="b", pch=i)
    var.emp = run.many.sgd(nsamples = var.samples, alpha.inv = a.inv, L = 1, nruns = nruns)
    labels = c(labels, sprintf("a = %6.2f (var = %4.2f)", 1/a.inv, var.emp))
  }
  legend(x=c(0.6 * nplotPoints), y=c(-0.6), legend = c(labels), 
         col=cols, pch=c(1:m), lwd=c(rep(1.5, m)))
}
  
run.sgd <- function(nruns=1000, alpha.inv, L, use.binom=T) {
  theta.sgd = c(0)
  e.obs.list = c()
  for(i in 1:nruns) {
    theta.old = tail(theta.sgd, 1)
    p.est <- logistic(theta.old)
    
    E.est = c()
    E.obs = c()
    if(!use.binom) {
      g.sim <- simulate(network(K) ~ edges, coef=c(p.est))
      g.real <- simulate(network(K) ~ edges, coef=c(p.true))
      E.est = as.numeric(summary(g.sim ~ edges)) / all.edges
      E.obs = as.numeric(summary(g.real ~ edges)) / all.edges
    } else {
      E.est = rbinom(L, size=all.edges, prob=p.est) / all.edges 
      E.obs =  rbinom(1, size=all.edges, prob=p.true) / all.edges
    }
    # print(sprintf("Obs=%.2f  Real=%.2f", Eold.obs, E.obs))
    ai =   1 / ( 1 + alpha.inv * i)
    theta.sgd = c(theta.sgd, theta.old + ai * (E.obs - mean(E.est)))
    e.obs.list = c(e.obs.list, E.obs)
  }
  return(theta.sgd)
}
