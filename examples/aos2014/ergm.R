## ERGM code
rm(list=ls())
library(ergm)
data("sampson")
K = 20
all.edges = 2 * choose(K, 2)
theta.true = -1.15
p.true = exp(theta.true) / (1 + exp(theta.true))
V = p.true * (1-p.true) / all.edges


theoretical.var <- function(L, a) (1 + 1/L) * a^2 * V  / (2 * a * V * all.edges - 1)
alpha.opt =  1 / (V * all.edges)
alpha.min = 0.6  / (V * all.edges)


run.many.sgd <- function(nsamples, alpha.inv=0.01, L=1, nruns=1000) {
  all.est <- c()
  pb = txtProgressBar(style=3)
  for(j in 1:nsamples) {
    x = run.sgd(nruns, alpha.inv, L=L, verbose=F)
    all.est <- c(all.est, tail(x, 1))
    setTxtProgressBar(pb, j/nsamples)
  }
  hist(all.est)
  abline(v=theta.true, col="red")
  a = 1 / alpha.inv
  theoretical.variance <- theoretical.var(L, a)
  print(sprintf("empirical variance=%.2f -- theoretical =%.2f", 
                var(all.est) * nruns,
                theoretical.variance))
}

logistic <- function(x) exp(x) / (1+exp(x))

run.sgd <- function(nruns=1000, alpha.inv, L, verbose=T) {
  theta.sgd = c(0)
  e.obs.list = c()
  for(i in 1:nruns) {
    theta.old = tail(theta.sgd, 1)
    p.est <- logistic(theta.old)
    
#     g.sim <- simulate(network(K) ~ edges, coef=c(p.est))
#     g.real <- simulate(network(K) ~ edges, coef=c(p.true))
#     E.est = as.numeric(summary(g.sim ~ edges)) / all.edges
#     E.obs = as.numeric(summary(g.real ~ edges)) / all.edges
    
    E.est = rbinom(L, size=all.edges, prob=p.est) / all.edges 
     E.obs =  rbinom(1, size=all.edges, prob=p.true) / all.edges
  
    # print(sprintf("Obs=%.2f  Real=%.2f", Eold.obs, E.obs))
    ai =   1 / ( 1 + alpha.inv * i)
    theta.sgd = c(theta.sgd, theta.old + ai * (E.obs - mean(E.est)))
    e.obs.list = c(e.obs.list, E.obs)
  }
  if(verbose) {
    plot(theta.sgd, type="l", xlab="iterations", ylab="fitted parameter (theta)",
         main="simple bernoulli example:  p(G) ~ exp(theta * #edges)",
         ylim=c(theta.true, -0.5))
    abline(h=theta.true, col="red")
    legend(x=c(0.8 * nruns), y=c(-0.6), legend = c("real value", "sgd"), 
           col=c("red", "black"), lwd=c(1, 1))
  }
  
  return(theta.sgd)
}
