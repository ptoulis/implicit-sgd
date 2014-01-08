## Utils. Sampling, distance metrix etc.

## returns how many 1's if total size=p
sparsity = function(p) return(max(as.integer(0.1 * p) , 2))
##  Samples a setup
## Through the max rate we can control the magnitude of the samples.
sample.betas = function(p, max.rate=5) {
    #t.betas=c( 0.35, 0.86, -0.49, 0.41)
    t.betas= 2 *  runif(p)-1
    t.betas[1] = 0 ## fix one value for identifiability
    active.no = sparsity(p)
    ## Rescale in order to control the maximum rate.
    m = sum(tail(sort(t.betas[2:p]), active.no))
    #scale = ifelse(m>log(max.rate), log(max.rate)/m, 1)
    add = log(max.rate) - m
    t.betas[2:p] = t.betas[2:p]+add/(active.no)
    return(as.matrix(t.betas, ncol=1))
}



## Samples a 1xp  matrix.
sample.x <- function(p) {
    how.many = sparsity(p)
    zeroes = rep(0, p)
    which.ones = sample(1:p, size=how.many)
    zeroes[which.ones] = 1
    return(t(as.matrix(zeroes)))
}

##   Samples the counts    
##   nx1   vector.
sample.y <- function(xt, betas) {
  mu = exp( xt %*% betas    )   ## linear predictor
  yt = rpois(1, lambda=mu)  # observations
  return(yt)
}   

## Computes MLE 
## Calls optim, might be slow.
mle <- function(Y, X) {
  b0 = rep(0, ncol(X))
  cost <- function(b.tmp) {
    return(-log.lik(Y,X, b.tmp))
  }
  fit <- optim(b0, cost, method="BFGS")
  return(fit$par)
}

### Log likelihood
log.lik <- function(yt, xt, fitted.b) {
  mus.tmp <- exp( xt %*% fitted.b)
  return(dpois(yt, lambda=mus.tmp, log=T))
}
std.dist <- function(y, x.center) {
    return(sum((x.center-y)^2) / sum(x.center^2))
}
####   Learning rate schedules.
typical.eta.t = function(t) {
  5 * min(0.01, 1/t)
}


###    Plots  here. 
## Distance to the real values.
plot.metric = function(what, sgd, impl) {
    par(mfrow=c(1,1))
    plot(cumsum(sgd$paths[[what]]), type="l", ylab=sprintf("%s", what), 
         main=sprintf("%s paths", what), col="red")
    lines(cumsum(impl$paths[[what]]), col="blue")
}







