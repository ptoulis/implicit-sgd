##   sgd.R
##  Implements SGD method for Poisson Regression model


##  SGD algorithm for Poisson regression model.
sgd.update <- function(beta.old, t, yt, xt) {
  hat.mut = as.numeric(exp(xt %*% beta.old))
  eta.t = typical.eta.t(t)
  beta.new = beta.old + eta.t * (yt - hat.mut) * t(xt)
  return(beta.new)
}
