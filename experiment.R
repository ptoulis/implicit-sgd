source("terminology.R")
library(mvtnorm)  # recall rmvnorm(n,...) returns n x p matrix.

empty.experiment <- function(niters) {
  return(list(theta.star=matrix(0, nrow=1, ncol=1),
              niters=niters,
              score.function=function(theta, datapoint) {},
              sample.dataset=function() {}))
}

get.experiment <- function(name="normal",
                           niters=1000) {
  # Creates an EXPERIMENT object (see terminology)
  # from the mnemonic name and the specified niters variable
  function.name = sprintf("%s.experiment", name)
  return(do.call(function.name, args=list(niters=niters)))
}

normal.experiment <- function(niters) {
  # Normal experiment (linear regression)
  # Defined in Xu (2011), Section 6.2, p.8
  #
  p = 100  # dimension of the parameter vector
  # 1. Define θ*
  experiment = empty.experiment(niters)
  experiment$theta.star = matrix(rep(1, p), ncol=1)  # all 1's
  A = diag(seq(0.01, 1, length.out=p))
  # 2. Define the sample dataset function.
  experiment$sample.dataset = function() {
    epsilon = matrix(rnorm(niters), ncol=1)
    X = rmvnorm(niters, mean=rep(0, p), sigma=A)
    Y = X %*% experiment$theta.star + epsilon
    CHECK_MU0(as.vector(Y), 0)
    CHECK_TRUE(nrow(X) == niters)
    return(list(X=X, Y=Y))
  }
  # 3. Define the score function
  experiment$score.function = function(theta, datapoint) {
    xt = datapoint$xt
    yt = datapoint$yt
    CHECK_EQ(length(xt), length(theta))
    CHECK_EQ(length(yt), 1, msg="one-dimensional")
    with(datapoint, matrix(yt - sum(xt * theta) * xt, ncol=1))
  }
  return(experiment)
}
