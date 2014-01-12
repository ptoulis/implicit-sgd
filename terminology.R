rm(list=ls())
source("../r-toolkit/checks.R")
# source("../r-toolkit/logs.R")

# Assume we have a parametric statistical model with a p x 1 parameter "theta"
#
# Define DATAPOINT to be a list(xt, yt) where xt = px1 vector, yt=1dim value
# Then a DATASET is a list (X, Y) where X=(niters x p) and Y = niters x 1
#
# Create an EXPERIMENT (stochastic approximation). This is comprised by 
#   "theta.star" = # (px1) vector of real values, 
#   niters = (#iterations), 
#   sample.dataset() = samples DATASET object
#   score.function = \nabla loglik(theta, data_t) = score function
#   learning.rate = function (t, ...) = gives the learning rate > 0
#
# An OnlineAlgorithm is defined by a function that accepts 
#   t = no. of iteration
#   theta.old = vector of current estimate
#   datapoint  = DATAPOINT object at t (xt vector, yt value)
#   experiment = EXPERIMENT object, has learning rate/score function
# The object returns an OnlineOutput object:
#   estimates = (p  x niters) matrix of estimate i.e. (theta_t)
#   last = last vector of estimates
#
# A Benchmark is defined by
#   OnlineAlgorithm(s) + EXPERIMENT
#   Risk(theta) = function that computes the expected risk for some a vector theta
#
# Some functions for datasets/datapoint

get.dataset.point <- function(dataset, t) {
  CHECK_TRUE(nrow(dataset$X) >= t && nrow(dataset$Y) >= t)
  return(list(xt=dataset$X[t, ],
              yt=dataset$Y[t, ]))
}

dataset.size <- function(dataset) {
  # Return a LIST with the #samples and the #parameters
  return(list(nsamples=nrow(dataset$X),
              p=ncol(dataset$X)))
}

empty.onlineOutput <- function(dataset)  {
  size = dataset.size(dataset)  
  return(list(estimates=matrix(0, nrow=size$p, ncol=size$nsamples)))
}

add.estimate.onlineOutput <- function(out, t, estimate) {
  # Adds an estimate in the OnlineOutput object
  CHECK_TRUE(t <= ncol(out$estimates), msg="t < #total samples")
  CHECK_EQ(length(estimate), nrow(out$estimates), msg="Correct p=#parameters")
  out$estimates[, t] <- estimate
  return(out)
}

onlineOutput.estimate <- function(out, t) {
  if(t==0) {
    warning("Default is to return 1-vector, for t=0")
    return(matrix(1, nrow=nrow(out$estimates), ncol=1))
  }
  CHECK_TRUE(t <= ncol(out$estimates), msg="t < #total samples")
  return(matrix(out$estimates[, t], ncol=1))
}

CHECK_dataset <- function(dataset) {
  CHECK_SETEQ(names(dataset), c("X", "Y"))
  CHECK_columnVector(dataset$Y)
  CHECK_EQ(nrow(dataset$X), nrow(dataset$Y))
  CHECK_numeric(dataset$X)
  CHECK_numeric(dataset$Y)
}

CHECK_onlineAlgorithm = function(algorihm) {
  CHECK_SETEQ(names(formals(algorithm)), c("t", "data.history", "experiment", "online.out"))
}

CHECK_columnVector <- function(x) {
  CHECK_GE(nrow(x), 1, msg="Has >=1 row")
  CHECK_EQ(ncol(x), 1, msg="Only one col")
}

CHECK_rowVector <- function(x) {
  CHECK_columnVector(t(x))
}

CHECK_numeric <- function(x) {
  # Checks that there are no NA's, Inf's etc
  CHECK_TRUE(all(!is.na(x)), msg="No NAs")
  CHECK_TRUE(all(!is.infinite(x)), msg="No Infs")
  CHECK_TRUE(is.numeric(x), msg="Should be numeric")
}

CHECK_experiment <- function(experiment) {
  CHECK_MEMBER(names(experiment),
               c("name", "p", "theta.star", "niters", 
                 "sample.dataset", "Vx",
                 "score.function",
                 "learning.rate",
                 "risk"), msg="Correct fields for the experiment")
  CHECK_columnVector(experiment$theta.star)
  CHECK_EQ(experiment$p, length(experiment$theta.star))
  D = experiment$sample.dataset()
  CHECK_dataset(D)
  point = get.dataset.point(D, 2)
  CHECK_SETEQ(names(formals(experiment$score.function)),
              c("theta", "datapoint"))
  p = length(experiment$theta.star)
  bad.theta = rep(0, length(experiment$theta.star) - 1)
  # we send a theta with wrong dimension
  CHECK_EXCEPTION(experiment$score.function(bad.theta, point))
  CHECK_columnVector(experiment$score.function(experiment$theta.star, point))
  
  CHECK_MEMBER("t", names(formals(experiment$learning.rate)))
  CHECK_EQ(experiment$risk(experiment$theta.star), 0, msg="Risk of ground-truth is 0")
  for(i in 1:100) {
    theta.random = runif(length(experiment$theta.star), min=-2, max=2)
    CHECK_TRUE(experiment$risk(theta.random) >= 0, msg="Loss has to be positive.")
  }
  # check that the covariance matrix is positive definite
  CHECK_TRUE(all(eigen(experiment$Vx)$values > 0), msg="Cov matrix is nonneg def")
}

CHECK_onlineOutput <- function(onlineOut) {
  CHECK_MEMBER(names(onlineOut), c("estimates", "last"))
  CHECK_numeric(onlineOut$estimates)
  CHECK_EQ(onlineOut$estimates[ , ncol(onlineOut$estimates)], onlineOut$last)
}

CHECK_benchmark <- function(benchmark) {
  CHECK_MEMBER(names(benchmark), c("onlineAlgorithms", "experiment"))
  e = benchmark$experiment
  CHECK_experiment(e)
  CHECK_EQ(benchmark$risk(e$theta.star), 0, msg="Risk for benchmark")
}
