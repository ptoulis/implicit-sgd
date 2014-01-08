source("../r-toolkit/checks.R")
source("../r-toolkit/logs.R")

# Assume we have a parametric statistical model with a p x 1 parameter "theta"
#
# Define DATAPOINT to be a list(xt, yt) where xt = px1 vector, yt=1dim value
# Then a DATASET is a list (X, Y) where X=(niters x p) and Y = niters x 1
#
# Create an EXPERIMENT. This is comprised by 
#   "theta.star" = # (px1) vector of real values, 
#   niters = (#iterations), 
#   sample.dataset() = samples DATASET object
#   score.function = \nabla loglik(theta, data_t) = score function
#
# An OnlineAlgorithm is defined by a function that accepts a DATASET 
# and returns an OnlineOutput object:
#   estimates = (niters x p) matrix of estimate i.e. (theta_t)
#
# An OnlineExperiment is defined by
#   OnlineAlgorithm + EXPERIMENT
#   Risk(theta) = function that computes the expected risk for some a vector theta
#
# Some functions for datasets/datapoint

get.dataset.point <- function(dataset, t) {
  CHECK_TRUE(nrow(dataset$X) >= t && nrow(dataset$Y) >= t)
  return(list(xt=dataset$X[t, ],
              yt=dataset$Y[t]))
}

CHECK_dataset <- function(dataset) {
  CHECK_SETEQ(names(dataset), c("X", "Y"))
  CHECK_columnVector(dataset$Y)
  CHECK_EQ(nrow(dataset$X), nrow(dataset$Y))
  CHECK_numeric(dataset$X)
  CHECK_numeric(dataset$Y)
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
  CHECK_MEMBER(names(experiment), c("theta.star", "niters", 
                                    "sample.dataset",
                                    "score.function"))
  CHECK_columnVector(experiment$theta.star)
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
}

CHECK_OnlineOutput <- function(onlineOut) {
  CHECK_MEMBER(names(onlineOut), c("estimates"))
  CHECK_numeric(onlineOut$estimates)
}

