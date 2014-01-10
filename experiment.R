library(mvtnorm)  # recall rmvnorm(n,...) returns n x p matrix.

base.learning.rate <- function(t, gamma0, alpha, c) {
  CHECK_TRUE(all(c(gamma0, alpha, c) > 0))
  CHECK_INTERVAL(c, 0, 1)
  exp(log(gamma0) - c * log(1 + alpha * gamma0 * t))
}

glm.score.function <- function(h.transfer, theta, datapoint) {
  xt = datapoint$xt
  yt = datapoint$yt
  CHECK_EQ(length(xt), length(theta))
  CHECK_EQ(length(yt), 1, msg="one-dimensional")
  yt.hat = h.transfer(sum(xt * theta))
  with(datapoint, matrix((yt - yt.hat) * xt, ncol=1))
}

empty.experiment <- function(niters) {
  # returns an empty EXPERIMENT object.
  # Useful for initialization and inspection.
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
  e = do.call(function.name, args=list(niters=niters))
  e$name = name
  return(e)
}

normal.experiment <- function(niters) {
  # Normal experiment (linear regression)
  # Defined in Xu (2011), Section 6.2, p.8
  #
  p = 10  # dimension of the parameter vector
  # 1. Define θ*
  experiment = empty.experiment(niters)
  experiment$theta.star = matrix(rep(1, p), ncol=1)  # all 1's
  experiment$p = p
  A = diag(seq(0.01, 1, length.out=p))
  # 2. Define the sample dataset function.
  experiment$sample.dataset = function() {
    epsilon = matrix(rnorm(niters), ncol=1)
    X = rmvnorm(niters, mean=rep(0, p), sigma=A)
    Y = X %*% experiment$theta.star + epsilon
    if(niters > 1000) {
      CHECK_MU0(as.vector(Y), 0)
    }
    CHECK_TRUE(nrow(X) == niters)
    return(list(X=X, Y=Y))
  }
  
  id.fn = function(x) x
  gamma0 = 1 / sum(diag(A))
  # 3. Define the score function
  experiment$score.function = function(theta, datapoint) {
    glm.score.function(h.transfer=id.fn, theta, datapoint)
  }
  
  # 4. Define the learning rate
  experiment$learning.rate <- function(t) {
    base.learning.rate(t, gamma0=gamma0, alpha=1, c=2/3)
  }
  
  # 5. Define the risk . This is usually the negative log-likelihood
  experiment$risk <- function(theta) {
    CHECK_EQ(length(theta), length(experiment$theta.star))
    truth = experiment$theta.star
    tmp = 0.5 * t(theta-truth) %*% A %*% (theta-truth)
    CHECK_EQ(nrow(tmp), 1)
    CHECK_EQ(ncol(tmp), 1)
    CHECK_TRUE(all(tmp >= 0))
    return(as.numeric(tmp))
  }
  
  return(experiment)
}
