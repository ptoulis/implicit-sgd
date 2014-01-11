# Online algorithms.
source("terminology.R")
source("experiment.R")

run.onlineAlgorithm <- function(dataset, experiment, algorithm) {
  # Implements the vanilla SGD algorithm.
  CHECK_dataset(dataset)
  CHECK_experiment(experiment)
  out = empty.onlineOutput(dataset)
  nsamples = dataset.size(dataset)$nsamples
  pb = txtProgressBar(style=3)
  algo.name = as.character(substitute(algorithm))
  cat(sprintf("Running algorithm %s, Experiment=%s, samples=%d \n",
              algo.name, experiment$name, experiment$niters))
  for(t in 1:nsamples) {
    ## Run all iterations here.
    # History has data from 1 to t
    # Recall that theta_t = estimate of theta AFTER seeing datapoint t.
    history = list(X=matrix(dataset$X[1:t, ], ncol=experiment$p),
                   Y=matrix(dataset$Y[1:t, ], ncol=1))
    # 1. Runs the online-algorithm step.
    theta.new = algorithm(t, online.out=out,
                          data.history=history,
                          experiment=experiment)
    out <- add.estimate.onlineOutput(out, t, estimate=theta.new)
    setTxtProgressBar(pb, value=t/nsamples)
  }

  if(length(grep("asgd", algo.name)) > 0) {
    ## This is ASGD need to rework the estimates
    warning("Transforming the ASGD estimates")
    estimates = out$estimates
    avg.estimates = matrix(0, nrow=nrow(estimates), ncol(estimates))
    for(t in 2:ncol(estimates)) {
      avg.estimates[,t] = (1-1/t) * avg.estimates[,t-1] + (1/t) * estimates[,t]
    }
    # y = t(apply(X, 1, function(s) cumsum(s)))
    out$estimates = avg.estimates
  }
  out$last = out$estimates[, nsamples]
  return(out)
}

sgd.onlineAlgorithm <- function(t, online.out, data.history, experiment) {
  # Implements the SGD algorithm
  #
  datapoint = get.dataset.point(dataset=data.history, t=t)
  at = experiment$learning.rate(t)
  theta.old = onlineOutput.estimate(online.out, t-1)
  score.t = experiment$score.function(theta.old, datapoint)
  theta.new = theta.old + at * score.t
  return(theta.new)
}

asgd.onlineAlgorithm <- function(t, online.out, data.history, experiment) {
  # Implements the ASGD algorithm (Polyak 1992)
  #
  return(sgd.onlineAlgorithm(t, online.out, data.history, experiment))
}

batch.onlineAlgorithm <- function(t, online.out, data.history, experiment) {
  fit = lm(Y ~ X, data=data.history)
  theta.new = tail(as.numeric(fit$coefficients), experiment$p)
  return(theta.new)
}
oracle.onlineAlgorithm <- function(t, online.out, data.history, experiment) {
  return(experiment$theta.star)
}

implicit.onlineAlgorithm <- function(t, online.out, data.history, experiment) {
  datapoint = get.dataset.point(dataset=data.history, t=t)
  at = experiment$learning.rate(t)
  xt = datapoint$xt
  yt = datapoint$yt
  theta.old = onlineOutput.estimate(online.out, t-1)
  get.score.coeff <- function(theta) {
    # this returns the value  yt - h(theta' xt)  -- for a GLM
    norm.score = experiment$score.function(theta, datapoint) / xt
    CHECK_TRUE(all(abs(norm.score - norm.score[1]) < 1e-5))  # all are the same
    return(norm.score[1])  
  }
  # 1. Define the search interval
  rt = at * get.score.coeff(theta.old)
  Bt = c(0, rt)
  if(rt < 0) {
    Bt <- c(rt, 0)
  }
  
  implicit.fn <- function(u) {
    u  - at * get.score.coeff(theta.old + u * xt)
  }
  # 2. Solve implicit equation
  xit = uniroot(implicit.fn, interval=Bt)$root
  theta.new = theta.old + xit * xt
  return(theta.new)
}