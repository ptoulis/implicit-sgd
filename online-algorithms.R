# Online algorithms.
source("terminology.R")
source("experiment.R")

baseline.learning.rate = function(t) min(0, 1/t)

run.online.algorithm <- function(dataset, experiment, algorithm) {
  # Implements the vanilla SGD algorithm.
  CHECK_dataset(dataset)
  CHECK_experiment(experiment)
  out = empty.onlineOutput(dataset)
  nsamples = dataset.size(dataset)$nsamples
  pb = txtProgressBar(style=3)
  cat(sprintf("Running SGD, Experiment=%s, samples=%d \n",
              experiment$name, experiment$niters))
  for(t in 2:nsamples) {
    ## Run all iterations here.
    point = get.dataset.point(dataset, t)
    theta.old = onlineOutput.estimate(out, t-1)
    # 1. Runs the online-algorithm step.
    theta.new = algorithm(t, online.out=out,
                          datapoint=point,
                          experiment=experiment)
    # 2. Debugging info
    # debugging info. To be removed.
    logdebug(sprintf("Iteration %d  = %.5f", t))
    logdebug("xt=")
    logdebug(round(point$xt, 3))
    logdebug(sprintf("yt=%.3f", point$yt))
    logdebug(sprintf("yt-xt * theta_t = %.3f", point$yt - sum(point$xt * theta.old)))
    logdebug("Score=")
    logdebug(round(t(score.t), 3))
    logdebug("Old theta vector (theta_t)")
    logdebug(round(t(theta.old), 3))
    logdebug("New theta vector (theta_t+1)")
    logdebug(round(t(theta.new), 3))
    logdebug("______________________________________________")
    # 3. Save the update to the out object (OnlineOutput)
    out <- add.estimate.onlineOutput(out, t, estimate=theta.new)
    setTxtProgressBar(pb, value=t/nsamples)
  }
  
  out$last = out$estimates[, nsamples]
  return(out)
}

sgd.onlineAlgorithm <- function(t, online.out, datapoint, experiment) {
  at = experiment$learning.rate(t)
  theta.old = onlineOutput.estimate(online.out, t-1)
  st = experiment$score.function(theta.old, datapoint)
  theta.new = theta.old + at * st
  return(theta.new)
}

asgd.onlineAlgorithm <- function(t, online.out, datapoint, experiment) {
  theta.t = sgd.onlineAlgorithm(t, online.out, datapoint, experiment)
  # Average out
  theta.new = (1/t) * theta.t + (1-1/t) * onlineOutput.estimate(online.out, t-1) 
  return(theta.new)
}

implicit <- function(dataset, experiment) {
  
}