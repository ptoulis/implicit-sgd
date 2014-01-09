## Test for terminology
source("terminology.R")
source("experiment.R")

test.dataset.operations <- function() {
  nsamples = rpois(1, lambda=1000)
  # NOTE: This picks the "normal" experiment by default
  e = get.experiment(niters=nsamples)
  d = e$sample.dataset()
  # 1. Test an exception occurs when index is out of bounds
  CHECK_EXCEPTION(get.dataset.point(d, 2 * nsamples), msg="out of range")
  random.t = sample(1:nsamples, 1)
  print(sprintf("Testing GET.DATASET.POINT: Random iteration =%d", random.t))
  point = get.dataset.point(d, random.t)
  # 2. Check the X, Y values are correct
  CHECK_EQ(point$xt, d$X[random.t, ])
  CHECK_EQ(point$yt, d$Y[random.t, ])
  
  CHECK_EQ(dataset.size(d)$nsamples, nsamples)
}

test.output.operations <- function() {
  nsamples = rpois(1, lambda=1000)
  # NOTE: This picks the "normal" experiment by default
  e = get.experiment(niters=nsamples)
  d = e$sample.dataset()
  out = empty.onlineOutput(dataset=d)
  CHECK_SETEQ(out$estimates[, 10], 0)
  CHECK_EQ(ncol(out$estimates), e$niters)
  CHECK_EQ(out$estimates[, 20], onlineOutput.estimate(out, 20))
  randt = sample(1:nsamples, 1)
  p = length(e$theta.star)
  x = rpois(p, lambda=100)
  out <- add.estimate.onlineOutput(out, randt, x)
  CHECK_EQ(out$estimates[, randt], x)
  for(i in 1:nsamples) {
    out <- add.estimate.onlineOutput(out, i, rep(i, p))
  }
  # Now this is a table like:
  # 1 2 3 4 ...
  # 1 2 3 4 ...
  # 1 2 3 4 ...  etc
  CHECK_EQ(colMeans(out$estimates), 1:nsamples)
  CHECK_EQ(unique(rowSums(out$estimates)), sum(1:nsamples))
}