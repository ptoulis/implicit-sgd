## Unit-tests for the experiment module.
source("terminology.R")
source("experiment.R")

test.base.learning.rate <- function() {
  rand = rpois(1, lambda=100)
  CHECK_EQ(base.learning.rate(t=rand, gamma0=rand, alpha=rand, c=0),
           rand, msg="If c=0 then rate=gamma0")
  CHECK_EQ(base.learning.rate(t=rand, gamma0=rand, alpha=0, c=1),
           rand, msg="If a=0 then rate=gamma0")
  CHECK_NEAR(base.learning.rate(t=1, gamma0=rand, alpha=1/rand, c=1),
             0.5 * rand, msg="If a *gamma0 * t = 1 then rate=0.5 gamma0")
}

test.glm.score.function <- function() {
  h0 = function(x) { 0 }
  h1 = function(x) { x }
  n = 1 + rpois(1, lambda=10)
  xt = rpois(n, lambda=20)
  yt = sum(xt)
  point = list(xt=xt, yt=yt)
  CHECK_EXCEPTION(glm.score.function(h, theta=c(0), datapoint=point),
                  msg="length(theta) should be = to length(xt)")
  CHECK_EQ(glm.score.function(h0, theta=rep(1, n), datapoint=point),
           yt * xt, msg="If h()=0 then return yt*xt")

  CHECK_NEAR(glm.score.function(h1, theta=rep(1, n), datapoint=point),
             rep(0, n), msg="If h(x)=x then return 0 since yt = sum(xt)")
  CHECK_NEAR(glm.score.function(h1, theta=rep(0, n), datapoint=point),
             yt * xt, msg="If h(x)=x then return 0 since yt = sum(xt)")
}

test.normal.experiment <- function() {
  n = 5 * 1000
  e = normal.experiment(niters=n)
  CHECK_experiment(e)
  CHECK_EQ(e$niters, n)
  d = e$sample.dataset()
  CHECK_EQ(nrow(d$X), n)
  CHECK_EQ(nrow(d$Y), n)
  residuals = d$Y - d$X %*% e$theta.star
  CHECK_TRUE(shapiro.test(residuals)$p.value > 0.05)
  qq = qqnorm(residuals)
  diff = abs(qq$x - qq$y)
  # Less than 10% have difference > 0.1 in qqplot
  CHECK_TRUE(length(which(diff > 0.1)) < 0.1 * n)
  rand.t = sample(1:n, 1)
  point = get.dataset.point(d, rand.t)
  theta = runif(e$p)
  CHECK_EQ(e$score.function(theta, point),
           (point$yt - sum(point$xt * theta)) * point$xt)
}



