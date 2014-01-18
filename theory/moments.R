# Copyright (c) 2013
# Panos Toulis, ptoulis@fas.harvard.edu
#
# Some simulations to back the theory.

#
variance.recursion <- function(experiment) {
  
}

test.product.series <- function() {
  alpha.vector = seq(0, 1, length.out=20)
  y = sapply(alpha.vector, product.series)
  i = sample(length(alpha.vector), 1)
  ralpha = alpha.vector[i]
  shouldBe = sapply(1:10^6, function(t) (1 - (ralpha / t)^2))
  CHECK_NEAR(prod(shouldBe), y[i], tol=1e-4)
}

product.series <- function(alpha, type, tol=1e-4) {
  # This is equal to sin (π α) / πα (Euler formula)
  maxT = 10^8
  entire.interval = 1:maxT
  chops = 5000
  current.interval = seq(1, chops)
  run.loop = T
  y.old = 1
  while(run.loop) {
    y.vector = NA
    if(type=="sgd") {
      y.vector = sapply(current.interval, function(t) (1 - alpha / t))
    } else {
      y.vector = sapply(current.interval, function(t) 1 / (1 + (alpha / t)))
    }
    y.new = prod(y.vector) * y.old
    diff = abs(y.new - y.old)

    if(diff > tol) {
      current.interval = seq(head(current.interval,1 ) + chops,
                             tail(current.interval, 1) + chops)
      y.old = y.new
    } else {
      run.loop = F
    }
  }
  print(sprintf("alpha=%.3f Last point was %d", alpha, tail(current.interval, 1)))
  return(y.old)
}

plot.bias.asymptotic <- function() {
  a = seq(0.05, 1, length.out=213)
  y = sin(pi * a) / (pi * a)
  plot(a, y, "l")
}