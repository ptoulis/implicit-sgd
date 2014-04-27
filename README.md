implicit-glms
=============

Fitting large-scale Generalized Linear Models (GLMs) with implicit updates.

Assume we have a model that produces observations 
    
    y_t  ~  f(theta*)
   
for t=1,2...  theta* = parameter vector in R^p
and y_t are one-dimensional observations, indexed by t.

General methods in online learning (such as stochastic gradient descent)
are using stochastic approximation methods that are of the form:

    theta_t = theta_{t-1} + a_t * S'(yt; theta_{t-1})     (1)

where S'(.) is the Fisher score function (gradient of log-likelihood).

The implicit method is a simple twist in (1) as

  theta_t = theta_{t-1} + a_t * S'(yt; theta_t)     (2)
  
where the Fisher score is evaluated in the future iterate.

Such updates have long been known to possess desirable stability properties in numerical
analysis (e.g. Crank-Nicholson method) and in signal processing (e.g. NLMS filter)
However, they haven't gained popularity because it is usually hard to calculate (2).

It is shown by Toulis et. al. (2014) that for the GLM family (2) can be written as

 theta_t = theta_{t-1} + a_t * (y_t - h(x_t' theta_t)     (2)

and that this is equivalent to

  ξ_t = 




# Assume we have a parametric statistical model with a p x 1 parameter "theta"
#
# ----  Define DATAPOINT to be a list(xt, yt) where xt = px1 vector, yt=1dim value
# Then a DATASET is a list (X, Y) where X=(niters x p) and Y = niters x 1
#
# ----  Create an EXPERIMENT (stochastic approximation). This is comprised by 
#   "theta.star" = # (px1) vector of real values, 
#    p = length(theta.star) = #parameters. That is a terrible name..
#   niters = (#iterations), 
#   sample.dataset() = samples DATASET object
#   score.function = \nabla loglik(theta, data_t) = score function
#   learning.rate = function (t, ...) = gives the learning rate > 0
#   Vx = covariance matrix of xt (covariates/features)
#   Sigma = theoretical covariance matrix for t θt  , t -> infty
#
# ----  An OnlineAlgorithm is defined by a function that accepts 
#   t = no. of iteration
#   onlineOutput = current OnlineOutput object.
#   data.history  = DATASET object from 1:t
#   experiment = EXPERIMENT object, has learning rate/score function
#
# The idea is that the algorithm will use the data up to t, and the current estimates
# to create a new estimate. Usually, it will only need xt, yt, θt, 
# i.e. only the data + estimate at the previous time point.
#
# ----  The online algorithm returns an OnlineOutput object:
#   estimates = (p  x niters) matrix of estimate i.e. (theta_t)
#   last = last vector of estimates
#
# Assume that we run the online algorithm for a specific experiment, k times.
# ----  A MultipleOnlineOutput object is the result of run.online.algorithm.many()
# and it is *not* a LIST of OnlineOutput objects. Rather it is a list
# with the following dimensions:
#
#   {algorithmName}{iteration} = matrix(p x nsamples)
# This holds the output like
#   out[[sgd.onlineAlgorithm]][[t]] = matrix(p x nsamples)
# has all the samples of θt    (nsamples)
#
# ----  A MultipleOnlineOutputParams object (mulOutParams) defines all arguments
# to run many samples run.online.algorithm:
# i.e. it is 
#   {experiment, nsamples, algos}
# 
# ----  A BENCHMARK is a LIST object {mulOut, lowHigh, experiment}:
#   mulOut = MultipleOnlineOutput object (all data)
#   lowHigh = LIST{algoName}{low/high} = [] vector of values
#   experiment = EXPERIMENT that generated the data
#   draw = OPTIONAL draw params
#
# ----  A processParams object defines the data transformation to multipleOnlineOutput.
#   It is a list {vapply, theta.fn} where vapply = {T,F} defines whether
#   we are transforming vectors of theta_t or F if we transform the entire 
#   matrix of theta_t samples.
#   NEW: theta.t.fn is disabled. We made the following hard-coding
#       if vapply = T then theta.fn = default.bias.dist()
#            ""   = F then theta.fn = default.var.dist()
#
# ----  An object BenchmarkParams defines what we need to create a BENCHMARK i.e, 
#   LIST{name, mulOutParams, processParams}
#
# ----  This is the output of generic.benchmark(). 
# A MultipleBenchmark is a LIST of BENCHMARK objects.
