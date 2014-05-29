# Introduction

## Research paper
This is the accompanying code implementation of the methods and algorithms 
presented
```
Panos Toulis, Jason Rennie, Edoardo Airoldi, 
"Statistical analysis of stochastic gradient methods for generalized linear models", 
ICML, Beijing, China, 2014.
```

We will refer to this paper as (Toulis et. al., 2014) in what follows.

## Model overview
Assume we have a model that produces observations 
```
    y_t  ~  f(theta*)
```   
for t = 1,2...  theta* = parameter vector in R^p
and y_t are one-dimensional observations, indexed by t.

General methods in online learning (such as stochastic gradient descent)
are using stochastic approximation methods that are of the form:
```
    theta_t = theta_{t-1} + a_t * S'(yt; theta_{t-1})     (1)
```

where S'(.) is the Fisher score function (gradient of log-likelihood).
The implicit method is a simple twist in (1) as
```
    theta_t = theta_{t-1} + a_t * S'(yt; theta_t)     (2)
```
where the Fisher score is evaluated in the future iterate.

Such updates have long been known to possess desirable stability properties in numerical
analysis (e.g. Crank-Nicholson method) and in signal processing (e.g. NLMS filter)
However, they haven't gained popularity because it is usually hard to calculate (2).

It is shown by Toulis et. al. (2014) that for the GLM family (2) can be written as
```
 theta_t = theta_{t-1} + a_t * (y_t - h(x_t' theta_t)     (3)
```
and that this is equivalent to
```
  ξ_t = g(ξ_t)   (4)
```

where ξ_t is a real variable. In other words, the implicit update in (3)
is reduced to a one-dimensional implicit equation (4).
Furthermore, solving (4) can be computationally easy because efficient search
bounds are obtainable.

## Main results

The main results in (Toulis et. al., 2014) can be summarized as follows:
* Implicit methods can be efficiently applied in the GLM family using Equation (4).

* SGD and Implicit methods (Eq. (2) + (3)) in the GLM family are asymptotically equivalent.
  Exact formulas for the bias and variance can be derived for both methods, showing that 
  the two methods have the same asymptotic bias and efficiency.

* The Implicit method is more biased in small samples, however it exhibits 
 smaller empirical variance of the estimates. This is because the Implicit method
  is "unconditionally stable", whereas the typical SGD is not.

# Code Terminology

There are several terms/keywords defined in this project.
To clarify, we use the terms as follows:

* ```DATAPOINT``` is defined as a ```LIST(xt, yt)``` where xt = (p x 1) vector, 
  and yt is real. Then a ```DATASET``` is a ```LIST(X, Y)``` where 
   X = (niters x p) and Y = (niters x 1).

* An ```EXPERIMENT``` is a list comprised by 
   * **theta.star** =  (p x 1) vector of real values, 
    * **p** = length(```theta.star```) = #parameters. (that is a terrible name..)
    * **niters** = #iterations in the experiment
    * **sample.dataset()** = samples ```DATASET``` object
    * **score.function()** = \nabla loglik(theta, data_t) = px1 vector
    * **learning.rate()** = function (t, ...) = gives the learning rate at t > 0
    * **Vx** = covariance matrix of xt (covariates/features) i.e. Vx = Cov(xt)
    * **J** = Fisher information = E(h'(θ'x) x x')

* An ```OnlineAlgorithm``` is a function with the following arguments 
   * **t** = no. of iteration
   * **onlineOutput** = current ```OnlineOutput``` object.
   * **data.history**  = ```DATASET``` object from 1:t.
   * **experiment** = EXPERIMENT object, has learning ratescore function

 The idea is that the algorithm will use the data up to t, and the current estimates
 to create a new estimate. Usually, it will only need xt, yt, θt, 
 i.e. only the data + estimate at the previous time point.

* The ```OnlineOutput``` is the output of an ```OnlineAlgorithm```.
  This is comprised by
    * **estimates** = (p  x niters) matrix of estimate i.e. (theta_t)
    * **last** = last vector of estimates (length p)

 Assume that we run the online algorithm for a specific experiment, k times.

* A ```MultipleOnlineOutput``` object is the result of ```run.online.algorithm.many()```
 and it is a ```LIST{algorithmName}{iteration}``` = matrix(p x nsamples)
 For example, something like
   ```out[[sgd.onlineAlgorithm]][[t]]``` = matrix(p x nsamples)
 having all the samples of θt    (nsamples)

* A ```MultipleOnlineOutputParam``` object (mulOutParams) defines all arguments
 to run many samples from ```run.online.algorithm```:
 i.e. it is a ```LIST{experiment, nsamples, algos}```.
 
* A ```BENCHMARK``` is a ```LIST{mulOut, lowHigh, experiment}``` where 
    * **mulOut** = ```MultipleOnlineOutput``` object (all data)
    * **lowHigh** = ```LIST{algoName}{lowhigh}``` = [] vector of values
    * **experiment** = ```EXPERIMENT``` that generated the data
    * **draw** = OPTIONAL drawing params

* A ```processParams``` object defines the data transformation to ```multipleOnlineOutput```.
   It is a ```LIST{vapply, theta.fn}``` where vapply = {T,F} defines whether
   we are transforming vectors of theta_t or F if we transform the entire 
   matrix of theta_t samples.
   NEW: theta.t.fn is disabled. We made the following hard-coding
       if vapply = T then theta.fn = default.bias.dist()
            ""   = F then theta.fn = default.var.dist()

* An object ```BenchmarkParams``` defines what we need to create a ```BENCHMARK``` i.e, 
   ```LIST{name, mulOutParams, processParams}```.
