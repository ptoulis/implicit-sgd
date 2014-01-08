rm(list=ls())
N = 10000

xs = runif(N,min=0.1,max=1)
beta.true = 22.5
eta.t = function(t) {
    10^6/(t+1)
}
sample.y = function(t, epsilon=0) 
  {
   s1 = rnorm(1, mean=xs[t] * beta.true, sd=2)
   s2 = rnorm(1, mean=200 * xs[t], sd=2)
   Ind = rbinom(1, size=1, prob=1-epsilon)
   return(Ind * s1 + (1-Ind) * s2)
 }
## SGD once.
both = function(epsilon=0, plot=F) {
  betas.sgd = rep(1, N)
  betas.impl = rep(1,N)
  ys = rep(0,N)
  ds = rep(0,N)
    for(t in 1:(N-1)) {
        yt = sample.y(t, epsilon=epsilon)
        etat = eta.t(t)
        ys[t] = yt
        ## sgd
        betas.sgd[t+1] = betas.sgd[t] + etat * (yt-betas.sgd[t] * xs[t])
        
        ## implicit
        xi.t = (yt-xs[t] * betas.impl[t]) * etat / (1+etat * xs[t]**2)
        betas.impl[t+1] = betas.impl[t] +xi.t * xs[t]
        ds[t] = betas.impl[t+1] - yt/xs[t]
        
      }
  start = as.integer(0.1 * N)
  x = list(sgd=betas.sgd[start:N],  impl = betas.impl[start:N], 
           ys=ys[start:(N-1)],
           ds = ds[start:N])
    if(plot){
        plot.both(x)
    }
    return(x )
}
plot.both = function(both.obj) {
  betas.sgd = both.obj$sgd
  betas.impl = both.obj$impl
  m = max(c(min(betas.sgd,na.rm=T), min(betas.impl)))
  M = min(c(max(betas.sgd, na.rm=T), max(betas.impl)))
  print(c(m,M))
  xaxis = 1:length(betas.sgd)
  plot(xaxis, betas.sgd, type="l", ylim=c(min(m, beta.true-1), 
                                          max(M, beta.true+1)))
  points(xaxis, betas.impl, type="l", col="green")
  points(xaxis, rep(beta.true, length(xaxis)), col="red", cex=0.2)
}

###   Poisson var.recur
sgd.var.recur =function(beta, u.t, sigma, eta.t, x.t ) {
  v1 = u.t^2 
  v2 = eta.t^2 * sigma^2
  
  a = exp(beta * x.t)
  b = exp(x.t^2 * u.t^2)
  
  v3 = eta.t^2 * a^2 * (b-1) * b
  v4 = -2 * eta.t * x.t * u.t^2 * a * sqrt(b)
  
  return(sum(c(v1,v2,v3,v4)))
}
impl.var.recur = function(beta, u.t, sigma, eta.t, x.t, verbose=F) {
    v1 = function(u) u^2
    v3= function(u) {
        a = exp(beta * x.t)
        b = exp(x.t^2 * u^2)
        
        eta.t^2 * a^2 * (b-1) * b
    }
    v4 = function(u) {
        a = exp(beta * x.t)
        b = exp(x.t^2 * u^2)
        2 * eta.t * x.t * u^2 * a * sqrt(b)
    }
   
    
    D = u.t^2 + eta.t^2 * sigma^2
    gn = function(u) v1(u)+v3(u)+v4(u)
    fn = function(u) gn(u)-D
    ## Now solve the implicit equation.
    sD =sqrt(D)
    ux = seq(0,sD, length.out=100)
    plot(ux, sapply(ux, fn), type="l")
    x = uniroot(fn,lower=0, upper=sD,tol=1e-8)
    if(verbose)
        print(sprintf("D=%.3f  Root=%.3f gval=%.5f", D, x$root, gn(x$root)))
    return(x$root)
}

##   Runs the variance recurrence formula.
## With these params, σ = 44   is the critical value (above that it cracks)
run.sgd.var.recur= function(runs=10, 
                        beta=1.0, 
                        u.t=0.5, 
                        sigma=10, 
                        x.t=1.0, 
                        etat= function(t){0.1/(t+1)}) {
    us = c(u.t, rep(NA, runs))
   
    
    for(i in 2:(runs+1)) {
        us[i] = sqrt(sum(sgd.var.recur(beta, us[i-1], sigma, eta.t=etat(i), x.t)))
    }
    return(us)
}
run.impl.var.recur= function(runs=10, 
                        beta=1.0, 
                        u.t=0.5, 
                        sigma=10, 
                        x.t=1.0, 
                        etat= function(t){0.1/(t+1)}) {
    us = c(u.t, rep(NA, runs))
    
    
    for(i in 2:(runs+1)) {
        us[i] = impl.var.recur(beta, us[i-1], sigma, eta.t=etat(i), x.t)
    }
    return(us)
}

find.crackpoint = function(method="sgd", 
                           eta.fn=function(t){0.1/(t+1)}) {
    sigma=1
    isok =T
    while(isok && sigma < 50) {
        ## Which method to run?
        x = NA
        if(method=="sgd") {
            x =run.sgd.var.recur(runs=100, sigma=sigma, etat=eta.fn)
        } else {
            x =run.impl.var.recur(runs=100, sigma=sigma, etat=eta.fn)
            sigma = sigma*2 ## do this because the method is slow.
        }
        isok = !is.na(tail(x,1))
        sigma = sigma+2 * runif(1)
        if(!isok) print("Cracked")
    }
    return(sigma)
}
## Given  ht ~ n^-a    find σ  that cracks the recurrence.
divergence.plot = function(method="sgd") {
    as = seq(0.2, 1, length.out=10) 
    cracks = rep(0, length(as))
    for (i in 1:length(as)) {
        a = as[i]
        eta.fn = function(t) 0.2*(t+1)^{-a}
        cracks[i]=find.crackpoint(method, eta.fn)
        print(sprintf("Done with a=%.1f crack=%.2f", a, cracks[i]))
        ## now sigma is the break point.
    }
    plot(as, cracks, main="Variance recurrence breakdown point", 
         xlab="Learning rate power a (eta_t = n^-a)",
         ylab="Outcome (y) variance breakdown",
         type="l")
    return(list(as=as, cs=cracks))
}
example.poisson.recur = function() {
  beta = 1.3
  u.t  = 1.12
  x.t = 1.0
  eta.t = 0.1
  N = 10^5
  beta.t = rnorm(N, mean=beta, sd=u.t)
  
  
  y.t = rnorm(N, sd=100)#rpois(100000, exp(beta * x.t))
  sigma = sd(y.t)
  
  v1 = var(beta.t)
  v2 = var(eta.t * y.t)
  v3 = var(eta.t * exp(beta.t * x.t))
  v4 = -2 * eta.t * cov(beta.t, exp(beta.t * x.t))
  
  vtheor =poisson.var.recur(beta=beta, 
                            u.t = u.t, 
                               sigma=sigma,
                               eta.t=eta.t, 
                               x.t = x.t)
  print("Empirical values")
  print(sprintf("V(b.t)=%.3f  V(ht yt)=%.3f  V(ht e^bt xt)=%.3f Cov(bt, exp)=%.3f",
                v1, v2, v3, v4))
  print("Theoretical values")
  print(sprintf("V(b.t)=%.3f  V(ht yt)=%.3f  V(ht e^bt xt)=%.3f Cov(bt, exp)=%.3f",
                vtheor[1], vtheor[2], vtheor[3], vtheor[4]))
  
  print(sprintf("Old var value %.3f New var value %.3f", u.t, sum(vtheor)))
  
}

