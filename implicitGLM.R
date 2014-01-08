##    GLMs and implicit updates.

## The ξ equation:  (g() = canonical link)
#  θ = βt * xt
#  ξ  =  η * (y - g^-1(θ + c * ξ))
solve.implicit <- function(beta.old,
                           eta.t, 
                           yt, 
                           xt,
                           g.fn,
                           g.inv.fn,
                           verbose=F) {
    
    ## 1. Compute the implicit Eq. params
    ct = as.numeric(xt %*% t(xt))
    theta.t = as.numeric( xt %*% beta.old)
    rt = eta.t * (yt - g.inv.fn(theta.t))
    h.fn <- function(x) {
             return(eta.t * (yt-g.inv.fn(theta.t+ct * x)))
    }
    ## 2. Define the object function to zero out
    fn.obj <- function(x) x-h.fn(x)
    
    cross.point = (g.fn(yt)-theta.t) / ct
    lower.x = ifelse(rt<0, max(rt, cross.point), 0)
    upper.x = ifelse(rt>0, min(rt, cross.point), 0)
    #print(sprintf("yt=%d ct=%s Lower = %s, upper=%s, rt = %s, cross=%s", yt, ct, lower.x, upper.x, rt,  cross.point))
    if(lower.x == upper.x)
        return(lower.x)
   ###  Something went wrong
    if(lower.x > upper.x) {
        print(sprintf("Eta.t=%f  rt=%f", eta.t, rt))
        warning(sprintf("Something wrong in implicit solving. lower=%f  upper=%f", lower.x, upper.x))
        print(sprintf("eta.t=%f, yt= %d", eta.t, yt))
        print("Beta")
        print(beta.old)
        print("xt")
        print(xt)
        error("A")
    }
    if(verbose) {
        print(sprintf("theta.t =%.1f  ct=%.1f", theta.t, ct))
        print(sprintf("Solving implicit:  rt=%.2f", rt)) 
        print(sprintf("X,Y  bounds [%3f, %3f]", lower.x, upper.x))
    }

    return(uniroot(fn.obj, lower=lower.x, upper=upper.x, tol=1e-20)$root)   
}

##   Implicit updates. 
implicit.update = function(beta.old, t, yt, xt) {
  eta.t = typical.eta.t(t)
  ## Solve implicit equation
  xi.t = solve.implicit(beta.old, 
                        eta.t, 
                        yt, xt, 
                        g.poisson, 
                        g.inv.poisson, 
                        verbose=F )
  
  ##  New beta vector.
  new.beta = beta.old + xi.t * t(xt)
  
  return(new.beta)
}



g.inv.logistic <- function(x) {
    return(exp(x) / (1+exp(x)))
}
g.poisson = function(x) log(x)
g.inv.poisson <- function(x) exp(x)
