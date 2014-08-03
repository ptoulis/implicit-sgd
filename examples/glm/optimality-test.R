# Can SGD be made optimal?
rm(list=ls())
theta.star = NA#

sample.x <- function(n, p, rho=0.2) {
  #    generate x's multivariate normal with equal corr rho
  # Xi = b Z + Wi, and Z, Wi are independent normal.
  # Then Var(Xi) = b^2 + 1
  #  Cov(Xi, Xj) = b^2  and so cor(Xi, Xj) = b^2 / (1+b^2) = rho
  z=rnorm(n)
  if(abs(rho)<1){
    beta=sqrt(rho/(1-rho))
    x0=matrix(rnorm(n*p),ncol=p)
    A = matrix(z, nrow=n, ncol=p, byrow=F)
    x= beta * A + x0
  }
  if(abs(rho)==1){ x=matrix(z,nrow=n,ncol=p,byrow=F)}
  
  return(x)
}

get.J <- function(p, rho=0.2) {
  beta=sqrt(rho/(1-rho))
  J = diag(p) + matrix(beta^2, nrow=p, ncol=p)
  # print(J)
  return(J)
}

best.alpha <- function(J) {
  lambdas = eigen(J)$values
  f = function(a) sum(a^2 * lambdas / (2 * a * lambdas-1))
  optim(par=c(0), fn=f, method="L-BFGS-B", lower=1/(2 * min(lambdas)) + 1e-4)$par
}

run.optimal.experiment <- function(N=1000, p=5, nsamples=100) {
  
  theta.star = sample(c(-1.2, -0.5, 0.05, 0.15, 0.75, 1.5), size=p, replace=T)
  p = length(theta.star)
  
  mle.est = matrix(NA, nrow=0, ncol=p)
  sgd.est.opt = matrix(NA, nrow=0, ncol=p)
  sgd.est = matrix(NA, nrow=0, ncol=p)
  
  for(i in 1:nsamples) {
    X = sample.x(N, p=length(theta.star))
    preds = X %*% theta.star
    J = get.J(p)
    J.inv = solve(J)
    best.a = best.alpha(J) - 0.0
    y = preds + rnorm(N)
    
    # update mle
    fit <- lm(y~X + 0)
    theta.mle.vector = coef(fit)
    mle.est = rbind(mle.est, theta.mle.vector)
    
    # Run SGD
    theta.sgd.opt = matrix(0, nrow=p, ncol=1) # "optimal" SGD
    theta.sgd     =  matrix(0, nrow=p, ncol=1) # other SGD
    for(j in 1:length(y)) {
      xj = matrix(X[j, ], ncol=1)
      yj = y[j]
      aj = 1/(1 + (1/best.a) * j)
      
      # Do classical SGD
      y.pred = sum(xj * theta.sgd)
      Delta = (yj - y.pred) * xj
      theta.sgd = theta.sgd + aj * Delta
      rm(list=c("y.pred", "Delta", "aj")) # just to be safe
      
      # Do the optimal now.
      y.pred.opt = sum(xj * theta.sgd.opt)
      Delta.opt = (yj - y.pred.opt) * xj
      theta.sgd.opt = theta.sgd.opt + min(0.1, 1/j) * J.inv %*% Delta.opt
    }
    
    sgd.est.opt = rbind(sgd.est.opt, as.numeric(theta.sgd.opt))
    sgd.est = rbind(sgd.est, as.numeric(theta.sgd))
  }
  rownames(mle.est) <- NULL
  V = solve(J)
  
  print.est <- function(M, str) {
    cat("\n\n=--==-=-=---===-=-\n")
    print(str)
    C = N * cov(M)
    if(nrow(C) < 1)
      print(C)
    theta.est = colMeans(M)
    print(sprintf(">> Estimated theta=%s", paste(round(theta.est, 2), collapse=" ")))
    
    cat(sprintf(">> Trace=%.3f  L2=%.5f", sum(diag(C)), sqrt(mean((theta.est-theta.star)^2))))
  }
  print(sprintf(">> - Real theta=%s", paste(theta.star, collapse=" ")))
  
  print.est(V, "Theoretical variance.")
  print.est(mle.est, "MLE estimate.")
  print.est(sgd.est, "SGD optimal rate.")
  print.est(sgd.est.opt, "SGD w/ optimal variance.")
}


