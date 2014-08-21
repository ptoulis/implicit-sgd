## Recursion
rm(list=ls())

init.matrix <- function(p) {
  # Creates the matrices to run the recursions.
  # I = identity
  # B = the iteration matrix (I - a_n * B) ...
  # V = the "noise" 
  I <- diag(p)
  R1 = matrix(runif(p^2), nrow=p, ncol=p)
  R2 = matrix(runif(p^2), nrow=p, ncol=p)
  
  B <-  I + R1 %*% t(R1)
  V  <- R2 %*% t(R2) 
  return(list(B=B, I=I, V=V))
}

best.learning.rate <- function(M) {
  lam = eigen(M)$values
  optim(par=c(0), lower = 1/min(lam)+1e-5, method="L-BFGS-B",
        fn = function(x) sum(x^2 / (x * lam - 1)))$par
}

recursion.first.order <- function(D=c(), alpha=0.2, constant.rate=F, 
                                  niters=1000,
                                  implicit=F) {
  if(length(D)==0) {
    warning("Data were empty. Creating default with p=10")
    D = init.matrix(p=10)
  }
  # Runs the recursion
  #   Xn = (I - an B) * Xn-1 + an * V
  B = D$B
  I = D$I
  V = D$V
  p = nrow(B)
  
  Y.lim = solve(B) %*% V  # limit
  
  if(any(eigen(Y.lim)$values < 0) | any(eigen(V)$values < 0)) {
    stop("Eigenvalues of B cannot be negative.")
  }
  if(any(eigen(I - alpha * B)$values < 0)) {
    warning(sprintf("Likely to be unstable. Try values for alpha < %.3f", 2 / max(eigen(B)$values)))
  }
  # Current iterate.
  Xi = matrix(0, nrow=p, ncol=p)
  # residuals.
  ei = c()
  
  for(i in 1:niters) {
    # 1. Compute learning rate
    ai  = alpha / (alpha + i)
    if(constant.rate) {
      ai = alpha
    }
    if(implicit) {
      Xi = solve(I + ai * B) %*% (Xi + ai * V)
    } else {
      Xi = (I - ai * B) %*% Xi + ai * V 
    }
    Yi = Xi
    ei = c(ei, norm(Yi - Y.lim, type = "F"))
  }
  if(p < 5) {
    print("Empirical limit")
    print(Xi)
    # print(Xi)
    print("Theoretical limit")
    print(Y.lim)
  } 
  plot(ei, type="l", ylim=c(-0.5, 10))
  print("Summary")
  print(summary(tail(ei, 100)))
}

recursion.second.order <- function(D, a=0.2, 
                                  niters=1000,
                                  implicit=F) {
  B = D$B
  I = D$I
  V = D$V
  p = nrow(B)
  
  Y.lim = solve(B - I/a) %*% V  # limit
  
  if(any(eigen(Y.lim)$values < 0) | any(eigen(V)$values < 0)) {
    stop("Eigenvalues of B cannot be negative.")
  }
  
  Xi = matrix(0, nrow=p, ncol=p)
  ei = c()
  plot.points = as.integer(seq(1, niters, length.out=niters/10000))
  
  for(i in 1:niters) {
    ai  =  a#a / i
    if(implicit) {
      Xi = solve(I + ai * B) %*% (Xi + ai^2 * V)
    } else {
      Xi = (I - ai * B) %*% Xi + ai^2 * V 
    }
    Yi = (1/ai) * Xi
    ei = c(ei, norm(Yi - Y.lim, type = "F"))
    if(i %in% plot.points) {
      print(sprintf("i = %d/%d", i, niters))
      ei.sub = ei[100:length(ei)]
      last.ei <- ei[100:length(ei)]
      plot(last.ei, type="l")
      print(summary(last.ei))
    }
  }
  print("Empirical limit")
  print(Yi)
  # print(Xi)
  print("Theoretical limit")
  print(Y.lim)
  print("Summary")
  print(summary(tail(ei, 100)))
}
