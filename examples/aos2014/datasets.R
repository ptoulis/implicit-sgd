# Copyright (c) 2013
# Panos Toulis, ptoulis@fas.harvard.edu
#
dataset.filename <- function(dim.p, dim.n, beta.file=F, is.big=F) {
  return(sprintf("datasets/%sDataset-p%1.1f-n%1.1f%s.csv", 
                 ifelse(is.big, "big/", ""),
                 log(dim.p, 10), log(dim.n, 10),
                 ifelse(beta.file, "-BetaFile", "")))
}

sample.beta <- function(dim.p) {
  # 1. Samples a beta vector (px1) of model parameters.
  # The model parameters are either "high", "medium" or "low" (levels)
  beta.levels = c(-1, -0.35, -0.001, 0.001, 0.35, 1)
  # How frequent are the levels
  beta.level.freq = c(0.05, 0.1, 0.25, 0.3, 0.2, 0.1)    
  beta = c()
  for(i in 1:length(beta.levels)) {
    beta <- c(beta, 
              rep(beta.levels[i], as.integer(dim.p * beta.level.freq[i])))
  }
  while(length(beta) != dim.p) {
    warning("length(beta) not equal to p. Truncating...")
    if(length(beta) > dim.p) {
      beta = beta[1:dim.p]
    } else {
      beta = c(beta, sample(beta.levels, size=dim.p-length(beta), replace=F))
    }
  }
  return(sample(beta))
}

# 2. Sample the covariates (nxp) sparse matrix.
sample.X <- function(dim.p, dim.n, verbose=F) {
  # Create a binary feature matrix
  L = dim.p * dim.n 
  if(F) {
    # TODO(ptoulis): Add support for sparse matrices.
    nnzeros = min(0.05 * L, L**0.7)
    nnzero.pos = sample(L, size=nnzeros, replace=F)
    nnzero.i = 1 + as.integer(nnzero.pos / dim.p)
    nnzero.j = nnzero.pos - (nnzero.i-1) * dim.p
    k =  which(nnzero.j==0)
    nnzero.i[k] <- nnzero.i[k] - 1
    nnzero.j[k] <- dim.p
    return(sparseMatrix(i=nnzero.i, j=nnzero.j, dims=c(dim.n, dim.p)))
  } else {
    # 1. Pick non-zeros (and positions in the matrix)
    nnzeros = 0.08 * L
    nnzero.pos = sample(L, size=nnzeros, replace=F)
    
    X = matrix(0, nrow=dim.n, ncol=dim.p)
    X[nnzero.pos] <- 1
    X[, 1] <- 1
    return(matrix(as.numeric(X), nrow=dim.n, ncol=dim.p))
  }
}

create.dataset <- function(dim.p=1e2, dim.n=1e3, is.big=F) {
  filename = dataset.filename(dim.p, dim.n, is.big=is.big)
  beta.filename = dataset.filename(dim.p, dim.n, beta.file=T, is.big=is.big)
  chunk = as.integer(dim.n / 50)
  beta = sample.beta(dim.p)
  write.table(matrix(beta, nrow=1), file=beta.filename, sep=",", 
              row.names=F, col.names=sapply(1:dim.p, function(i) sprintf("b%d", i)))
  print("> Saved beta vector to file..")
  cnames = sapply(1:dim.p, function(i) sprintf("X%d", i))
  cnames = c(cnames, "Y")
  
  imax = as.integer(dim.n / chunk)
  print(sprintf("> Need to iterate save-file %d times", imax))
  pb = txtProgressBar(style=3)
  sigma = 1
  rows.saved = 0
  for(i in 1:imax) {
    if(i==imax) {
      chunk = dim.n - rows.saved
    }
    X = sample.X(dim.p=dim.p, dim.n=chunk)
    y = round(X %*% beta + rnorm(chunk, sd=sigma), 3)
    X = cbind(X, y)
    colnames(X) <- cnames
    if(i==1) {
      # save the betas as the second line.
      write.table(X, file=filename, row.names=F, col.names=T, sep=",")
    } else {
      write.table(X, file=filename, row.names=F, col.names=F, append=T, sep=",")
    }
    setTxtProgressBar(pb, value=i/imax)
    rows.saved = rows.saved + chunk
  }
}

remove.dataset <- function(dim.p, dim.n) {
  filename = dataset.filename(dim.p, dim.n)
  if(file.exists(filename)) {
    unlink(filename)
    beta.file = dataset.filename(dim.p, dim.n, beta.file=T)
    unlink(beta.file)
  }
}

vector.dist <- function(x, y) {
  if(length(x) != length(y))
    stop("Vectors should have equal length to calculate distance.")
  # Compute MSE.
  sqrt(sum((x-y)^2))
}

CHECK <- function(claim, msg) {
  if(!claim) 
    stop(msg)
}

load.dataset <- function(dim.p, dim.n) {
  filename = dataset.filename(dim.p, dim.n)
  print(sprintf("> Loading dataset..%s.", filename))
  df = read.csv(filename)
  print(sprintf("> Dataset loaded. no.obs=%d no.covariates=%d",
                nrow(df), ncol(df)))
  return(df)
}

load.beta <- function(dim.p, dim.n, is.big=F) {
  beta.filename = dataset.filename(dim.p, dim.n, beta.file=T, is.big=is.big)
  print(sprintf("> Loading ground-truth (beta vector) ..%s.", beta.filename))
  return(scan(beta.filename, skip=1, sep=",", nlines=1, quiet=T))
}
