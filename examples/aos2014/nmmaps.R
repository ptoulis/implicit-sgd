# Copyright (c) 2013
# Panos Toulis, ptoulis@fas.harvard.edu
#
# Using Implicit SGD to fit the NMMAPS data

DATA_DIR = "~/A/data/nmmaps/"
data.files <- list.files(path=DATA_DIR, full.names=T, pattern="bz2")

###############################################################################
## Fit the core NMMAPS model with PM10 and mortality
## Copyright (C) 2004, Roger D. Peng <rpeng@jhsph.edu> and Aidan McDermott
##     
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
###############################################################################
 
## NOTE: Before the function `fitSingleCity' can be used, the data
## must be preprocessed using the `basicNMMAPS' function in the
## `NMMAPSdata' package.  Type `?basicNMMAPS' after loading the
## `NMMAPSdata' package into R for more information.  `fitSingleCity'
## will not work correctly with the full/raw database.
rm(list=ls())
library(ggplot2)
library(splines)

merge.all.data <- function() {
  data.nmmaps = NULL
  files = list.files("~/A/data/nmmaps/", pattern="bz2", full.names=T)
  print(sprintf("There are %d files.", length(files)))
  for(i in 1:length(files)) {
    file = files[i]
    m = regexec(pattern="(\\w+)\\.rda", text=file)
    df.name = regmatches(file, m)[[1]][2]
    print(sprintf("Loading object %s -- Total files %d/%d", df.name, i, length(files)))
    # 1. Load the .rda file for a specific city
    load(file)
    # 2a. Fix the dataset
    df = get(df.name)
    df = fixData(df)
    if(!is.null(df)) {
      # 2b. Bind to big data file if it is ok.
      data.nmmaps = rbind(data.nmmaps, df)
    } else {
      print(sprintf("Dataset %s was NOT added..", file))
    }
    # 3. Clear some memory
    rm(df)
    rm(df.name)
    gc()
    print(object.size(data.nmmaps), units="Mb")
  } 
  return(data.nmmaps)
}

as.design <- function(df) {
  X = matrix(0, nrow=nrow(df), ncol=0)
  #days = seq(2, 7)
  #X = cbind(X, matrix(sapply(df$dow, function(i) days==i), nrow=length(df$dow), byrow=T))
  #X = cbind(X, df$agecat==1)
  #X = cbind(X, df$agecat==2)
  cities = unique(df$city)
  print(cities)
  save(cities, file="cities.rda")
  pollute = df$o3median * df$tmpd
  print("> Fitting cubic spline to pollute var...")
  ns.pollute = ns(pollute, 100)
  # Monitor progress
  print("> Binding time cubic spline (dimension 400)...")
  X = cbind(X, ns(df$time, 400))
  print("> Binding Age 1 cubic spline (dimension 100)...")
  X = cbind(X, ns.pollute * (1-df$Age2Ind) * (1-df$Age3Ind))
  print("> Binding Age 2 cubic spline (dimension 100)...")
  X = cbind(X, ns.pollute * df$Age2Ind)
  print("> Binding Age 3 cubic spline (dimension 100)...")
  X = cbind(X, ns.pollute * df$Age3Ind)
  print(sprintf("> Binding Cities (dimension %d)...", length(cities)))
  X = cbind(X, matrix(sapply(df$city, function(i) as.numeric(cities==i)), 
                      nrow=length(df$city), byrow=T))
  X = cbind(X, df$death)
  print("> Setting column names..")
  cols <-   c(paste("t", 1:400, sep=""),
              paste("age1_pollute", 1:100, sep=""),
              paste("age2_pollute", 1:100, sep=""),
              paste("age3_pollute", 1:100, sep=""),
              paste("city", 1:length(cities), sep=""),
              "death")
  colnames(X) = cols
  return(X)
}
fixData <- function(dataframe) {
  if (all(is.na(dataframe[, "pm10tmean"])))
    return(NULL)
  is.na(dataframe[, "death"]) <- as.logical(dataframe[, "markdeath"])
  is.na(dataframe[, "cvd"]) <- as.logical(dataframe[, "markcvd"])
  is.na(dataframe[, "resp"]) <- as.logical(dataframe[, "markresp"])
  Age2Ind <- as.numeric(dataframe[, "agecat"] == 2)
  Age3Ind <- as.numeric(dataframe[, "agecat"] == 3)
  dataframe[, "dow"] <- as.factor(dataframe[, "dow"])
  dataframe[, "agecat"] <- as.factor(dataframe[, "agecat"])
  #   varList <- c("cvd", "death", "resp", "tmpd", "rmtmpd", "dptp",
  #                "rmdptp", "time", "agecat", "dow", "pm10tmean", "o3mean",
  #                paste("l", 1:7, "pm10tmean", sep = ""))
  varList <- c("city", "death", "tmpd", "o3median", "time")
  df = data.frame(dataframe[, varList], Age2Ind = Age2Ind, Age3Ind = Age3Ind)
  ma <- function(x, n=4) { 
    filter(x, rep(1/n,n), sides=1)
    x[1:(n-1)] <- x[n]
  }
  filter.vars = c("o3mean", "tmpd")
  # simple imputation2
  for(f in names(df)) {
    y = df[[f]]
    x = which(is.na(y))
    if(length(x) > 0) {
      # 
      df[[f]][x] <- mean(y, na.rm=T)
      # if temperature or O3, take a running average.
      if(f %in% filter.vars) {
        z = df[[f]]
        df[[f]] <- (z-mean(z))/sd(z)
      }
    }
  }
  if(sum(is.na(df))) {
    warning("Some columns still have NA values..")
    return(NULL)
  }
  return(df)
}


fit.all.implicit <- function(X) {
  n = nrow(X)
  p = ncol(X) - 1 # last column is response.
  cov.i = 1:p
  theta.hat = rep(0, p)
  a.optimal = 1800.0
  solve.implicit <- function(y, theta, a, x, x.norm, B) {
    if(B[1] == B[2]) return(B[1])
    eta = sum(theta * x)
    f = function(z) {
      z - a * (y - exp(eta + x.norm * z))
    }
    uniroot(f, interval=B)$root
  }
  pb = txtProgressBar(style=3)
  for(i in 1:n) {
    xi = as.numeric(X[i, cov.i])
    yi = as.numeric(X[i, p+1])
    yi.pred = exp(sum(theta.hat * xi))
    ai = 1 / (1 + (1/a.optimal) * i)
    
    ri = ai * (yi - yi.pred)
    Bi = c(0, ri)

    if(ri <= 0)
      Bi = c(ri, 0)
    
    ksi = solve.implicit(yi, theta.hat, ai, xi, sum(xi^2), Bi)
    theta.hat = theta.hat + ksi * xi
    setTxtProgressBar(pb, i/n)
  }
  return(theta.hat)
}
mse <- function(x, y) sqrt(mean((x-y)^2))
run.experiment.nmmaps <- function(data.size="small") {
  print(sprintf("Loading data...(size=%s)", data.size))
  if(data.size=="small") {
    print("> Running small dataset..")
    load("datasets/X.small.rda")
    print(sprintf("> Design matrix X: %d obs. %d covariates. Memory Size=%.1f Mb", 
                  nrow(X.small), ncol(X.small), object.size(X.small) / (1024 * 1024)))
    
    df = as.data.frame(X.small)
    df$death = as.integer(df$death)
    print("> Running glm()...")
    fit = glm(death ~ . + 0, family=poisson, data=df)
    glm.estimates = as.numeric(coef(fit))
    print("> Running implicit...")
    sgd.estimates = fit.all.implicit(X.small)
    print(sprintf("MSE = %.3f", mse(sgd.estimates, glm.estimates)))
    df = data.frame(glm=glm.estimates, sgd=sgd.estimates)
    g = ggplot() + geom_point(data=df, aes(x=glm, y=sgd)) + 
      geom_line(data=df, aes(x=glm, y=glm), lty=2)
    g
    #  + geom_line(aes(x=glm.estimates, y=glm_estimates))
    # lines(b1, b1, col="red", lty=3, lwd=0.5) 
  } else if(data.size=="big") {
    load("datasets/X.all.rda")
    print(sprintf("> Design matrix X: %d obs. %d covariates. Memory Size=%.1f Mb", 
                  nrow(X.all), ncol(X.all), object.size(X.all) / (1024 * 1024)))
    print("> Fit using the implicit sgd method.")
    f = system.time({ betas = fit.all.implicit(X.all) })
    print(sprintf("Elapsed time %.2f seconds..", f[["elapsed"]]))
  } else {
    stop("Data size should be in {small, big}")
  }
}
