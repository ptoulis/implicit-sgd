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


library(splines)

fitSingleCity <- function(data, pollutant = "l1pm10tmean", cause = "death",
                          dfyr.Time = 7, pdfyr.time = 0.15, df.Temp = 6,
                          df.Dew = 3, extractors = NULL) {
  ## Argument checking
  stopifnot(is.character(pollutant), is.character(cause), length(cause) == 1)
  
  ## Modify degrees of freedom based on missingness of data
  sub <- data[, c("time", "agecat", "dow", "tmpd", "rmtmpd", "dptp",
                  "rmdptp", cause, pollutant)]
  subset <- complete.cases(sub)
  df.Time <- round( numdf(subset, dfyr.Time) )
  df.time <- round( df.Time * pdfyr.time )
  
  ## Don't setup smooth function of time where there are incomplete cases
  is.na(data[, "time"]) <- !subset
  
  modelFormula <- setupFormula(cause, pollutant, df.Time, df.time,
                               df.Temp, df.Dew)
  print(modelFormula)
  ## Fit the model!
  fit <- glm(modelFormula, family = quasipoisson, data = data,
             control = glm.control(epsilon = 1e-10, maxit = 1000),
             na.action = na.omit)
  
  ## Extract information from the fitted glm model object using the
  ## list of functions in `extractors'.  If no extractors are
  ## specified, just return the entire fitted model object.
  rval <- if(is.null(extractors))
    fit
  else 
    lapply(extractors, function(f) f(fit))
  invisible(rval)
  return(fit)
}

setupFormula <- function(cause, pollutant, df.Time, df.time, df.Temp, df.Dew) {
  covariates.f <- paste(cause, "~ dow + agecat")
  
  ## Smooth functions of temperature and dew points (w/running means)
  weather.f <- paste(c(paste("ns(tmpd,", df.Temp, ")"),
                       paste("ns(rmtmpd,", df.Temp, ")"),
                       paste("ns(dptp,", df.Dew, ")"),
                       paste("ns(rmdptp,", df.Dew, ")")), 
                     collapse = "+")
  poll.f <- paste(pollutant, collapse = "+")
  ## Smooth function(s) of time; separate ones for age categories 2, 3
  time.f <- paste(c(paste("ns(time,", df.Time, ")"),
                    paste("I(ns(time,", df.time, ")*Age2Ind)"),
                    paste("I(ns(time,", df.time, ")*Age3Ind)")),
                  collapse = "+")
  form.str <- paste(c(covariates.f, time.f, weather.f, poll.f),
                    collapse = "+")
  as.formula(form.str)
}

## Return a discounted df (based on 12 consecutive days of missings)

numdf <- function(usedata, num = 7){
  n <- length(usedata)
  use <- usedata[1:(n/3)]        
  ll  <- round(length(use)/12)
  
  ## this is to eliminate the warning message the length of use is
  ## not a multiple of 12
  usenew <- use[1:(ll*12)]
  
  mat <- matrix(usenew, ncol = 12,byrow = TRUE)
  m   <- sum(ceiling(apply(mat, 1, sum, na.rm = TRUE)/12)) ##-365.25/12   
  df  <- round(12*m/365.25*num)   
  max(1, df)
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
  varList <- c("cvd", "death", "resp", "tmpd", "rmtmpd", "dptp",
               "rmdptp", "time", "agecat", "dow", "pm10tmean", paste("l",
                                                                     1:7, "pm10tmean", sep = ""))
  data.frame(dataframe[, varList], Age2Ind = Age2Ind, Age3Ind = Age3Ind)
}