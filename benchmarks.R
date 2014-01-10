# Benchmarks. Contains code to make benchmark comparisons
# across algorithms. Defines evaluation metrics and creates plots.
source("online-algorithms.R")

run.normal.benchmark <- function(niters=10000) {
  e = get.experiment(name="normal", niters=niters)
  d = e$sample.dataset()
  algos =  list()
  algos[[1]] <- sgd.onlineAlgorithm
  algos[[2]] <- asgd.onlineAlgorithm
  algos[[3]] <-  implicit.onlineAlgorithm
  cols = c("red", "black", "magenta")
  
  for(i in 1:length(algos)) {
    out = NA
    if(i==1) {
      out = run.onlineAlgorithm(dataset=d, experiment=e, algorithm=sgd.onlineAlgorithm)
    } else if(i==2) {
      out = run.onlineAlgorithm(dataset=d, experiment=e, algorithm=asgd.onlineAlgorithm)
    } else {
      out = run.onlineAlgorithm(dataset=d, experiment=e, algorithm=implicit.onlineAlgorithm)
    }
    x = log(1:ncol(out$estimates), base=10)
    y = log(apply(out$estimates, 2, e$risk), base=10)
    if(i==1) {
      plot(x, y, type="l", col=cols[i])
      legend(0.5, -1, legend=c("sgd", "asgd", "implicit"), col=cols, lty=rep(1, 3))
    } else {
      lines(x, y, col=cols[i])
    }
  }
}