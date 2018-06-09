library(Rcpp)
library(RcppDE)

set.seed(1234)

Rosenbrock <- function(x) {
  
  # briefly, calls to functions generated with Rcpp Attributes
  # could break optimizers due to incorrect sync of RNG state
  RcppDE:::putFunPtrInXPtr("dummy")
  
  x1 <- x[1]
  x2 <- x[2]
  100 * (x2 - x1 * x1) ^ 2 + (1 - x1) ^ 2
}

lower <- c(-10, -10)
upper <- -lower
cont <- DEoptim.control(trace = 0)

fit <- DEoptim(Rosenbrock, lower, upper, control = cont)

if (packageVersion("Rcpp") >= "0.12.17.1")
  stopifnot(all.equal(fit$optim$bestmem, c(par1 = 1, par2 = 1)))
