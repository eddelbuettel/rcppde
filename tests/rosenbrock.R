library(Rcpp)
library(RcppDE)

set.seed(1234)

nothing <- function() {}
sourceCpp(code = '
#include <Rcpp.h>
          
// [[Rcpp::export]]
void nothing() {}
')

Rosenbrock <- function(x) {
  nothing()
  x1 <- x[1]
  x2 <- x[2]
  100 * (x2 - x1 * x1) ^ 2 + (1 - x1) ^ 2
}

lower <- c(-10, -10)
upper <- -lower
cont <- DEoptim.control(trace = 0)

fit <- DEoptim(Rosenbrock, lower, upper, control = cont)

stopifnot(all.equal(fit$optim$bestmem, c(par1 = 1, par2 = 1)))
