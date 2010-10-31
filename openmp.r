#!/usr/bin/r -t
#
# with OpenMP we do not get the same uniform random number draws as we act in parallel, so just compare results (and timings)

suppressMessages(library(DEoptim)) 	# the original, currently 2.0.7
suppressMessages(library(RcppDE))    	# the contender

Genrose <- function(x) { 	## One generalization of the Rosenbrock banana valley function (n parameters)
    n <- length(x)
    1.0 + sum (100 * (x[-n]^2 - x[-1])^2 + (x[-1] - 1)^2)
}


maxIt <- 500                           # not excessive but so that we get some run-time on simple problems
n <- 20

set.seed(42)
print(system.time( {
    res <- RcppDE::DEoptim(fn=Genrose, lower=rep(-25, n), upper=rep(25, n), control=list(NP=10*n, itermax=maxIt, trace=FALSE))
    print(res[[1]])
}))

set.seed(42)
print(system.time( {
    res <- DEoptim::DEoptim(fn=Genrose, lower=rep(-25, n), upper=rep(25, n), control=list(NP=10*n, itermax=maxIt, trace=FALSE))
    print(res[[1]])
}))

