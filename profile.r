#!/usr/bin/r -t

suppressMessages(library(RcppDE))

Genrose <- function(x) { 	## One generalization of the Rosenbrock banana valley function (n parameters)
    n <- length(x)
    1.0 + sum (100 * (x[-n]^2 - x[-1])^2 + (x[-1] - 1)^2)
}


maxIt <- 2500                            # not excessive but so that we get some run-time on simple problems
n <- 100

RcppDE::DEoptim(fn=Genrose, lower=rep(-25, n), upper=rep(25, n), control=list(NP=10*n, itermax=maxIt, trace=FALSE))
