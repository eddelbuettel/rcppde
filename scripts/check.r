#!/usr/bin/r -t

Wild <- function(x) { 		## 'Wild' function, global minimum at about -15.81515
    sum(10 * sin(0.3 * x) * sin(1.3 * x^2) + 0.00001 * x^4 + 0.2 * x + 80)/length(x)
}

Rastrigin <- function(x) {
    sum(x+2 - 10 * cos(2*pi*x)) + 20
}

Genrose <- function(x) { 	## One generalization of the Rosenbrock banana valley function (n parameters)
    n <- length(x)
    1.0 + sum (100 * (x[-n]^2 - x[-1])^2 + (x[-1] - 1)^2)
}

#maxIt <- 100                            # not excessive but so that we get some run-time on simple problems
n <- 20
maxIt <- 50
useBS <- TRUE
storeFrom <- maxIt+1
strat <- 6                              # TODO fix segfault when strat==6

suppressMessages(library(DEoptim)) 	# the original, currently 2.0.7
suppressMessages(library(RcppDE))    	# the contender

ctrl <- DEoptim::DEoptim.control(NP=10*n,
                                 itermax=maxIt,
                                 trace=FALSE,
                                 bs=useBS,
                                 storepopfrom=storeFrom,
                                 strategy=strat)

basicDE <- function(n, maxIt, fun) DEoptim::DEoptim(fn=fun, lower=rep(-25, n), upper=rep(25, n), control=ctrl)
cppDE <- function(n, maxIt, fun) RcppDE::DEoptim(fn=fun, lower=rep(-25, n), upper=rep(25, n), control=ctrl)

set.seed(42)
valBasic <- basicDE(n, maxIt, function(...) Rastrigin(...))
#print(str(valBasic[[2]]))
set.seed(42)
valCpp <- cppDE(n, maxIt, function(...) Rastrigin(...))
#print(str(valCpp[[2]]))
stopifnot( all.equal(valBasic[[1]], valCpp[[1]]) )

cat("# Done ", format(Sys.time()), "\n")







