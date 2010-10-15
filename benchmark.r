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

maxIt <- 250                            # not excessive but so that we get some run-time on simple problems

suppressMessages(library(DEoptim)) 	# the original, currently 2.0.7
suppressMessages(library(RcppDE))    	# the contender

basicDE <- function(n, maxIt, fun) DEoptim::DEoptim(fn=fun, lower=rep(-25, n), upper=rep(25, n),
                                                    control=list(NP=10*n, itermax=maxIt, trace=FALSE))
cppDE <- function(n, maxIt, fun) RcppDE::DEoptim(fn=fun, lower=rep(-25, n), upper=rep(25, n),
                                                 control=list(NP=10*n, itermax=maxIt, trace=FALSE))

runPair <- function(n, maxIt, fun) {

    set.seed(42)
    valBasic <- basicDE(n, maxIt, fun)
    set.seed(42)
    valCpp <- cppDE(n, maxIt, fun)
    stopifnot( all.equal(valBasic, valCpp) )

    gc()
    set.seed(42)
    bt <- mean(replicate(20, system.time(invisible(basicDE(n, maxIt, fun)))[3]), trim=0.05)

    gc()
    set.seed(42)
    ct <- mean(replicate(20, system.time(invisible(cppDE(n, maxIt, fun)))[3]), trim=0.05)

    return(data.frame(DEoptim=bt, RcppDE=ct))
}

res <- rbind(runPair(2, maxIt, function(...) Rastrigin(...)),
             runPair(5, maxIt, function(...) Rastrigin(...)),
             runPair(20, maxIt, function(...) Rastrigin(...)),
             runPair(2, maxIt, function(...) Wild(...)),
             runPair(5, maxIt, function(...) Wild(...)),
             runPair(20, maxIt, function(...) Wild(...)),
             runPair(2, maxIt, function(...) Genrose(...)),
             runPair(5, maxIt, function(...) Genrose(...)),
             runPair(20, maxIt, function(...) Genrose(...)),
             runPair(50, maxIt, function(...) Genrose(...))
#             runPair(100, maxIt, function(...) Genrose(...))
             )
res <- rbind(res, colMeans(res))
rownames(res) <- c("Rastrigin2", "Rastrigin5", "Rastrigin20",
                   "Wild2", "Wild5", "Wild20",
                   "Genrose2", "Genrose5", "Genrose20", "Genrose50", #"Genrose100",
                   "MEANS")
res$ratioRcppToBasic <- res[,2]/res[,1]
res$pctGainOfRcpp <- (1-res[,2]/res[,1])*100

print(res)

