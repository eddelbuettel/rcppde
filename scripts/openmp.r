#!/usr/bin/r -t
#
# with OpenMP we do not get the same uniform random number draws as we act in parallel, so just compare results (and timings)

suppressMessages(library(DEoptim)) 	# the original, currently 2.0.7
suppressMessages(library(RcppDE))    	# the contender
suppressMessages(library(inline))    	# helps with

Genrose <- function(x) { 		## One generalization of the Rosenbrock banana valley function (n parameters)
    n <- length(x)
    1.0 + sum (100 * (x[-n]^2 - x[-1])^2 + (x[-1] - 1)^2)
}


#maxIt <- 500                           # not excessive but so that we get some run-time on simple problems
#n <- 40

maxIt <- 500
n <- 20

inc <- 'double genrose(SEXP xs) {
                Rcpp::NumericVector x(xs);
                int n = x.size();
                double sum = 1.0;
                for (int i=1; i<n; i++) {
                   sum += 100*( pow(x[i-1]*x[i-1] - x[i], 2)) + (x[i] - 1)*(x[i] - 1);
                }
                return(sum);
             }

             double wild(SEXP xs) {
                Rcpp::NumericVector x(xs);
                int n = x.size();
                double sum = 0.0;
                for (int i=0; i<n; i++) {
                   sum += 10 * sin(0.3 * x[i]) * sin(1.3 * x[i]*x[i]) + 0.00001 * x[i]*x[i]*x[i]*x[i] + 0.2 * x[i] + 80;
                }
                sum /= n;
                return(sum);
             }

             double rastrigin(SEXP xs) {
                Rcpp::NumericVector x(xs);
                int n = x.size();
                double sum = 20.0;
                for (int i=0; i<n; i++) {
                   sum += x[i]+2 - 10*cos(2*M_PI*x[i]);
                }
                return(sum);
             }

             '

## now via a class returning external pointer
src.xptr <- 'std::string fstr = Rcpp::as<std::string>(funname);
	         typedef double (*funcPtr)(SEXP);
                 if (fstr == "genrose")
                     return(XPtr<funcPtr>(new funcPtr(&genrose)));
                 else if (fstr == "wild")
                     return(XPtr<funcPtr>(new funcPtr(&wild)));
                 else
                     return(XPtr<funcPtr>(new funcPtr(&rastrigin)));
                 '

create_xptr <- cxxfunction(signature(funname="character"), body=src.xptr, inc=inc, plugin="Rcpp")
xptr <- create_xptr("genrose")

set.seed(42)
print(system.time( {
    res <- RcppDE::DEoptim(fn=xptr, lower=rep(-25, n), upper=rep(25, n), control=list(NP=10*n, itermax=maxIt, trace=FALSE))
    print(res[[1]][["bestval"]])
}))

set.seed(42)
print(system.time( {
    res <- DEoptim::DEoptim(fn=Genrose, lower=rep(-25, n), upper=rep(25, n), control=list(NP=10*n, itermax=maxIt, trace=FALSE))
    print(res[[1]][["bestval"]])
}))

