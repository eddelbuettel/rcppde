// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// Port of DEoptim (2.0.7) by Ardia et al to Rcpp/RcppArmadillo/Armadillo
// Copyright (C) 2010  Dirk Eddelbuettel <edd@debian.org>
//
// DEoptim is Copyright (C) 2009 David Ardia and Katharine Mullen

#ifndef Rcpp_DE_evaluate_h_
#define Rcpp_DE_evaluate_h_

#include <Rcpp.h>

namespace Rcpp {
    namespace DE {
	double genrose(SEXP xs) {	// genrose function in C++
	    Rcpp::NumericVector x(xs);
	    int n = x.size();
	    double sum = 1.0;
	    for (int i=1; i<n; i++) {
		sum += 100*( ::pow(x[i-1]*x[i-1] - x[i], 2)) + (x[i] - 1)*(x[i] - 1);
	    }
	    return(sum);
	}

	double wild(SEXP xs) {		// wild function in C++
	    Rcpp::NumericVector x(xs);
	    int n = x.size();
	    double sum = 0.0;
	    for (int i=0; i<n; i++) {
		sum += 10 * ::sin(0.3 * x[i]) * ::sin(1.3 * x[i]*x[i]) + 0.00001 * x[i]*x[i]*x[i]*x[i] + 0.2 * x[i] + 80;
	    }
	    sum /= n;
	    return(sum);
	}

	double rastrigin(SEXP xs) {	// rastrigin function in C++
	    Rcpp::NumericVector x(xs);
	    int n = x.size();
	    double sum = 20.0;
	    for (int i=0; i<n; i++) {
		sum += x[i]+2 - 10*::cos(2*M_PI*x[i]);
	    }
	    return(sum);
	}

#if 0
class Evaluator {		// class to evaluate a given function at a parameter
public:
    typedef double (Evaluator::*FunctionPointer)(SEXP);
    Evaluator( FunctionPointer funptr_ ) : funptr(funptr_) {};
    Evaluator(SEXP fcall_, SEXP env_) {
	if (TYPEOF(env_) == ENVSXP) { 		// standard mode: env_ is an env, fcall_ is a function 
	    //REprintf("In env case, default\n");
	    fcall = fcall_;
	    env = env_;
	    funptr = &Evaluator::defaultfun;
	} else {
	    REprintf("NOT in env case, trying something new -- does not work yet\n");
	    fcall = fcall_;
	    env = env_;
	    funptr = &Evaluator::defaultfun;
	}
    };
    double eval(SEXP par) {
	return ((*this).*(funptr))(par); 	// isn't the syntax to eval a function pointer easy :) 
    }
    inline FunctionPointer get() { return funptr ; }
private:
    SEXP fcall, env;
    double defaultfun(SEXP par) {
	SEXP fn = ::Rf_lang2(fcall, par); 			   // this could be done with Rcpp 
	SEXP sexp_fvec = ::Rf_eval(fn, env);		           // but is still a lot slower right now
	double f_result = REAL(sexp_fvec)[0];
	if (ISNAN(f_result)) 
	    ::Rf_error("NaN value of objective function! \nPerhaps adjust the bounds.");
	return(f_result); 
    }
    FunctionPointer funptr;
};
#endif

	class Fun {			// class to wrap an external pointer around eval. function
	public:
	    typedef double (*FunctionPointer)(SEXP);
	    Fun( FunctionPointer ptr_ ) : ptr(ptr_) {};
	    inline FunctionPointer get() { return ptr ; }
	private:
	    FunctionPointer ptr ;
	};

	class EvalBase {
	public:
	    virtual double eval(SEXP par) = 0;
	};

	class EvalStandard : public EvalBase {
	public:
	    //typedef double (EvalStandard::*FunctionPointer)(SEXP);
	    EvalStandard(SEXP fcall_, SEXP env_) : fcall(fcall_), env(env_) { 
		//funptr = &EvalStandard::defaultfun;
	    } 
	    double eval(SEXP par) {
		//return ((*this).*(funptr))(par); 	// isn't the syntax to eval a function pointer easy :) 
		return defaultfun(par);
	    }
	private:
	    SEXP fcall, env;
	    double defaultfun(SEXP par) {
		SEXP fn = ::Rf_lang2(fcall, par); 			   // this could be done with Rcpp 
		SEXP sexp_fvec = ::Rf_eval(fn, env);		           // but is still a lot slower right now
		double f_result = REAL(sexp_fvec)[0];
		if (ISNAN(f_result)) 
		    ::Rf_error("NaN value of objective function! \nPerhaps adjust the bounds.");
		return(f_result); 
	    }
	    //FunctionPointer funptr;
	};

	class EvalCompiled : public EvalBase {
	public:
	    EvalCompiled( Rcpp::XPtr<Fun> xptr ) {
		funptr = xptr->get();
	    };
	    EvalCompiled( SEXP xps ) {
		Rcpp::XPtr<Fun> xptr(xps);
		funptr = xptr->get();
	    };
	    double eval(SEXP par) {
		return funptr(par);
	    }
	private:
	    Fun::FunctionPointer funptr;
	};

	RcppExport SEXP putFunPtrInXPtr(SEXP funname) {
	    std::string fstr = Rcpp::as<std::string>(funname);
	    if (fstr == "genrose")
		return(Rcpp::XPtr<Fun>(new Fun(&genrose)));
	    else if (fstr == "wild")
		return(Rcpp::XPtr<Fun>(new Fun(&wild)));
	    else
		return(Rcpp::XPtr<Fun>(new Fun(&rastrigin)));
	}

    }

}

#endif
