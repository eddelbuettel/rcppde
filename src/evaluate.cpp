// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// Port of DEoptim (2.0.7) to Rcpp/RcppArmadillo/Armadillo
// Copyright (C) 2010  Dirk Eddelbuettel <edd@debian.org>

#include <Rcpp.h>

/*------objective function---------------------------------------*/

RcppExport double evaluate(long &l_nfeval, double *param, SEXP parS, SEXP fcall, SEXP env) {
    Rcpp::NumericVector par(parS);
    memcpy(par.begin(), param, par.size() * sizeof(double));

    l_nfeval++;  		// increment function evaluation count 

    SEXP fn = ::Rf_lang2(fcall, par); 			// this can be done with Rcpp 
    SEXP sexp_fvec = ::Rf_eval(fn, env);		// but is still a lot slower right now
    double f_result = Rcpp::as<double>(sexp_fvec);
   
    if (ISNAN(f_result)) 
	::Rf_error("NaN value of objective function! \nPerhaps adjust the bounds.");
   
   return(f_result); 
}

