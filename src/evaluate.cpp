// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// Port of DEoptim (2.0.7) by Ardia et al to Rcpp/RcppArmadillo/Armadillo
// Copyright (C) 2010  Dirk Eddelbuettel <edd@debian.org>
//
// DEoptim is Copyright (C) 2009 David Ardia and Katharine Mullen

#include <RcppArmadillo.h>

//RcppExport double evaluate(long &l_nfeval, const arma::rowvec & param, SEXP parS, SEXP fcall, SEXP env) {
RcppExport double evaluate(long &l_nfeval, const double *param, SEXP parS, SEXP fcall, SEXP env) {
    Rcpp::NumericVector par(parS); 			// access parS as numeric vector to fill it
    //memcpy(par.begin(), param.memptr(), par.size() * sizeof(double));
    //std::copy(param.begin(), param.end(), par.begin()); // STL way of copying
    memcpy(par.begin(), param, par.size() * sizeof(double));
    SEXP fn = ::Rf_lang2(fcall, par); 			// this could be done with Rcpp 
    SEXP sexp_fvec = ::Rf_eval(fn, env);		// but is still a lot slower right now
    double f_result = Rcpp::as<double>(sexp_fvec);
    if (ISNAN(f_result)) 
	::Rf_error("NaN value of objective function! \nPerhaps adjust the bounds.");
    l_nfeval++;  					// increment function evaluation count 
   return(f_result); 
}

