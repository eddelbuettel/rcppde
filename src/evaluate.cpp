// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// Port of DEoptim (2.0.7) by Ardia et al to Rcpp/RcppArmadillo/Armadillo
// Copyright (C) 2010  Dirk Eddelbuettel <edd@debian.org>
//
// DEoptim is Copyright (C) 2009 David Ardia and Katharine Mullen

#include <RcppArmadillo.h>

// Slighly accelerated version of evaluate to evaluate function fcall with parameters param in environment env
// Uses externally allocated par() vector into which param are copied
//
double evaluate(long &l_nfeval, const double *param, SEXP par, SEXP fcall, SEXP env) {
    // -- safer approach: cast to NumericVector, fill it and eval
    //Rcpp::NumericVector parvec(par); 			// access parS as numeric vector to fill it
    //memcpy(parvec.begin(), param, parvec.size() * sizeof(double));
    //SEXP fn = ::Rf_lang2(fcall, parvec); 			// this could be done with Rcpp 

    // -- faster: direct access _assuming_ numeric vector
    memcpy(REAL(par), param, Rf_nrows(par) * sizeof(double));

    SEXP fn = ::Rf_lang2(fcall, par); 			// this could be done with Rcpp 
    SEXP sexp_fvec = ::Rf_eval(fn, env);		// but is still a lot slower right now

    double f_result = Rcpp::as<double>(sexp_fvec);
    if (ISNAN(f_result)) 
	::Rf_error("NaN value of objective function! \nPerhaps adjust the bounds.");
    l_nfeval++;  					// increment function evaluation count 
   return(f_result); 
}

