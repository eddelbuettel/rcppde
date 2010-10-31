// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// Port of DEoptim (2.0.7) by Ardia et al to Rcpp/RcppArmadillo/Armadillo
// Copyright (C) 2010  Dirk Eddelbuettel <edd@debian.org>
//
// DEoptim is Copyright (C) 2009 David Ardia and Katharine Mullen

#ifdef USE_OPENMP
#include <RcppArmadillo.h>

// Slighly accelerated version of evaluate to evaluate function fcall with parameters param in environment env
// Uses externally allocated par() vector into which param are copied
double evaluate(long &l_nfeval, const double *param, SEXP par, SEXP fcall, SEXP env) {
    memcpy(REAL(par), param, Rf_nrows(par) * sizeof(double));      // -- faster: direct access _assuming_ numeric vector
    Rcpp::NumericVector x(par);
    int n = x.size();
    double sum = 1.0;
    for (int i=1; i<n; i++) {
	sum += 100*( pow(x[i-1]*x[i-1] - x[i], 2)) + (x[i] - 1)*(x[i] - 1);
    }
    l_nfeval++;  						   // increment function evaluation count 
    return(sum);
}
#endif

