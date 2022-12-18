
// Port of DEoptim (2.0.7) by Ardia et al to Rcpp/RcppArmadillo/Armadillo
// Copyright (C) 2010 - 2022  Dirk Eddelbuettel <edd@debian.org>
//
// DEoptim is Copyright (C) 2009 David Ardia and Katharine Mullen

#include <Rcpp/Lightest>

double genrose(Rcpp::NumericVector x) {       // genrose function in C++
    const double a = 1.0;
    const double b = 100.0;
    int n = x.size();
    double sum = 1.0;
    for (int i=1; i<n; i++) {
        sum += b*( ::pow(x[i-1]*x[i-1] - x[i], 2)) + (x[i] - a)*(x[i] - a);
    }
    return(sum);
}

double wild(Rcpp::NumericVector x) {          // wild function in C++
    int n = x.size();
    double sum = 0.0;
    for (int i=0; i<n; i++) {
        double xsq = x[i]*x[i];
        sum += 10 * ::sin(0.3 * x[i]) * ::sin(1.3 * xsq) + 0.00001 * xsq*xsq + 0.2 * x[i] + 80;
    }
    sum /= n;
    return(sum);
}

double rastrigin(Rcpp::NumericVector x) {     // rastrigin function in C++
    int n = x.size();
    double sum = 20.0;
    for (int i=0; i<n; i++) {
        sum += x[i]+2 - 10*::cos(M_2PI*x[i]);
    }
    return(sum);
}

typedef double (*funcPtr)(Rcpp::NumericVector);

// [[Rcpp::export]]
SEXP putFunPtrInXPtr(const std::string& fstr) { 			// needed for tests/
    if (fstr == "genrose")
        return(Rcpp::XPtr<funcPtr>(new funcPtr(&genrose)));
    else if (fstr == "wild")
        return(Rcpp::XPtr<funcPtr>(new funcPtr(&wild)));
    else
        return(Rcpp::XPtr<funcPtr>(new funcPtr(&rastrigin)));
}
