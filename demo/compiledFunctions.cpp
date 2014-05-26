
#include <Rcpp.h>

// [[Rcpp::interfaces(r, cpp)]]

// Note that the interface is double fun(SEXP) as it comes from a time
// before Rcpp Attributes. We could rewrite it for NumericVector or
// Armadillo vector arguments but would then have to do so consistenly
// in the RcppDE sources as well where the XPtr gets unwrapped.

double genrose(SEXP xs) {
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

// [[Rcpp::export]]
SEXP create_xptr(std::string fstr) {
  typedef double (*funcPtr)(SEXP);
  if (fstr == "genrose")
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&genrose)));
  else if (fstr == "wild")
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&wild)));
  else
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&rastrigin)));
}  
