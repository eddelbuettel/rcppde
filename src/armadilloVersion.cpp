
#include <RcppArmadillo/Lighter>

// [[Rcpp::export]]
Rcpp::List armadilloVersion() {
    // create a vector of major, minor, patch
    auto ver = Rcpp::IntegerVector::create(ARMA_VERSION_MAJOR, ARMA_VERSION_MINOR, ARMA_VERSION_PATCH);
    // and place it in a list (as e.g. packageVersion() in R returns)
    auto lst = Rcpp::List::create(ver);
    // and class it as 'package_version' accessing print() etc methods
    lst.attr("class") = Rcpp::CharacterVector::create("package_version", "numeric_version");
    return lst;
}
