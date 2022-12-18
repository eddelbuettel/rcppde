
// Port of DEoptim (2.0.7) by Ardia et al to Rcpp/RcppArmadillo/Armadillo
// Copyright (C) 2010 - 2022  Dirk Eddelbuettel <edd@debian.org>
//
// DEoptim is Copyright (C) 2009 David Ardia and Katharine Mullen

#ifndef Rcpp_DE_evaluate_h_
#define Rcpp_DE_evaluate_h_

#include <Rcpp.h>

namespace Rcpp {
    namespace DE {

        class EvalBase {
        public:
            EvalBase() : neval(0) {};
            virtual double eval(NumericVector par) = 0;
            unsigned long getNbEvals() { return neval; }
        protected:
            unsigned long int neval;
        };

        class EvalStandard : public EvalBase {
        public:
            EvalStandard(SEXP fcall_, SEXP env_) : fcall(fcall_), env(env_) {}
            double eval(NumericVector par) {
                neval++;
                return defaultfun(par);
            }
        private:
            SEXP fcall, env;
            double defaultfun(NumericVector par) {        // essentialy same as the old evaluate
                Rcpp::Shield<SEXP> fn(::Rf_lang3(fcall, par, R_DotsSymbol));
                Rcpp::Shield<SEXP> sexp_fvec(::Rf_eval(fn, env));
                double f_result = REAL(sexp_fvec)[0];
                if (ISNAN(f_result))
                    ::Rf_error("NaN value of objective function! \nPerhaps adjust the bounds.");
                return(f_result);
            }
        };

        typedef double (*funcPtr)(NumericVector);

        class EvalCompiled : public EvalBase {
        public:
            EvalCompiled(Rcpp::XPtr<funcPtr> xptr, SEXP env_) {
                funptr = *(xptr);
                env = env_;
            };
            EvalCompiled(SEXP xps, SEXP env_) {
                Rcpp::XPtr<funcPtr> xptr(xps);
                funptr = *(xptr);
                env = env_;
            };
            double eval(NumericVector par) {
                neval++;
                return funptr(par); //, env);
            }
        private:
            funcPtr funptr;
            SEXP env;
        };

    }

}


#endif
