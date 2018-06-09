## RcppDE [![Build Status](https://travis-ci.org/eddelbuettel/rcppde.svg)](https://travis-ci.org/eddelbuettel/rcppde) [![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](https://www.gnu.org/licenses/gpl-2.0.html) [![CRAN](http://www.r-pkg.org/badges/version/RcppDE)](https://cran.r-project.org/package=RcppDE) [![Downloads](http://cranlogs.r-pkg.org/badges/RcppDE?color=brightgreen)](http://www.r-pkg.org/pkg/RcppDE)

Rcpp port of Differential Evolution

### About

The package provides global optimization by differential evolution.

It uses an efficient C++ based implementation of the DEoptim
function which performs global optimization by differential evolution.  Its
creation was motivated by trying to see if the old approximation "easier,
shorter, faster: pick any two" could in fact be extended to achieving all
three goals while moving the code from plain old C to modern C++.  The
initial version did in fact do so, but a good part of the gain was due to an
implicit code review which eliminated a few inefficiencies which have since
been eliminated in DEoptim.

## Author 

Dirk Eddelbuettel extending DEoptim by David Ardia, Katharine Mullen, Brian
Peterson and Joshua Ulrich, which itself is based on DE-Engine by Rainer
Storn.
