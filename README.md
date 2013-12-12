# RcppDE: Rcpp port of Differential Evolution

[![Build Status](https://travis-ci.org/eddelbuettel/rcppde.png)](https://travis-ci.org/eddelbuettel/rcppde)

## Global optimization by differential evolution in C++

This package provides an efficient C++ based implementation of the DEoptim
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
