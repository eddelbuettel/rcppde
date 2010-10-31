// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// Port of DEoptim (2.0.7) by Ardia et al to Rcpp/RcppArmadillo/Armadillo
// Copyright (C) 2010  Dirk Eddelbuettel <edd@debian.org>
//
// DEoptim is Copyright (C) 2009 David Ardia and Katharine Mullen
// and based on DE-Engine v4.0, Rainer Storn, 2004  
// (http://www.icsi.berkeley.edu/~storn/DeWin.zip)

#include <RcppArmadillo.h>

// Function       : void permute(int ia_urn2[], int i_urn2_depth)
// Author         : Rainer Storn (w/bug fixes contributed by DEoptim users)
// Description    : Generates i_urn2_depth random indices ex [0, i_NP-1]
//                  which are all distinct. This is done by using a
//                  permutation algorithm called the "urn algorithm"
//                  which goes back to C.L.Robinson.
// Functions      : -
// Globals        : -
// Parameters     : ia_urn2       (O)    array containing the random indices
//                  i_urn2_depth  (I)    number of random indices (avoided index included)
//                  i_NP          (I)    range of indices is [0, i_NP-1]
//                  i_avoid       (I)    is the index to avoid and is located in ia_urn2[0].
//                  ia_urn1       (I)    additional temp vector
// Preconditions  : # Make sure that ia_urn2[] has a length of i_urn2_depth.
//                  # i_urn2_depth must be smaller than i_NP.
// Postconditions : # the index to be avoided is in ia_urn2[0], so fetch the
//                   indices from ia_urn2[i], i = 1, 2, 3, ..., i_urn2_depth.
// Return Value   : -
void permute(int ia_urn2[], int i_urn2_depth, int i_NP, int i_avoid, int ia_urn1[]) {
    GetRNGstate();
    int k = i_NP;
    int i_urn1 = 0;
    int i_urn2 = 0;
    for (int i = 0; i < i_NP; i++)
	ia_urn1[i] = i; 		   /* initialize urn1 */

    i_urn1 = i_avoid;                      /* get rid of the index to be avoided and place it in position 0. */
    while (k > i_NP - i_urn2_depth) {      /* i_urn2_depth is the amount of indices wanted (must be <= NP) */
	ia_urn2[i_urn2] = ia_urn1[i_urn1]; /* move it into urn2 */
	ia_urn1[i_urn1] = ia_urn1[k-1];    /* move highest index to fill gap */
	k = k - 1;                         /* reduce number of accessible indices */
	i_urn2 = i_urn2 + 1;               /* next position in urn2 */
	i_urn1 = static_cast<int>(::unif_rand() * k);   /* choose a random index */
    }
    PutRNGstate();
}
