// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// Port of DEoptim (2.0.7) by Ardia et al to Rcpp/RcppArmadillo/Armadillo
// Copyright (C) 2010  Dirk Eddelbuettel <edd@debian.org>
//
// DEoptim is Copyright (C) 2009 David Ardia and Katharine Mullen
// and based on DE-Engine v4.0, Rainer Storn, 2004  
// (http://www.icsi.berkeley.edu/~storn/DeWin.zip)

#ifdef USE_OPENMP
#include <RcppArmadillo.h>	// declarations for both Rcpp and RcppArmadillo offering Armadillo classes
#include <omp.h>		// OpenMP for compiler-generated multithreading

void permute(int ia_urn2[], int i_urn2_depth, int i_NP, int i_avoid, int ia_urntmp[]);
double evaluate(long &l_nfeval, const double *param, SEXP parS, SEXP fcall, SEXP env);

void devol(double VTR, double f_weight, double f_cross, int i_bs_flag,
           arma::colvec & fa_minbound, arma::colvec & fa_maxbound, SEXP fcall, SEXP rho, int i_trace,
           int i_strategy, int i_D, int i_NP, int i_itermax, arma::mat & initialpopm, 
	   int i_storepopfrom, int i_storepopfreq, int i_specinitialpop, int i_check_winner, int i_av_winner,
           arma::mat &ta_popP, arma::mat &ta_oldP, arma::mat &ta_newP, arma::colvec & t_bestP, 
           arma::colvec & ta_popC, arma::colvec & ta_oldC, arma::colvec & ta_newC, double & t_bestC,
           arma::colvec & t_bestitP, arma::colvec & t_tmpP, 
           arma::mat &d_pop, Rcpp::List &d_storepop, arma::mat & d_bestmemit, arma::colvec & d_bestvalit,
           int & i_iterations, double i_pPct, long & l_nfeval) {

    //ProfilerStart("/tmp/RcppDE.prof");
    const int urn_depth = 5;   			// 4 + one index to avoid 
    Rcpp::NumericVector par(i_D);		// initialize parameter vector to pass to evaluate function 
    arma::icolvec::fixed<urn_depth> ia_urn2; 	// fixed-size vector for urn draws
    arma::icolvec ia_urntmp(i_NP); 		// so that we don't need to re-allocated each time in permute
    arma::mat initialpop(i_D, i_NP); 
    int i_nstorepop = ceil((i_itermax - i_storepopfrom) / i_storepopfreq);
    int p_NP = round(i_pPct * i_NP);  		// choose at least two best solutions 
    p_NP = p_NP < 2 ? 2 : p_NP;
    arma::icolvec sortIndex(i_NP); 		// sorted values of ta_oldC 
    if (i_strategy == 6) {
	for (int i = 0; i < i_NP; i++) 
	    sortIndex[i] = i; 
    }
    GetRNGstate();

    initialpop.zeros();		 		// initialize initial popuplation 
    d_bestmemit.zeros();    			// initialize best members
    d_bestvalit.zeros();			// initialize best values 
    d_pop.zeros();				// initialize best population
    i_nstorepop = (i_nstorepop < 0) ? 0 : i_nstorepop;
      
    if (i_specinitialpop > 0) {    		// if initial population provided, initialize with values 
	initialpop = trans(initialpopm);	// transpose as we prefer columns for population members here
    }

    for (int i = 0; i < i_NP; i++) {		// ------Initialization-----------------------------
	if (i_specinitialpop <= 0) { 		// random initial member 
	    for (int j = 0; j < i_D; j++) {
		ta_popP.at(j,i) = fa_minbound[j] + ::unif_rand() * (fa_maxbound[j] - fa_minbound[j]);
	    }
	} else { 				// or user-specified initial member 
	    ta_popP.col(i) = initialpop.col(i);
	} 
	ta_popC[i] = evaluate(l_nfeval, ta_popP.colptr(i), par, fcall, rho);
	if (i == 0 || ta_popC[i] <= t_bestC) {
	    t_bestC = ta_popC[i];
	    t_bestP = ta_popP.unsafe_col(i);
	}
    }

    ta_oldP = ta_popP.cols(0, i_NP-1);		// ---assign pointers to current ("old") population---
    ta_oldC = ta_popC.rows(0, i_NP-1);
  
    int i_iter = 0;				// ------Iteration loop--------------------------------------------
    int popcnt = 0;
    int i_xav = 1;
  
    for (i_iter=0; (i_iter < i_itermax) && (t_bestC > VTR); i_iter++) { // main loop ====================================
	if (i_iter % i_storepopfreq == 0 && i_iter >= i_storepopfrom) {  	// store intermediate populations
	    d_storepop[popcnt++] = Rcpp::wrap( trans(ta_oldP) );
	} // end store pop 

	d_bestmemit.col(i_iter) = t_bestP;	// store the best member
	d_bestvalit[i_iter] = t_bestC;		// store the best value 
	t_bestitP = t_bestP;
	//i_iter++;				// increase iteration counter
     
	double f_dither = f_weight + ::unif_rand() * (1.0 - f_weight);	// ----computer dithering factor -----------------
      
	if (i_strategy == 6) {			// ---DE/current-to-p-best/1 -----------------------------------------------------
	    arma::colvec temp_oldC = ta_oldC;					// create a copy of ta_oldC to avoid changing it 
	    rsort_with_index( temp_oldC.memptr(), sortIndex.begin(), i_NP );  	// sort temp_oldC to use sortIndex later 
	}

#pragma omp parallel for shared(ta_oldP,ta_newP,ta_newC,t_tmpP) private(i) schedule(dynamic)
	for (int i = 0; i < i_NP; i++) {	// ----start of loop through ensemble------------------------

	    t_tmpP = ta_oldP.col(i);		// t_tmpP is the vector to mutate and eventually select

	    permute(ia_urn2.memptr(), urn_depth, i_NP, i, ia_urntmp.memptr()); // Pick 4 random and distinct 
	    int k = 0;				// loop counter used in all strategies below 

	    // ===Choice of strategy=======================================================
	    switch (i_strategy) {

	    case 1: {				// ---classical strategy DE/rand/1/bin-----------------------------------------
		int j = static_cast<int>(::unif_rand() * i_D); 	// random parameter 
		do {				// add fluctuation to random target 
		    t_tmpP[j] = ta_oldP.at(j,ia_urn2[1]) + f_weight * (ta_oldP.at(j,ia_urn2[2]) - ta_oldP.at(j,ia_urn2[3]));
		    j = (j + 1) % i_D;
		} while ((::unif_rand() < f_cross) && (++k < i_D));
		break;
	    }
	    case 2: {				// ---DE/local-to-best/1/bin---------------------------------------------------
		int j = static_cast<int>(::unif_rand() * i_D); 	// random parameter 
		do {				// add fluctuation to random target 
		    t_tmpP[j] = t_tmpP[j] + f_weight * (t_bestitP[j] - t_tmpP[j]) + f_weight * (ta_oldP.at(j,ia_urn2[2]) - ta_oldP.at(j,ia_urn2[3]));
		    j = (j + 1) % i_D;
		} while ((::unif_rand() < f_cross) && (++k < i_D));
		break;
	    }
	    case 3: {				// ---DE/best/1/bin with jitter------------------------------------------------
		int j = static_cast<int>(::unif_rand() * i_D); 	// random parameter 
		do {				// add fluctuation to random target 
		    double f_jitter = 0.0001 * ::unif_rand() + f_weight; 
		    t_tmpP[j] = t_bestitP[j] + f_jitter * (ta_oldP.at(j,ia_urn2[1]) - ta_oldP.at(j,ia_urn2[2]));
		    j = (j + 1) % i_D;
		} while ((::unif_rand() < f_cross) && (++k < i_D));
		break;
	    }
	    case 4: {				// ---DE/rand/1/bin with per-vector-dither-------------------------------------
		int j = static_cast<int>(::unif_rand() * i_D); 	// random parameter 
		do {				// add fluctuation to random target *
		    t_tmpP[j] = ta_oldP.at(j,ia_urn2[1]) + (f_weight + ::unif_rand()*(1.0 - f_weight))* (ta_oldP.at(j,ia_urn2[2]) - ta_oldP.at(j,ia_urn2[3]));
		    j = (j + 1) % i_D;
		} while ((::unif_rand() < f_cross) && (++k < i_D));
		break;
	    }
	    case 5: {				// ---DE/rand/1/bin with per-generation-dither---------------------------------
		int j = static_cast<int>(::unif_rand() * i_D); 	// random parameter 
		do {				// add fluctuation to random target 
		    t_tmpP[j] = ta_oldP.at(j,ia_urn2[1]) + f_dither * (ta_oldP.at(j,ia_urn2[2]) - ta_oldP.at(j,ia_urn2[3]));
		    j = (j + 1) % i_D;
		} while ((::unif_rand() < f_cross) && (++k < i_D));
		break;
	    }
	    case 6: {				// ---DE/current-to-p-best/1 (JADE)--------------------------------------------
		int i_pbest = sortIndex[static_cast<int>(::unif_rand() * p_NP)]; // select from [0, 1, 2, ..., (pNP-1)] 
		int j = static_cast<int>(::unif_rand() * i_D); 	// random parameter 
		do {				// add fluctuation to random target 
		    t_tmpP[j] = ta_oldP.at(j,i) + f_weight * (ta_oldP.at(j,i_pbest) - ta_oldP.at(j,i)) + f_weight * (ta_oldP.at(j,ia_urn2[1]) - ta_oldP.at(j,ia_urn2[2]));
		    j = (j + 1) % i_D;
		} while ((::unif_rand() < f_cross) && (++k < i_D));
		break;
	    }
	    default: {				// ---variation to DE/rand/1/bin: either-or-algorithm--------------------------
		int j = static_cast<int>(::unif_rand() * i_D); 	// random parameter 
		if (::unif_rand() < 0.5) { 	// differential mutation, Pmu = 0.5 
		    do {			// add fluctuation to random target */
			t_tmpP[j] = ta_oldP.at(j,ia_urn2[1]) + f_weight * (ta_oldP.at(j,ia_urn2[2]) - ta_oldP.at(j,ia_urn2[3]));
			j = (j + 1) % i_D;
		    } while ((::unif_rand() < f_cross) && (++k < i_D));

		} else { 			// recombination with K = 0.5*(F+1) -. F-K-Rule 
		    do {			// add fluctuation to random target */
			t_tmpP[j] = ta_oldP.at(j,ia_urn2[1]) + 0.5 * (f_weight + 1.0) * (ta_oldP.at(j,ia_urn2[2]) + ta_oldP.at(j,ia_urn2[3]) - 2 * ta_oldP.at(j,ia_urn2[1]));
			j = (j + 1) % i_D;
		    } while ((::unif_rand() < f_cross) && (++k < i_D));
		}
		break;
	    }
	    } // end switch (i_strategy) ...
	
	    for (int j = 0; j < i_D; j++) {		// ----boundary constraints, bounce-back method was not enforcing bounds correctly
		if (t_tmpP[j] < fa_minbound[j]) {
		    t_tmpP[j] = fa_minbound[j] + ::unif_rand() * (fa_maxbound[j] - fa_minbound[j]);
		}
		if (t_tmpP[j] > fa_maxbound[j]) {
		    t_tmpP[j] = fa_maxbound[j] - ::unif_rand() * (fa_maxbound[j] - fa_minbound[j]);
		}
	    }

	    // ------Trial mutation now in t_tmpP-----------------
	    double t_tmpC = evaluate(l_nfeval, t_tmpP.memptr(), par, fcall, rho);	// Evaluate mutant in t_tmpP[]
	    if (t_tmpC <= ta_oldC[i] || i_bs_flag) {	    		// i_bs_flag means that we will choose best NP vectors from old and new population later
		ta_newP.col(i) = t_tmpP;				// replace target with mutant 
		ta_newC[i] = t_tmpC;
		if (t_tmpC <= t_bestC) {
		    t_bestP = t_tmpP;
		    t_bestC = t_tmpC;
		}
	    } else {
		ta_newP.col(i) = ta_oldP.col(i);
		ta_newC[i] = ta_oldC[i];
	    }
	} // End mutation loop through pop., ie the "for (i = 0; i < i_NP; i++)"

	if (i_bs_flag) {	// examine old and new pop. and take the best NP members into next generation 
	    
	    ta_popP.cols(0, i_NP-1) = ta_oldP;
	    ta_popC.rows(0, i_NP-1) = ta_oldC;

	    ta_popP.cols(i_NP, 2*i_NP-1) = ta_newP;
	    ta_popC.rows(i_NP, 2*i_NP-1) = ta_newC;

	    int i_len = 2 * i_NP;
	    int step = i_len, done;	// array length 
	    while (step > 1) {
		step /= 2;   		// halve the step size 
		do {
		    done = 1;
		    int bound  = i_len - step;
		    for (int j = 0; j < bound; j++) {
			int i = j + step + 1;
			if (ta_popC[j] > ta_popC[i-1]) {
			    ta_popP.swap_cols(j, i-1);
			    ta_popC.swap_rows(j, i-1);
			    done = 0;
			}  // if 
		    }  // for 
		} while (!done); // while
	    } // while (step > 1) 
	    ta_newP = ta_popP.cols(0, i_NP-1);	// now the best NP are in first NP places in gta_pop, use them
	    ta_newC = ta_popC.rows(0, i_NP-1);
	} // i_bs_flag

	ta_oldP = ta_newP;			// have selected NP mutants move on to next generation 
	ta_oldC = ta_newC;

	if (i_check_winner)  {			// check if the best stayed the same, if necessary 
	    int same = 1;
	    for (int j = 0; j < i_D; j++) {
		if (t_bestitP[j] != t_bestP[j]) {
		    same = 0;
		}
	    }
	    if (same && i_iter > 1)  {
		i_xav++;
		double tmp_best = evaluate(l_nfeval, t_bestP.memptr(), par, fcall, rho);	// if re-evaluation of winner 
		
		if (i_av_winner)		//  possibly letting the winner be the average of all past generations 
		    t_bestC = ((1/(double)i_xav) * t_bestC) + ((1/(double)i_xav) * tmp_best) + 
			(d_bestvalit[i_iter-1] * ((double)(i_xav - 2))/(double)i_xav);
		else
		    t_bestC = tmp_best;
	    } else {
		i_xav = 1;
	    }
	}
	t_bestitP = t_bestP;

	if ( (i_trace > 0)  &&  ((i_iter % i_trace) == 0) ) {
	    Rprintf("Iteration: %d bestvalit: %f bestmemit:", i_iter, t_bestC);
	    for (int j = 0; j < i_D; j++)
		Rprintf("%12.6f", t_bestP[j]);
	    Rprintf("\n");
	}
    } // end loop through generations 
    
    d_pop = ta_oldP;
    i_iterations = i_iter;

    PutRNGstate();   
    // ProfilerStop();
}
#endif
