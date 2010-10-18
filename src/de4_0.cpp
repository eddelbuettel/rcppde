// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// Port of DEoptim (2.0.7) to Rcpp/RcppArmadillo/Armadillo
// Copyright (C) 2010  Dirk Eddelbuettel <edd@debian.org>
//
// DEoptim is Copyright (C) 2009 David Ardia and Katharine Mullen
// and based on DE-Engine v4.0, Rainer Storn, 2004  
// (http://www.icsi.berkeley.edu/~storn/DeWin.zip)

#include <RcppArmadillo.h>
//#include <google/profiler.h>

RcppExport SEXP DEoptimC(SEXP lower, SEXP upper, SEXP fn, SEXP control, SEXP rho);
void devol(double VTR, double f_weight, double fcross, int i_bs_flag, 
           arma::rowvec & lower, arma::rowvec & upper, SEXP fcall, SEXP rho, int i_trace,
           int i_strategy, int i_D, int i_NP, int i_itermax,
           arma::mat & initpopm, int i_storepopfreq, int i_storepopfrom,
           int i_specinitialpop, int i_check_winner, int i_av_winner,
           arma::mat    & ta_popP, arma::mat    & ta_oldP, arma::mat    & ta_newP, arma::rowvec & t_bestP, 
	   arma::rowvec & ta_popC, arma::rowvec & ta_oldC, arma::rowvec & ta_newC, double       & t_bestC,	
           arma::rowvec & t_bestitP, arma::rowvec & t_tmpP, arma::rowvec & tempP,
           arma::mat & d_pop, Rcpp::List & d_storepop, arma::mat & d_bestmemit, arma::rowvec & d_bestvalit,
           int & i_iterations, double i_pPct, long & l_nfeval);
//RcppExport void permute(arma::icolvec & ia_urn2, int i_avoid, arma::icolvec & ia_urntmp);
//RcppExport double evaluate(long & l_nfeval, const arma::rowvec & param, SEXP par, SEXP fcall, SEXP env);
void permute(int ia_urn2[], int i_urn2_depth, int i_NP, int i_avoid, int ia_urntmp[]);
RcppExport double evaluate(long &l_nfeval, const double *param, SEXP parS, SEXP fcall, SEXP env);

RcppExport SEXP DEoptimC(SEXP lowerS, SEXP upperS, SEXP fnS, SEXP controlS, SEXP rhoS) {
    //ProfilerStart("/tmp/RcppDE.prof");
    BEGIN_RCPP ;	// macro to fill in try part of try/catch exception handler

    Rcpp::Function fn(fnS);						// function to mininise
    Rcpp::Environment rho(rhoS); 					// environment to do it in
    Rcpp::NumericVector f_lower(lowerS), f_upper(upperS); 		// User-defined bounds
    Rcpp::List          control(controlS); 				// named list of params

    double VTR           = Rcpp::as<double>(control["VTR"]);		// value to reach
    int i_strategy       = Rcpp::as<int>(control["strategy"]);    	// chooses DE-strategy
    int i_itermax        = Rcpp::as<int>(control["itermax"]);		// Maximum number of generations
    long l_nfeval        = 0;//Rcpp::as<int>(control["nfeval"]);	// number of function evaluations    
    int i_D              = Rcpp::as<int>(control["npar"]);		// Dimension of parameter vector
    int i_NP             = Rcpp::as<int>(control["NP"]);		// Number of population members
    int i_storepopfrom   = Rcpp::as<int>(control["storepopfrom"]) - 1;  // When to start storing populations 
    int i_storepopfreq   = Rcpp::as<int>(control["storepopfreq"]);  	// How often to store populations 
    int i_specinitialpop = Rcpp::as<int>(control["specinitialpop"]);  	// User-defined inital population 
    Rcpp::NumericMatrix initialpopm = Rcpp::as<Rcpp::NumericMatrix>(control["initialpop"]);
    double f_weight      = Rcpp::as<double>(control["F"]);  		// stepsize 
    double f_cross       = Rcpp::as<double>(control["CR"]);  		// crossover probability 
    int i_bs_flag        = Rcpp::as<int>(control["bs"]);   		// Best of parent and child 
    int i_trace          = Rcpp::as<int>(control["trace"]);  		// Print progress? 
    int i_check_winner   = Rcpp::as<int>(control["checkWinner"]); 	// Re-evaluate best parameter vector? 
    int i_av_winner      = Rcpp::as<int>(control["avWinner"]);  	// Average 
    double i_pPct        = Rcpp::as<double>(control["p"]); 		// p to define the top 100p% best solutions 

    arma::rowvec minbound(f_lower.begin(), f_lower.size(), false); 	// convert three Rcpp vectors to arma vectors
    arma::rowvec maxbound(f_upper.begin(), f_upper.size(), false);
    arma::mat initpopm(initialpopm.begin(), initialpopm.rows(), initialpopm.cols(), false);

    arma::mat ta_popP(i_NP*2, i_D);    					// Data structures for parameter vectors 
    arma::mat ta_oldP(i_NP,   i_D);
    arma::mat ta_newP(i_NP,   i_D);
    arma::rowvec t_bestP(i_D); 

    arma::rowvec ta_popC(i_NP*2);  				    	// Data structures for obj. fun. values associated with par. vectors 
    arma::rowvec ta_oldC(i_NP);
    arma::rowvec ta_newC(i_NP);
    double t_bestC; 

    arma::rowvec t_bestitP(i_D);
    arma::rowvec t_tmpP(i_D); 
    arma::rowvec tempP(i_D);

    int i_nstorepop = ceil((i_itermax - i_storepopfrom) / i_storepopfreq);
    arma::mat d_pop(i_NP, i_D); 
    Rcpp::List d_storepop(i_nstorepop);
    arma::mat d_bestmemit(i_itermax, i_D);       
    arma::rowvec d_bestvalit(i_itermax); 	 
    int i_iter = 0;

    // call actual Differential Evolution optimization given the parameters
    devol(VTR, f_weight, f_cross, i_bs_flag, minbound, maxbound, Rcpp::wrap(fn), Rcpp::wrap(rho), i_trace,
	  i_strategy, i_D, i_NP, i_itermax, initpopm, i_storepopfrom, i_storepopfreq, i_specinitialpop, i_check_winner, i_av_winner,
	  ta_popP, ta_oldP, ta_newP, t_bestP, ta_popC, ta_oldC, ta_newC, t_bestC, t_bestitP, t_tmpP, tempP,
	  d_pop, d_storepop, d_bestmemit, d_bestvalit, i_iter, i_pPct, l_nfeval);

    // and return a named list to R
    return Rcpp::List::create(Rcpp::Named("bestmem")   = t_bestP,	// sexp_bestmem,
			      Rcpp::Named("bestval")   = t_bestC,       // sexp_bestval,
			      Rcpp::Named("nfeval")    = l_nfeval,   	// sexp_nfeval,
			      Rcpp::Named("iter")      = i_iter,	// sexp_iter,
			      Rcpp::Named("bestmemit") = d_bestmemit, 	// sexp_bestmemit,
			      Rcpp::Named("bestvalit") = d_bestvalit,	// sexp_bestvalit,
			      Rcpp::Named("pop")       = d_pop,		// sexp_pop,
			      Rcpp::Named("storepop")  = d_storepop); 	// sexp_storepop)
    END_RCPP   // macro to fill in catch() part of try/catch exception handler
}


void devol(double VTR, double f_weight, double f_cross, int i_bs_flag,
           arma::rowvec & fa_minbound, arma::rowvec & fa_maxbound, SEXP fcall, SEXP rho, int i_trace,
           int i_strategy, int i_D, int i_NP, int i_itermax, arma::mat & initialpopm, 
	   int i_storepopfrom, int i_storepopfreq, int i_specinitialpop, int i_check_winner, int i_av_winner,
           arma::mat &ta_popP, arma::mat &ta_oldP, arma::mat &ta_newP, arma::rowvec & t_bestP, 
           arma::rowvec & ta_popC, arma::rowvec & ta_oldC, arma::rowvec & ta_newC, double & t_bestC,
           arma::rowvec & t_bestitP, arma::rowvec & t_tmpP, arma::rowvec & tempP,
           arma::mat &d_pop, Rcpp::List &d_storepop, arma::mat & d_bestmemit, arma::rowvec & d_bestvalit,
           int & i_iterations, double i_pPct, long & l_nfeval) {

    const int urn_depth = 5;   			// 4 + one index to avoid 
    //arma::rowvec parvec(i_D);  		// initialize parameter vector to pass to evaluate function 
    //SEXP par = Rcpp::wrap(parvec);
    Rcpp::NumericVector par(i_D);
    int i, j, k, i_r1, i_r2, i_r3, i_r4;  	// counting variables and placeholders for random indexes
    
    int ia_urn2[urn_depth];
    //Rcpp::IntegerVector ia_urntmp(i_NP); 	// so that we don't need to re-allocated each time in permute
    arma::irowvec ia_urntmp(i_NP);
    //arma::icolvec::fixed<urn_depth> ia_urn2; 	// fixed-size vector for urn draws
    //arma::icolvec ia_urntmp(i_NP); 		// so that we don't need to re-allocated each time in permute

    int i_nstorepop = ceil((i_itermax - i_storepopfrom) / i_storepopfreq);
    int i_xav, popcnt, bestacnt, same; 		// lazy cnters 
    double f_jitter, f_dither, t_bestitC, t_tmpC, tmp_best; 
    
    arma::mat initialpop(i_NP, i_D); 

    int i_pbest;    				// vars for DE/current-to-p-best/1 
    int p_NP = round(i_pPct * i_NP);  		// choose at least two best solutions 
    p_NP = p_NP < 2 ? 2 : p_NP;
    arma::irowvec sortIndex(i_NP); 		// sorted values of ta_oldC 
    for(i = 0; i < i_NP; i++) sortIndex[i] = i;

    int i_len, done, step, bound;    		// vars for when i_bs_flag == 1 */
    double tempC;

    GetRNGstate();

    initialpop.zeros();		 		// initialize initial popuplation 
    d_bestmemit.zeros();    			// initialize best members
    d_bestvalit.zeros();			// initialize best values 
    d_pop.zeros();				// initialize best population
    i_nstorepop = (i_nstorepop < 0) ? 0 : i_nstorepop;
      
    if (i_specinitialpop > 0) {    		// if initial population provided, initialize with values 
	initialpop = initialpopm;
    }
    //l_nfeval = 0;    				// already init'ed in main function

    for (i = 0; i < i_NP; i++) {		// ------Initialization-----------------------------
	if (i_specinitialpop <= 0) { 		// random initial member 
	    for (j = 0; j < i_D; j++) {
		ta_popP.at(i,j) = fa_minbound[j] + ::unif_rand() * (fa_maxbound[j] - fa_minbound[j]);
	    }
	} else { 				// or user-specified initial member 
	    ta_popP.row(i) = initialpop.row(i);
	} 
	arma::rowvec r = ta_popP.row(i);
	ta_popC[i] = evaluate(l_nfeval, r.memptr(), par, fcall, rho);
	if (i == 0 || ta_popC[i] <= t_bestC) {
	    t_bestC = ta_popC[i];
	    t_bestP = ta_popP.row(i);
	}
    }

    ta_oldP = ta_popP;				// ---assign pointers to current ("old") population---
    ta_oldC = ta_popC;
  
    int i_iter = 0;				// ------Iteration loop--------------------------------------------
    popcnt = 0;
    bestacnt = 0;
    i_xav = 1;
  
    while ((i_iter < i_itermax) && (t_bestC > VTR)) {    // main loop ====================================
	if (i_iter % i_storepopfreq == 0 && i_iter >= i_storepopfrom) {  	// store intermediate populations
	    d_storepop[popcnt++] = Rcpp::wrap(ta_oldP);
	} // end store pop 
      
	d_bestmemit.row(i_iter) = t_bestP;	// store the best member
	d_bestvalit[i_iter] = t_bestC;		// store the best value 
	t_bestitP = t_bestP;
	t_bestitC = t_bestC;
	i_iter++;				// increase iteration counter
     
	f_dither = f_weight + ::unif_rand() * (1.0 - f_weight);	// ----computer dithering factor -----------------
      
	if (i_strategy == 6) {			// ---DE/current-to-p-best/1 -----------------------------------------------------
	    arma::rowvec temp_oldC = ta_oldC;					// create a copy of ta_oldC to avoid changing it 
	    rsort_with_index( temp_oldC.memptr(), sortIndex.begin(), i_NP );  	// sort temp_oldC to use sortIndex later 
	}

	for (i = 0; i < i_NP; i++) {		// ----start of loop through ensemble------------------------

	    t_tmpP = ta_oldP.row(i);		// t_tmpP is the vector to mutate and eventually select
	    t_tmpC = ta_oldC[i];

	    permute(ia_urn2, urn_depth, i_NP, i, ia_urntmp.memptr()); // Pick 4 random and distinct 
	    //permute(ia_urn2, i, ia_urntmp); 	// Pick 4 random and distinct 
	    i_r1 = ia_urn2[1];  		// population members 
	    i_r2 = ia_urn2[2];
	    i_r3 = ia_urn2[3];
	    i_r4 = ia_urn2[4];
	    k = 0;				// loop counter used in all strategies below 

	    // ===Choice of strategy=======================================================
	    switch (i_strategy) { 		// and putting default value one first

	    case 2:				// ---DE/local-to-best/1/bin---------------------------------------------------
		j = static_cast<int>(::unif_rand() * i_D); 	// random parameter 
		do {				// add fluctuation to random target 
		    t_tmpP[j] = t_tmpP[j] + f_weight * (t_bestitP[j] - t_tmpP[j]) + f_weight * (ta_oldP.at(i_r2,j) - ta_oldP.at(i_r3,j));
		    j = (j + 1) % i_D;
		} while ((::unif_rand() < f_cross) && (++k < i_D));
		break;

	    case 1:				// ---classical strategy DE/rand/1/bin-----------------------------------------
		j = static_cast<int>(::unif_rand() * i_D); 	// random parameter 
		do {				// add fluctuation to random target 
		    t_tmpP[j] = ta_oldP.at(i_r1,j) + f_weight * (ta_oldP.at(i_r2,j) - ta_oldP.at(i_r3,j));
		    j = (j + 1) % i_D;
		} while ((::unif_rand() < f_cross) && (++k < i_D));
		break;

	    case 3:				// ---DE/best/1/bin with jitter------------------------------------------------
		j = static_cast<int>(::unif_rand() * i_D); 	// random parameter 
		do {				// add fluctuation to random target 
		    f_jitter = 0.0001 * ::unif_rand() + f_weight; 
		    t_tmpP[j] = t_bestitP[j] + f_jitter * (ta_oldP.at(i_r1,j) - ta_oldP.at(i_r2,j));
		    j = (j + 1) % i_D;
		} while ((::unif_rand() < f_cross) && (++k < i_D));
		break;

	    case 4:				// ---DE/rand/1/bin with per-vector-dither-------------------------------------
		j = static_cast<int>(::unif_rand() * i_D); 	// random parameter 
		do {				// add fluctuation to random target *
		    t_tmpP[j] = ta_oldP.at(i_r1,j) + (f_weight + ::unif_rand()*(1.0 - f_weight))* (ta_oldP.at(i_r2,j) - ta_oldP.at(i_r3,j));
		    j = (j + 1) % i_D;
		} while ((::unif_rand() < f_cross) && (++k < i_D));
		break;

	    case 5:				// ---DE/rand/1/bin with per-generation-dither---------------------------------
		j = static_cast<int>(::unif_rand() * i_D); 	// random parameter 
		do {				// add fluctuation to random target 
		    t_tmpP[j] = ta_oldP.at(i_r1,j) + f_dither * (ta_oldP.at(i_r2,j) - ta_oldP.at(i_r3,j));
		    j = (j + 1) % i_D;
		} while ((::unif_rand() < f_cross) && (++k < i_D));
		break;

	    case 6:				// ---DE/current-to-p-best/1 (JADE)--------------------------------------------
		i_pbest = sortIndex[static_cast<int>(::unif_rand() * p_NP)]; // select from [0, 1, 2, ..., (pNP-1)] 
		j = static_cast<int>(::unif_rand() * i_D); 	// random parameter 
		do {				// add fluctuation to random target 
		    t_tmpP[j] = ta_oldP.at(i,j) + f_weight * (ta_oldP.at(i_pbest,j) - ta_oldP.at(i,j)) + f_weight * (ta_oldP.at(i_r1,j) - ta_oldP.at(i_r2,j));
		    j = (j + 1) % i_D;
		} while ((::unif_rand() < f_cross) && (++k < i_D));
		break;

	    default:				// ---variation to DE/rand/1/bin: either-or-algorithm--------------------------
		j = static_cast<int>(::unif_rand() * i_D); 	// random parameter 
		if (::unif_rand() < 0.5) { 	// differential mutation, Pmu = 0.5 
		    do {
			// add fluctuation to random target */
			t_tmpP[j] = ta_oldP.at(i_r1,j) + f_weight * (ta_oldP.at(i_r2,j) - ta_oldP.at(i_r3,j));
			j = (j + 1) % i_D;
		    } while ((::unif_rand() < f_cross) && (++k < i_D));

		} else { 			// recombination with K = 0.5*(F+1) -. F-K-Rule 
		    do {
			// add fluctuation to random target */
			t_tmpP[j] = ta_oldP.at(i_r1,j) + 0.5 * (f_weight + 1.0) * (ta_oldP.at(i_r2,j) + ta_oldP.at(i_r3,j) - 2 * ta_oldP.at(i_r1,j));
			j = (j + 1) % i_D;
		    } while ((::unif_rand() < f_cross) && (++k < i_D));
		}
		break;
	    } // end switch (i_strategy) ...
	
	    for (j = 0; j < i_D; j++) {		// ----boundary constraints, bounce-back method was not enforcing bounds correctly
		if (t_tmpP[j] < fa_minbound[j]) {
		    t_tmpP[j] = fa_minbound[j] + ::unif_rand() * (fa_maxbound[j] - fa_minbound[j]);
		}
		if (t_tmpP[j] > fa_maxbound[j]) {
		    t_tmpP[j] = fa_maxbound[j] - ::unif_rand() * (fa_maxbound[j] - fa_minbound[j]);
		}
	    }

	    // ------Trial mutation now in t_tmpP-----------------
	    t_tmpC = evaluate(l_nfeval, t_tmpP.memptr(), par, fcall, rho);	// Evaluate mutant in t_tmpP[]

	    // note that i_bs_flag means that we will choose the best NP vectors from the old and new population later
	    if (t_tmpC <= ta_oldC[i] || i_bs_flag) {
		ta_newP.row(i) = t_tmpP;				// replace target with mutant 
		ta_newC[i] = t_tmpC;
		if (t_tmpC <= t_bestC) {
		    t_bestP = t_tmpP;
		    t_bestC = t_tmpC;
		}
	    } else {
		ta_newP.row(i) = ta_oldP.row(i);
		ta_newC[i] = ta_oldC[i];
	    }
	} // End mutation loop through pop., ie the "for (i = 0; i < i_NP; i++)"

	if (i_bs_flag) {	// examine old and new pop. and take the best NP members into next generation 
	    
	    ta_popP.rows(0, i_NP-1) = ta_oldP;
	    ta_popC.rows(0, i_NP-1) = ta_oldC;

	    ta_popP.rows(i_NP, 2*i_NP-1) = ta_newP;
	    ta_popC.rows(i_NP, 2*i_NP-1) = ta_newC;

	    i_len = 2 * i_NP;
	    step = i_len;  		// array length 
	    while (step > 1) {
		step /= 2;   		// halve the step size 
		do {
		    done = 1;
		    bound  = i_len - step;
		    for (j = 0; j < bound; j++) {
			i = j + step + 1;
			if (ta_popC[j] > ta_popC[i-1]) {
			    tempP = ta_popP.row(i-1);
			    tempC = ta_popC[i-1];
			    ta_popP.row(i-1) = ta_popP.row(j);
			    ta_popC[i-1] = ta_popC[j];
			    ta_popP.row(j) = tempP;
			    ta_popC[j] = tempC;
			    done = 0; 
			}  // if 
		    }  // for 
		} while (!done);   // do .. while 
	    } // while (step > 1) 
	    ta_newP = ta_popP;			// now the best NP are in first NP places in gta_pop, use them
	    ta_newC = ta_popC;
	} // i_bs_flag

	ta_oldP = ta_newP;			// have selected NP mutants move on to next generation 
	ta_oldC = ta_newC;

	if (i_check_winner)  {			// check if the best stayed the same, if necessary 
	    same = 1;
	    for (j = 0; j < i_D; j++)
		if (t_bestitP[j] != t_bestP[j]) {
		    same = 0;
		}
	    if (same && i_iter > 1)  {
		i_xav++;
		tmp_best = evaluate(l_nfeval, t_bestP.memptr(), par, fcall, rho);			// if re-evaluation of winner 
		
		if (i_av_winner)		//  possibly letting the winner be the average of all past generations 
		    t_bestC = ((1/(double)i_xav) * t_bestC) + ((1/(double)i_xav) * tmp_best) + (d_bestvalit[i_iter-1] * ((double)(i_xav - 2))/(double)i_xav);
		else
		    t_bestC = tmp_best;
	    } else {
		i_xav = 1;
	    }
	}
	t_bestitP = t_bestP;
	t_bestitC = t_bestC;

	if ( (i_trace > 0)  &&  ((i_iter % i_trace) == 0) ) {
	    Rprintf("Iteration: %d bestvalit: %f bestmemit:", i_iter, t_bestC);
	    for (j = 0; j < i_D; j++)
		Rprintf("%12.6f", t_bestP[j]);
	    Rprintf("\n");
	}
    } // end loop through generations 
    
    d_pop = ta_oldP;
    i_iterations = i_iter;

    PutRNGstate();
    //ProfilerStop();
}

inline void permute(int ia_urn2[], int i_urn2_depth, int i_NP, int i_avoid, int ia_urn1[])
//RcppExport inline void permute(arma::icolvec & ia_urn2, int i_avoid, arma::icolvec & ia_urn1)
/********************************************************************
 ** Function       : void permute(int ia_urn2[], int i_urn2_depth)
 ** Author         : Rainer Storn (w/bug fixes contributed by DEoptim users)
 ** Description    : Generates i_urn2_depth random indices ex [0, i_NP-1]
 **                  which are all distinct. This is done by using a
 **                  permutation algorithm called the "urn algorithm"
 **                  which goes back to C.L.Robinson.
 ** Functions      : -
 ** Globals        : -
 ** Parameters     : ia_urn2       (O)    array containing the random indices
 **                  i_urn2_depth  (I)    number of random indices (avoided index included)
 **                  i_NP          (I)    range of indices is [0, i_NP-1]
 **                  i_avoid       (I)    is the index to avoid and is located in
 **                                       ia_urn2[0].
 ** Preconditions  : # Make sure that ia_urn2[] has a length of i_urn2_depth.
 **                  # i_urn2_depth must be smaller than i_NP.
 ** Postconditions : # the index to be avoided is in ia_urn2[0], so fetch the
 **                   indices from ia_urn2[i], i = 1, 2, 3, ..., i_urn2_depth.
 ** Return Value   : -
 *********************************************************************/
{
    GetRNGstate();
    int k = i_NP;
    //int i_urn2_depth = ia_urn2.n_elem; 		// ie const urn_depth
    //int i_NP = ia_urn1.n_elem;			
    //int k = ia_urn1.n_elem;			
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
