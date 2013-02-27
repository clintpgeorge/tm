#ifndef LDAHP_H
#define LDAHP_H

#include <RcppArmadillo.h>

using namespace Rcpp ;
using namespace arma ;
using namespace std ;

RNGScope scope;

unsigned int sample_multinomial (vec theta);
vec sample_dirichlet (unsigned int num_elements, vec alpha);
rowvec sample_dirichlet_row_vec (unsigned int num_elements, rowvec alpha);

RcppExport SEXP lda_full(SEXP K_, SEXP V_, SEXP doc_lengths_, SEXP docs_,
		SEXP zid_, SEXP alpha_v_, SEXP eta_,
		SEXP max_iter_, SEXP burn_in_, SEXP spacing_);

RcppExport SEXP lda_full2(SEXP K, SEXP V, SEXP docN, SEXP wid,
		SEXP zid, SEXP alphav, SEXP etac,
		SEXP maxiter, SEXP burnin, SEXP spacing);

RcppExport SEXP lda_z_theta_fixed_beta(SEXP docN, SEXP wid,
		SEXP zid, SEXP alphav, SEXP beta,
		SEXP maxiter, SEXP burnin, SEXP spacing);

#endif
