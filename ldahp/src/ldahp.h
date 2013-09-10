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

RcppExport SEXP lda_full(SEXP num_topics_, SEXP vocab_size_, SEXP doc_lengths_, SEXP docs_,
  SEXP topic_assignments_, SEXP alpha_v_, SEXP eta_,
	SEXP max_iter_, SEXP burn_in_, SEXP spacing_, SEXP store_dirichlet_);

RcppExport SEXP lda_full2(SEXP num_topics_, SEXP vocab_size_,
    SEXP doc_lengths_, SEXP word_ids_,
  SEXP topic_assignments_, SEXP alpha_v_, SEXP eta_,
	SEXP max_iter_, SEXP burn_in_, SEXP spacing_, SEXP store_dirichlet_);

RcppExport SEXP lda_z_theta_fixed_beta(SEXP doc_lengths_, SEXP word_ids_,
  	SEXP topic_assignments_, SEXP alpha_v_, SEXP beta_,
		SEXP max_iter_, SEXP burn_in_, SEXP spacing_, SEXP store_dirichlet_) ;
    
RcppExport SEXP lda_collapsed_gibbs(SEXP num_topics_, SEXP vocab_size_,
  	SEXP doc_lengths_, SEXP word_ids_,
		SEXP topic_assignments_, SEXP alpha_v_, SEXP eta_,
		SEXP max_iter_, SEXP burn_in_, SEXP spacing_, SEXP store_dirichlet_);

#endif
