#include "ldahp.h"

/**
 * Samples an integer from [1, K] uniformly at random
 * 
 * Arguments: 
 * 		K - the upper interval 
 * Returns:
 * 		the sampled integer  
 */
size_t sample_uniform_int(size_t K){
	return (size_t) (runif(1)(0) * (double)K); // To speedup
}

/**
 * Samples from a given multimomial probability vector
 *
 * Arguments:
 * 		theta - the Multinomial probability vector
 * Returns:
 * 		t - the sampled index
 *
 */
unsigned int sample_multinomial (vec theta) {

	size_t t = 0;
	double total_prob = accu(theta);
	double u = runif(1)(0) * total_prob;
	double cumulative_prob = theta(0);

	while(u > cumulative_prob){
		t++;
		cumulative_prob += theta(t);
	}

	return t;

}

/**
 * Samples from a Dirichlet distribution from a 
 * given set of hyperparameters
 * 
 * Aruguments:
 * 		num_elements - the dimentionality of the Dirichlet distribution 
 * 		alpha - the hyperparameter vector which is in the column vector format  
 * Returns: 
 * 		the Dirichlet sample in the column vector format   
 */
vec sample_dirichlet(unsigned int num_elements, vec alpha){

	vec dirichlet_sample = zeros<vec>(num_elements);

	for ( register unsigned int i = 0; i < num_elements; i++ )
		dirichlet_sample(i) = rgamma(1, alpha(i), 1.0)(0);

	dirichlet_sample /= accu(dirichlet_sample);

	return dirichlet_sample;

}

/**
 * Samples from a Dirichlet distribution from a 
 * given set of hyperparameters
 * 
 * Aruguments:
 * 		num_elements - the dimentionality of the Dirichlet distribution 
 * 		alpha - the hyperparameter vector which is in the row vector format  
 * Returns: 
 * 		the Dirichlet sample in the row vector format   
 */
rowvec sample_dirichlet_row_vec (unsigned int num_elements, rowvec alpha){

	rowvec dirichlet_sample = zeros<rowvec>(num_elements);

	for ( register unsigned int i = 0; i < num_elements; i++ )
		dirichlet_sample(i) = rgamma(1, alpha(i), 1.0)(0);

	dirichlet_sample /= accu(dirichlet_sample);

	return dirichlet_sample;

}

/**
 * Samples random permutations for a given count
 *
 * Arguments:
 * 		n - the number of samples
 * Return:
 * 		order - a vector of indices that represents
 * 				the permutations of numbers in [1, n]
 **/
uvec randperm(int n) //
{
	uvec order = zeros<uvec>(n);
	int k, nn, takeanumber, temp;
	for (k=0; k<n; k++) order(k) = k;
	nn = n;
	for (k=0; k<n; k++) {
		takeanumber = sample_uniform_int(nn); // take a number between 0 and nn-1
		temp = order(nn-1);
		order(nn-1) = order(takeanumber);
		order(takeanumber) = temp;
		nn--;
	}
	return order;
}


vec log_gamma_vec(vec x_vec){

	vec lgamma_vec = zeros<vec>(x_vec.n_elem);

	for (size_t i = 0; i < x_vec.n_elem; i++)
		lgamma_vec(i) = lgamma(x_vec(i));

	return lgamma_vec;

}

/*
 * This function calculates log marginal posterior
 * of the LDA model for a given corpus
 *
 */
double calc_log_marginal_posterior(
		int num_topics,
		int num_docs,
		int vocab_size,
		mat theta_counts,
		mat beta_counts) {

	double lmp = 0.0;
	vec ln_gamma_Nj = zeros<vec> (num_topics);
	vec ln_gamma_Mj = zeros<vec> (num_topics);

	// \sum_{d = 0}^{D} ln_gamma (n_j,d + alpha_j)
	for (size_t d = 0; d < num_docs; d++)
		ln_gamma_Nj += log_gamma_vec(theta_counts.col(d)); // ignores the smoothing alpha_j, because it's already added to \theta 

	// \sum_{t = 0}^{V} ln_gamma (\sum_{d=0}^{D}m_{jt,d} + eta)
	for (size_t t = 0; t < vocab_size; t++)
		ln_gamma_Mj += log_gamma_vec(beta_counts.col(t)); // ignores the smoothing eta, because it's already added to \beta 

	lmp = accu(ln_gamma_Nj + ln_gamma_Mj); // sum over all topics

	return lmp;

}





/**
 *  The LDA full gibbs sampler:
 * 		This function samples z, \theta, and \beta. The only difference
 * 		between the succeeding function and this one is that the input
 * 		format of the corpus documents. This function assumes that the
 * 		documents are read from Blei corpus format. In addition,
 * 		this function is optimized for \beta (topic) sampling.
 *
 * 	Arguments:
 * 		num_topics_ - number of topics
 * 		vocab_size_ - vocabulary size
 * 		doc_lengths_ - contains the number of words in each document as a vector
 * 		docs_ - documents that are read from a Blei corpus and in list format
 * 		topic_assignments_ - initial topic assignments for each word in the corpus
 * 		alpha_v_ - the hyperparameter vector for the document topic Dirichlet
 * 		eta_ - the hyperparameter for the topic Dirichlet
 * 		max_iter_ - max number of Gibbs iterations
 * 		burn_in_ - burn in period
 * 		spacing_ - spacing between samples that are stored
 *
 * 	Returns:
 * 		thetas - the sampled thetas after burn in period
 * 		betas - the sampled betas after burn in period
 * 		Z - the sampled word topic assignments after burn in period
 *
 */
RcppExport SEXP lda_full(SEXP num_topics_, SEXP vocab_size_, SEXP doc_lengths_, SEXP docs_,
		SEXP topic_assignments_, SEXP alpha_v_, SEXP eta_,
		SEXP max_iter_, SEXP burn_in_, SEXP spacing_, SEXP store_dirichlet_) {

	// Variable from the R interface  

	uvec z = as<uvec>(topic_assignments_);
	vec alpha_v = as<vec>(alpha_v_);
	uvec doc_lengths = as<uvec>(doc_lengths_);

	double eta = as<double>(eta_);
	int num_topics = as<int>(num_topics_);
	int vocab_size = as<int>(vocab_size_);
	int max_iter = as<int>(max_iter_);
	int burn_in = as<int>(burn_in_);
	int spacing = as<int>(spacing_);
	int store_dirichlet = as<int>(store_dirichlet_);

	// Function variables 

	int num_docs = doc_lengths.n_elem;
	int num_word_instances = accu(doc_lengths);
	int valid_samples = ceil((max_iter - burn_in) / (double) spacing);
	cube thetas;
	cube betas;
	if (store_dirichlet == 1){
		thetas = cube(num_topics, num_docs, valid_samples);
		betas = cube(num_topics, vocab_size, valid_samples);
	}
	umat Z = zeros<umat>(num_word_instances, valid_samples);
	vec log_marginal = zeros<vec>(valid_samples);

	mat prior_beta_samples = zeros<mat>(num_topics, vocab_size);
	mat beta_counts = zeros<mat>(num_topics, vocab_size);
	unsigned int d, i, k, iter, count = 0, instances = 0, c = 0, word_id, word_count, num_unique_words;
	uvec word_ids = zeros<uvec>(num_word_instances);
	vector < vector < unsigned int > > doc_unique_words, doc_word_indices;
	vector < unsigned int > unique_words, word_idx;


	// Calculates the document word indices

	cout << "Loading documents....";

	for (d = 0; d < num_docs; d++){

		umat document = as<umat>(VECTOR_ELT(docs_, d));
		vector < unsigned int > word_idx;
		vector < unsigned int > uw;

		for (c = 0; c < document.n_cols; c++){
			word_id = document(0,c);
			word_count = document(1,c);
			for (i = 0; i < word_count; i++){
				word_ids(instances) = word_id;
				word_idx.push_back(instances);
				instances++;
			}
			uw.push_back(word_id);
		}
		doc_word_indices.push_back(word_idx);
		doc_unique_words.push_back(uw);

	}

	cout << "DONE." << endl;

	// Initilizes beta

	prior_beta_samples.fill(1e-10);
	beta_counts.fill(eta); // initializes with the smoothing parameter
	for (i = 0; i < num_word_instances; i++)
		beta_counts(z(i), word_ids(i)) += 1;
	mat prior_beta_counts = beta_counts;

	// The Gibbs sampling loop

	cout << "Gibbs sampling..." << endl;

	for (iter = 0; iter < max_iter; iter++){

		if (iter % 1000 == 0)
			cout << "gibbs iter# " << iter + 1;

		mat prior_theta_samples = zeros<mat>(num_topics, num_docs);
		mat prior_theta_counts = zeros<mat>(num_topics, num_docs);

		for (d = 0; d < num_docs; d++){ // for each document 

			word_idx = doc_word_indices[d];
			unique_words = doc_unique_words[d];
			num_unique_words = unique_words.size();
			mat doc_beta_samples = zeros<mat>(num_topics, num_unique_words + 1);

			// samples \beta (with selected samples)
			for (i = 0; i < num_unique_words; i++)
				doc_beta_samples.col(i) = beta_counts.col(unique_words[i]);
			doc_beta_samples.col(num_unique_words) = sum(beta_counts, 1) - sum(doc_beta_samples, 1); // to avoid error

			for(k = 0; k < num_topics; k++)
				doc_beta_samples.row(k) = sample_dirichlet_row_vec(num_unique_words + 1, doc_beta_samples.row(k));

			for (i = 0; i < num_unique_words; i++)
				prior_beta_samples.col(unique_words[i]) = doc_beta_samples.col(i);


			// samples \theta
			vec partition_counts = alpha_v; // initializes with the smoothing parameter
			for (i = 0; i < doc_lengths(d); i++)
				partition_counts(z(word_idx[i])) += 1;
			vec theta_d = sample_dirichlet(num_topics, partition_counts);

			prior_theta_samples.col(d) = theta_d;
			prior_theta_counts.col(d) = partition_counts;

			// samples z
			for(i = 0; i < doc_lengths(d); i++)
				beta_counts(z(word_idx[i]), word_ids(word_idx[i])) -= 1; // excludes document d's word-topic counts
			for (i = 0; i < doc_lengths(d); i++)
				z(word_idx[i]) = sample_multinomial(theta_d % prior_beta_samples.col(word_ids(word_idx[i])));
			for(i = 0; i < doc_lengths(d); i++)
				beta_counts(z(word_idx[i]), word_ids(word_idx[i])) += 1; // includes document d's word-topic counts


		}

		if ((iter >= burn_in) && (iter % spacing == 0)){ // handles burn in period 
			Z.col(count) = z;
			if (store_dirichlet == 1){
				thetas.slice(count) = prior_theta_samples;
				betas.slice(count) = prior_beta_samples;
			}
			log_marginal(count) = calc_log_marginal_posterior(num_topics, num_docs, vocab_size, prior_theta_counts, prior_beta_counts); 
			if (iter % 1000 == 0)
				cout << " lmp: " << log_marginal(count);

			count++;
		}

		prior_beta_counts = beta_counts;
		if (iter % 1000 == 0)
			cout << endl;

	}

	cout << "END Gibbs sampling..." << endl;

	if (store_dirichlet == 1){
		return List::create(
				Named("thetas") = wrap(thetas),
				Named("betas") = wrap(betas),
				Named("Z") = wrap(Z),
				Named("lmp") = wrap(log_marginal));
	}
	else {
		return List::create(
				Named("Z") = wrap(Z),
				Named("lmp") = wrap(log_marginal));
	}


}

/**
 *  The LDA full gibbs sampler:
 * 		This function samples z, \theta, and \beta.
 * 		This function assumes that the document word
 * 		instances are in Gibbs sampling format.
 *
 * 	Arguments:
 * 		num_topics_ - number of topics
 * 		vocab_size_ - vocabulary size
 * 		doc_lengths_ - a vector that contains the number of words in each document
 * 		word_ids_ - a vector of document word instances
 * 		topic_assignments_ - a vector of initial topic assignments for each word in the corpus
 * 		alpha_v_ - a vector of hyperparameters of the document Dirichlet
 * 		eta_ - the hyperparameter of the topic Dirichlet
 * 		max_iter_ - the max number of Gibbs iterations to be performed
 * 		burn_in_ - burn in period
 * 		spacing_ - the spacing between samples that are stored
 *
 * 	Returns:
 * 		thetas - the sampled thetas after burn in period
 * 		betas - the sampled betas after burn in period
 * 		Z - the sampled word topic assignments after burn in period
 *
 */
RcppExport SEXP lda_full2(SEXP num_topics_, SEXP vocab_size_, 
		SEXP doc_lengths_, SEXP word_ids_,
		SEXP topic_assignments_, SEXP alpha_v_, SEXP eta_,
		SEXP max_iter_, SEXP burn_in_, SEXP spacing_, SEXP store_dirichlet_) {

	// Variable from the R interface  

	uvec doc_lengths = as<uvec>(doc_lengths_);
	uvec word_ids = as<uvec>(word_ids_);
	uvec z = as<uvec>(topic_assignments_);
	vec alpha_v = as<vec>(alpha_v_);

	double eta = as<double>(eta_);
	int num_topics = as<int>(num_topics_);
	int vocab_size = as<int>(vocab_size_);
	int max_iter = as<int>(max_iter_);
	int burn_in = as<int>(burn_in_);
	int spacing = as<int>(spacing_);
	int store_dirichlet = as<int>(store_dirichlet_);

	// Function variables 

	int num_docs = doc_lengths.n_elem;
	int num_word_instances = word_ids.n_elem;
	int valid_samples = ceil((max_iter - burn_in) / (double) spacing);
  
  // cout << "number of saved samples: " << valid_samples << endl; 
  
	cube thetas;
	cube betas;
	if (store_dirichlet == 1){
		thetas = cube(num_topics, num_docs, valid_samples);
		betas = cube(num_topics, vocab_size, valid_samples);
	}
	umat Z = zeros<umat>(num_word_instances, valid_samples);
	vec log_marginal = zeros<vec>(valid_samples);

	mat prior_beta_samples = zeros<mat>(num_topics, vocab_size);
	mat beta_counts = zeros<mat>(num_topics, vocab_size);
	vector < vector < size_t > > document_word_indices;
	unsigned int d, i, k, iter, count = 0, instances = 0;
	rowvec eta_v = zeros<rowvec>(vocab_size);
	for(k = 0; k < vocab_size; k++)
		eta_v(k) = eta;

	// Calculates the indices for each word in a document as a vector
	for (d = 0; d < num_docs; d++){
		vector < size_t > word_idx;
		for (i = 0; i < doc_lengths(d); i++){
			word_idx.push_back(instances);
			instances++;
		}
		document_word_indices.push_back(word_idx);
	}

	// Initilizes beta

	beta_counts.fill(eta); // initializes with the smoothing parameter
	for (i = 0; i < num_word_instances; i++)
		beta_counts(z(i), word_ids(i)) += 1;

	// The Gibbs sampling loop

	uvec prior_z = z;
	mat prior_beta_counts = beta_counts;

	for (iter = 0; iter < max_iter; iter++){ 

		if (iter % 1000 == 0)
			cout << "gibbs iter# " << iter + 1;

		mat prior_theta_samples = zeros<mat>(num_topics, num_docs);
		mat prior_theta_counts = zeros <mat>(num_topics, num_docs);

		for (d = 0; d < num_docs; d++){ // for each document

			vector < size_t > word_idx = document_word_indices[d];

			// samples \beta
			for(k = 0; k < num_topics; k++)
				prior_beta_samples.row(k) = sample_dirichlet_row_vec(vocab_size, beta_counts.row(k));

			// samples \theta
			vec partition_counts = alpha_v; // initializes with the smoothing parameter
			for (i = 0; i < doc_lengths(d); i++)
				partition_counts(z(word_idx[i])) += 1;
			vec theta_d = sample_dirichlet(num_topics, partition_counts);
			prior_theta_samples.col(d) = theta_d;
			prior_theta_counts.col(d) = partition_counts;


			// samples z
			for(i = 0; i < doc_lengths(d); i++)
				beta_counts(z(word_idx[i]), word_ids(word_idx[i])) -= 1; // excludes document d's word-topic counts
			for (i = 0; i < doc_lengths(d); i++)
				z(word_idx[i]) = sample_multinomial(theta_d % prior_beta_samples.col(word_ids(word_idx[i])));
			for(i = 0; i < doc_lengths(d); i++)
				beta_counts(z(word_idx[i]), word_ids(word_idx[i])) += 1; // includes document d's word-topic counts

		}

		if ((iter >= burn_in) && (iter % spacing == 0)){ // Handles burn in period 
			Z.col(count) = prior_z;
			if (store_dirichlet == 1){
				thetas.slice(count) = prior_theta_samples; // theta_counts still has the counts from old z
				betas.slice(count) = prior_beta_samples;
			}

			log_marginal(count) = calc_log_marginal_posterior(num_topics, num_docs, vocab_size, prior_theta_counts, prior_beta_counts);
			if (iter % 1000 == 0)
				cout << " lmp: " << log_marginal(count);

			count++;
		}

		prior_z = z;
		prior_beta_counts = beta_counts;
		if ((iter+1) % 1000 == 0)
			cout << endl;

	} // The end of the Gibbs loop
  
  cout << "number of saved samples: " << count << endl; 
  
	if (store_dirichlet == 1){
		return List::create(
				Named("thetas") = wrap(thetas),
				Named("betas") = wrap(betas),
				Named("Z") = wrap(Z),
				Named("lmp") = wrap(log_marginal));
	}
	else {
		return List::create(
				Named("Z") = wrap(Z),
				Named("lmp") = wrap(log_marginal));
	}


}


/**
 *  The LDA collapsed Gibbs sampler:
 * 		This function samples only z. This function assumes
 * 		that the document word instances are in Gibbs sampling format.
 *
 *	References:
 *		1. Finding scientific topics by Griffiths and Steyvers
 *		2. LDA collapsed Gibbs sampler implementation by David Newman
 *
 *
 * 	Arguments:
 * 		num_topics_ - number of topics
 * 		vocab_size_ - vocabulary size
 * 		doc_lengths_ - contains the number of words in each document as a vector
 * 		word_ids_ - the document word instances
 * 		topic_assignments_ - initial topic assignments for each word in the corpus
 * 		alpha_v_ - the hyperparameter vector for the document topic Dirichlet
 * 		eta_ - the hyperparameter for the topic Dirichlet
 * 		max_iter_ - max number of Gibbs iterations
 * 		burn_in_ - burn in period
 * 		spacing_ - spacing between samples that are stored
 *
 * 	Returns:
 * 		thetas - the sampled thetas after burn in period
 * 		betas - the sampled betas after burn in period
 * 		Z - the sampled word topic assignments after burn in period
 *
 */
RcppExport SEXP lda_collapsed_gibbs(SEXP num_topics_, SEXP vocab_size_,
		SEXP doc_lengths_, SEXP word_ids_,
		SEXP topic_assignments_, SEXP alpha_v_, SEXP eta_,
		SEXP max_iter_, SEXP burn_in_, SEXP spacing_, SEXP store_dirichlet_) {

	// variable declarations

	uvec doc_lengths = as<uvec>(doc_lengths_); // number of words in each document
	uvec word_ids = as<uvec>(word_ids_); // word indices
	uvec z = as<uvec>(topic_assignments_); // we get this because, we wanna use a given starting point for Gibbs
	vec alpha_v = as<vec>(alpha_v_); // hyperparameters for the document Dirichlets

	double eta = as<double>(eta_); // hyperparameter for the topic Dirichlets
	int num_topics = as<int>(num_topics_);
	int vocab_size = as<int>(vocab_size_);
	int max_iter = as<int>(max_iter_);
	int burn_in = as<int>(burn_in_);
	int spacing = as<int>(spacing_);
	int store_dirichlet = as<int>(store_dirichlet_);
	int num_docs = doc_lengths.n_elem; // number of documents in the corpus
	int num_word_instances = word_ids.n_elem; // total number of words in the corpus
	int valid_samples = ceil((max_iter - burn_in) / (double) spacing);
	cube thetas;
	cube betas;
	if (store_dirichlet == 1){
		thetas = cube(num_topics, num_docs, valid_samples);
		betas = cube(num_topics, vocab_size, valid_samples);
	}
	umat Z = zeros<umat>(num_word_instances, valid_samples);
	vec log_marginal = zeros<vec>(valid_samples);

	mat beta_counts = zeros<mat>(num_topics, vocab_size);
	mat theta_counts = zeros <mat>(num_topics, num_docs);
	vec topic_counts = zeros<vec> (num_topics);
	uvec doc_ids = zeros<uvec>(num_word_instances);
	unsigned int d, i, k, iter, count = 0, instances = 0, wid, did, topic, new_topic, idx;
	double doc_denom;
	vec prob;

	// Gets a random permutation of indices
	// this may improve mixing
	// Reference: Dr. Newman's implementation
	// uvec porder = randperm (num_word_instances);


	// Gets the document index for each word instance

	for (d = 0; d < num_docs; d++){
		for (i = 0; i < doc_lengths(d); i++){
			doc_ids(instances) = d;
			instances++;
		}
	}


	// Initializes beta, theta, and topic counts

	beta_counts.fill(eta); // initializes with the smoothing parameter
	for (d = 0; d < num_docs; d++)
		theta_counts.col(d) = alpha_v; // initializes with the smoothing parameter

	for (i = 0; i < num_word_instances; i++){
		beta_counts(z(i), word_ids(i)) += 1;
		topic_counts (z(i)) += 1;
		theta_counts (z(i), doc_ids(i)) += 1;
	}


	// The Gibbs sampling loop

	for (iter = 0; iter < max_iter; iter++){ // for each Gibbs iteration

		if (iter % 1000 == 0)
			cout << "gibbs iter# " << iter + 1;

		for (i = 0; i < num_word_instances; i++){ // for each word instance

			idx = i; //  porder(i); // permutation
			wid = word_ids(idx); // word index
			did = doc_ids(idx); // document index
			topic = z(idx); // old topic
			prob = zeros <vec> (num_topics); // init. probability vector

			// decrements the counts by one, to ignore the current sampling word

			beta_counts(topic, wid)--;
			theta_counts(topic, did)--;
			topic_counts(topic)--;

			doc_denom = doc_lengths(did) - 1 + accu(alpha_v); // a constant for a term

			for (k = 0; k < num_topics; k++){ // for each topic compute P(z_i == j | z_{-i}, w)
				prob(k) = (theta_counts(k, did) / doc_denom) *  (beta_counts(k, wid) /  (topic_counts(k) + vocab_size * eta));
			}
			new_topic = sample_multinomial(prob); // new topic

			// increments the counts by one

			beta_counts(new_topic, wid)++;
			theta_counts(new_topic, did)++;
			topic_counts(new_topic)++;
			z(idx) = new_topic;

		} // end of the word topic sampling loop

		if ((iter >= burn_in) && (iter % spacing == 0)){ // handles the burn in period
			Z.col(count) = z;
			if (store_dirichlet == 1){
				thetas.slice(count) = theta_counts;
				betas.slice(count) = beta_counts;
			}

			log_marginal(count) = calc_log_marginal_posterior(num_topics, num_docs, vocab_size, theta_counts, beta_counts);
			if (iter % 1000 == 0)
				cout << " lmp: " << log_marginal(count);

			count++;
		}

		if (iter % 1000 == 0)
			cout << endl;

	} // The end of the Gibbs loop

	if (store_dirichlet == 1){
		return List::create(
				Named("thetas") = wrap(thetas),
				Named("betas") = wrap(betas),
				Named("Z") = wrap(Z),
				Named("lmp") = wrap(log_marginal));
	}
	else {
		return List::create(
				Named("Z") = wrap(Z),
				Named("lmp") = wrap(log_marginal));
	}

}



/**
 *  The LDA z-theta fixed beta gibbs sampler:
 * 		This function samples z, \theta, with a fixed \beta
 * 		This function assumes that the document word
 * 		instances are in Gibbs sampling format (mainly,
 * 		used for testing).
 *
 * 	Arguments:
 * 		doc_lengths_ - contains the number of words in each document as a vector
 * 		word_ids_ - the document word instances
 * 		topic_assignments_ - initial topic assignments for each word in the corpus
 * 		alpha_v_ - the hyperparameter vector for the document topic Dirichlet
 * 		beta_ - the fixed beta matrix
 * 		max_iter_ - max number of Gibbs iterations
 * 		burn_in_ - burn in period
 * 		spacing_ - spacing between samples that are stored
 *
 * 	Returns:
 * 		thetas - the sampled thetas after burn in period
 * 		betas - the sampled betas after burn in period
 * 		Z - the sampled word topic assignments after burn in period
 *
 */
RcppExport SEXP lda_z_theta_fixed_beta(SEXP doc_lengths_, SEXP word_ids_,
		SEXP topic_assignments_, SEXP alpha_v_, SEXP beta_,
		SEXP max_iter_, SEXP burn_in_, SEXP spacing_, SEXP store_dirichlet_) {

	uvec doc_N = as<uvec>(doc_lengths_);
	uvec word_ids = as<uvec>(word_ids_);
	uvec z = as<uvec>(topic_assignments_);
	vec alpha_v = as<vec>(alpha_v_);
	mat beta_m = as<mat>(beta_);
	int max_iter = as<int>(max_iter_);
	int burn_in = as<int>(burn_in_);
	int sample_spacing = as<int>(spacing_);
	int store_dirichlet = as<int>(store_dirichlet_);

	int num_topics = beta_m.n_rows;
	int vocab_size = beta_m.n_cols;
	int num_docs = doc_N.n_elem;
	int num_word_instances = word_ids.n_elem;
	int valid_samples = ceil((max_iter - burn_in) / (double) sample_spacing);

	mat prior_theta_samples = zeros<mat>(num_topics, num_docs);
	cube thetas;
	if (store_dirichlet == 1)
		thetas = cube(num_topics, num_docs, valid_samples);
	umat Z = zeros<umat>(num_word_instances, valid_samples);

	vector < vector < unsigned int > > document_word_indices;
	unsigned int d, i, iter, instances = 0, count = 0;

	// Calculates the document word indices

	for (d = 0; d < num_docs; d++){
		vector < unsigned int > word_idx;
		for (i = 0; i < doc_N(d); i++){
			word_idx.push_back(instances);
			instances++;
		}
		document_word_indices.push_back(word_idx);
	}

	uvec prior_z = z;

	// The Gibbs sampling loop

	for (iter = 0; iter < max_iter; iter++){

		if (iter % 1000 == 0)
			cout << "gibbs iter# " << iter + 1;

		for (d = 0; d < num_docs; d++){

			vector < unsigned int > word_idx = document_word_indices[d];

			vec partition_counts = alpha_v;
			for (i = 0; i < doc_N(d); i++)
				partition_counts(z(word_idx[i])) += 1;
			vec theta_d = sample_dirichlet(num_topics, partition_counts);
			prior_theta_samples.col(d) = theta_d; // partition_counts; //

			for (i = 0; i < doc_N(d); i++)
				z(word_idx[i]) = sample_multinomial(theta_d % beta_m.col(word_ids(word_idx[i])));

		}

		if ((iter >= burn_in) && (iter % sample_spacing == 0)){
			Z.col(count) = prior_z;
			if (store_dirichlet == 1)
				thetas.slice(count) = prior_theta_samples;
			count++;
		}

		prior_z = z;
		if (iter % 1000 == 0)
			cout << endl;
	}

	if (store_dirichlet == 1) {
		return List::create(
				Named("thetas") = wrap(thetas),
				Named("Z") = wrap(Z));
	}
	else {
		return List::create(Named("Z") = wrap(Z));
	}
}

