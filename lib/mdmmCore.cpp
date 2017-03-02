// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <stdexcept>
#define log2pi 1.8378770664093453


// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Function Definitions
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

// n_k
std::vector<double> z_sum_i(NumericMatrix& z, int& nrow, int& K_clust)	{
	std::vector<double> out (K_clust,0);
	for (int k = 0; k < K_clust; k++) {
		for (int i = 0; i < nrow; i++) {
			out[k] += z(i,k);
		}
	}
	return out;
}

// gaussian mean estimation function
NumericMatrix muEM( NumericMatrix& x, NumericMatrix& z, 
					std::vector<double>& sum_zk,
				    int& nrow, int& ncol, int& K_clust)	{
	NumericMatrix out( K_clust, ncol);
	for(int k = 0; k < K_clust; k++)	{
		for(int j = 0; j < ncol; j++)	{
			for(int i = 0; i < nrow; i++)	{
				out(k,j) += x(i,j) * z(i,k);
			}
			out(k,j) = out(k,j) / sum_zk[k];
		}
	}
	return out;
}

// gaussian covariance estimation function
List covEM( NumericMatrix& x, NumericMatrix& z, NumericMatrix& mu, 
			std::vector<double>& sum_zk, int& nrow, int& ncol, int& K_clust)	{
	List out(K_clust);
	// iterating through each cluster
	for( int k = 0; k < K_clust; k++ ) {
		NumericMatrix out_k(ncol,ncol);
		// iterating through each data point
		for( int h = 0; h < nrow; h++ ) {
			// Covariance Matrix for each data point
			for( int i = 0; i < ncol; i++ ) {
				for( int j = i; j < ncol; j++ ) {
					out_k(i,j) += z(h,k)*(x(h,i)-mu(k,i))*(x(h,j)-mu(k,j));
					out_k(j,i) = out_k(i,j);
				}
			}
		}
		for( int i = 0; i < ncol; i++ ) {
			for( int j = i; j < ncol; j++ ) {
				out_k(i,j) = out_k(i,j) / sum_zk[k];
				out_k(j,i) = out_k(i,j);
			}
		}
		out(k) = out_k;
	}
	return out;
}

// gaussian diagonal covariance estimation function
List covEM_diag(NumericMatrix& x, NumericMatrix& z, NumericMatrix& mu, 
				std::vector<double>& sum_zk, int& nrow, int& ncol, int& K_clust)	{
	List out(K_clust);
	// iterating through each cluster
	for( int k = 0; k < K_clust; k++ ) {
		NumericMatrix out_k(ncol,ncol);
		// iterating through each data point
		for( int h = 0; h < nrow; h++ ) {
			// Covariance Matrix for each data point
			for( int i = 0; i < ncol; i++ ) {
					out_k(i,i) += z(h,k)*(x(h,i)-mu(k,i))*(x(h,i)-mu(k,i));
			}
		}
		for( int i = 0; i < ncol; i++ ) {
				out_k(i,i) = out_k(i,i) / sum_zk[k];
		}
		out(k) = out_k;
	}
	return out;
}

// categorical/multinomial probability estimation
List catEM( List& x, NumericMatrix& z, 
			std::vector<double>& sum_zk, int& nvar, int& nrow,
			IntegerVector& ncol, int& K_clust)	{
	List out(nvar);
	for ( int l = 0; l < nvar; l++ )	{
		NumericMatrix x_l = x[l], out_l( K_clust, ncol[l]);
		for ( int k = 0; k < K_clust; k++ )	{
			for ( int i = 0; i < nrow; i++ ) {
				for ( int j = 0; j < ncol[l]; j++ ) {
					out_l(k,j) += z(i,k)*x_l(i,j);
				}
			}
			for ( int j = 0; j < ncol[l]; j++ ) {
				out_l(k,j) = out_l(k,j) / sum_zk[k];
			}
		}
		out(l) = out_l;
	}
	return out;
}

// multinomial probability estimation
List multEM(List& x, NumericMatrix& sum_x, 
			NumericMatrix& z, int& nvar, int& nrow,
			IntegerVector& ncol, int& K_clust)	{
	List out(nvar);
	double out_den = 0;
	for ( int l = 0; l < nvar; l++ )	{
		NumericMatrix x_l = x[l], out_l( K_clust, ncol[l]);
		for ( int k = 0; k < K_clust; k++ )	{
			for ( int i = 0; i < nrow; i++ ) {
				for ( int j = 0; j < ncol[l]; j++ ) {
					out_l(k,j) += z(i,k)*x_l(i,j);
				}
				out_den += z(i,k)*sum_x(i,l);
			}
			for ( int j = 0; j < ncol[l]; j++ ) {
				out_l(k,j) = out_l(k,j) / out_den;
			}
			out_den = 0;
		}
		out(l) = out_l;
	}
	return out;
}

// cluster probability estimation
NumericVector clust_prob(std::vector<double>& sum_zk, int& K_clust)	{
	NumericVector out(K_clust);
	double sum_z = 0;
	for (int k = 0; k < K_clust; k++ ) {
		sum_z += sum_zk[k];
	}
	for (int k = 0; k < K_clust; k++ ) {
		out[k] = sum_zk[k] / sum_z;
	}
	return out;
}

// gaussian log likelihood
NumericMatrix lmvgauss(NumericMatrix& x, NumericMatrix& mu, List& sig,
					     int& nrow, int& ncol, int& K_clust)	{
	
	std::ostream nullstream(0);
	arma::set_stream_err2(nullstream);
	NumericMatrix out(nrow, K_clust);
	arma::mat sig_inv = arma::zeros(ncol,ncol);
	bool flag = false;
	for (int k = 0; k < K_clust; k++ ) {
		double sign, lsig_det;
		arma::mat sig_k = sig[k];
		flag = arma::inv_sympd(sig_inv,sig_k);
		if(!flag) {
			Rcout << "\nNote: Covariance matrix for cluster " << k + 1 << " could not be inverted. Please" << std::endl;
			Rcout << "\t[1]repeat the function call or \n\t[2]reduce the number of clusters if the problem persists." << std::endl;
			throw std::runtime_error("matrix inversion failed");
		}
		log_det(lsig_det, sign, sig_k);
		for (int h = 0; h < nrow; h++ ) {
			for (int i = 0; i < ncol; i++ ) {
				for (int j = 0; j < ncol; j++ ) {
					out(h,k) += (x(h,i)-mu(k,i)) * sig_inv(i,j) * (x(h,j)-mu(k,j));
				}
			}
			out(h,k) = -0.5 * (ncol*log2pi + lsig_det + out(h,k));
		}
	}
	return out;
}

// categorical log likelihood
NumericMatrix lcat(List& x, List& p, int& nvar, 
				     int& nrow, IntegerVector& ncol, int& K_clust)	{
	NumericMatrix out(nrow,K_clust);
	for ( int l = 0; l < nvar; l++ )	{
		NumericMatrix x_l = x[l], p_l = p[l];
		for ( int k = 0; k < K_clust; k++ )	{
			for ( int i = 0; i < nrow; i++ ) {
				for ( int j = 0; j < ncol[l]; j++ ) {
					out(i,k) += x_l(i,j)*std::log(p_l(k,j));
				}
			}
		}
	}
	return out;
}

// multinomial sum_j x_ij for each i
NumericMatrix mult_xsum(List& x, int& nvar, int& nrow,
						  IntegerVector& ncol) {
	NumericMatrix out(nrow,nvar);
	for ( int l = 0; l < nvar; l++ ) {
		NumericMatrix x_l = x[l];
		for ( int i = 0; i < nrow; i++ ) {
			for ( int j = 0; j < ncol[l]; j++ ) {
				out(i,l) += x_l(i,j);
			}
		}
	}
	return out;
}

// log multinomial constant
NumericMatrix lmultconst(List& x, int& nvar, int& nrow,
						  IntegerVector& ncol) {
	NumericMatrix out(nrow,nvar);
	double sum_x_l =0;
	for ( int l = 0; l < nvar; l++ ) {
		NumericMatrix x_l = x[l];
		for ( int i = 0; i < nrow; i++ ) {
			for ( int j = 0; j < ncol[l]; j++ ) {
				out(i,l) -= lgamma(x_l(i,j)+1);
				sum_x_l += x_l(i,j);
			}
			out(i,l) += lgamma(sum_x_l+1);
			sum_x_l = 0;
		}
	}
	return out;
}

// multinomial log likelihood
NumericMatrix lmult(List& x, NumericMatrix& lmultconst_x, List& p, int& nvar, 
			           int& nrow, IntegerVector& ncol, int& K_clust)	{
	NumericMatrix out(nrow,K_clust);
	for ( int l = 0; l < nvar; l++ )	{
		NumericMatrix x_l = x[l], p_l = p[l];
		for ( int k = 0; k < K_clust; k++ )	{
			for ( int i = 0; i < nrow; i++ ) {
				for ( int j = 0; j < ncol[l]; j++ ) {
					out(i,k) += x_l(i,j)*std::log(p_l(k,j));
				}
				out(i,k) += lmultconst_x(i,l);
			}
		}
	}
	return out;
}

// sum of log likelihoods
NumericMatrix loglik( NumericMatrix& lgauss_k, NumericMatrix& lcat_k, 
					  NumericMatrix& lmult_k, NumericVector& Pi,
					  int& nrow, int& K_clust)	{
	NumericMatrix out(nrow, K_clust);
	double lpi[K_clust];
	for ( int k = 0; k < K_clust; k++ )	{
		lpi[k] = std::log(Pi[k]);
	}
	for ( int k = 0; k < K_clust; k++ )	{
		for ( int i = 0; i < nrow; i++ ) {
			out(i,k) += lgauss_k(i,k) + lcat_k(i,k) + lmult_k(i,k) + lpi[k];			
		}
	}
	return out;
}

// data point cluster probability
NumericMatrix data_prob(NumericMatrix& loglik_k, int& nrow, int& K_clust)	{
	NumericMatrix out(nrow, K_clust);
	NumericVector max_i(nrow);
	max_i = as<NumericVector>(wrap(max(as<arma::mat>(loglik_k),1)));
	long double logsumexp = 0;
	for ( int i = 0; i < nrow; i++ ) {
		for ( int k = 0; k < K_clust; k++ )	{
			logsumexp += std::exp(loglik_k(i,k)-max_i[i]);
		}
		logsumexp = max_i[i] + std::log(logsumexp);
		for ( int k = 0; k < K_clust; k++ )	{
			out(i,k) = std::exp(loglik_k(i,k) - logsumexp) ;			
		}
		logsumexp = 0;
	}
	return out;
}

//random dirichlet distribution
std::vector<double> rDirichlet(int& nrow, std::vector<double>& alp)	{
	
	int K = alp.size();
	std::vector<double> out (nrow*K,1);
	double sum_k = 0;
	
	for (int i = 0; i < nrow; i++)	{
		sum_k = 0;
		for (int j = 0; j < K; j++)	{
			out[j+(i*K)] = R::rgamma(alp[j],1);
			sum_k += out[j+(i*K)];
		}
		for (int j = 0; j < K; j++)	{
			out[j+(i*K)] = out[j+(i*K)]/sum_k;
		}
	}
	return out;
}

// transform rDirichlet to R matrix form
NumericMatrix cDirichlet2(int& nrow, NumericVector& alp) {
	
	int K = alp.size();
	std::vector<double> alpvec = as< std::vector<double> >(alp);
	std::vector<double> in = rDirichlet(nrow, alpvec);
	NumericMatrix out(nrow,K);
		
	for(int i = 0; i < nrow; i++)  {
        for(int j = 0; j < K; j++)  {
             out(i,j) = in[j + i*K];
        }
	}
	return out;
}
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo





// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Main loop
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//[[Rcpp::export]]
List mdmmCore(NumericMatrix& Xg, int& Xg_ncol,
			  List& Xc, IntegerVector& Xc_ncol, int& Xc_nvar,
			  List& Xm, IntegerVector& Xm_ncol, int& Xm_nvar,
			  int& nrow, int& K_clust, double& tol, bool& est_cov,
			  int& attempts, std::string& criterion)	{

	// Defining constants and repetitively utilized calculations
	NumericMatrix logmultconst(nrow,Xm_nvar), multsum(nrow,Xm_nvar);
	if (Xm_nvar>0)
		logmultconst = lmultconst(Xm, Xm_nvar, nrow, Xm_ncol);
		multsum = mult_xsum(Xm, Xm_nvar, nrow, Xm_ncol);
	std::vector<double> sum_zk (K_clust,0);
	
	// distribution variables
	NumericMatrix mu_k( K_clust, Xg_ncol);
	List sig_k(K_clust), pc_k(Xc_nvar), pm_k(Xm_nvar);
	NumericVector pi_k(K_clust);
	// initialize likelihood variables
	NumericMatrix lmvgauss_k(nrow,K_clust), lcat_k(nrow,K_clust);
	NumericMatrix lmult_k(nrow,K_clust), loglik_k(nrow,K_clust);
	
	// initializing loop
	List out;
	int iter_attempts = 0;
	
	// Parameter count for Goodness of fit statistics
	int N_param = 0;
	N_param += Xg_ncol;
	if(est_cov == true)	{
		N_param += Xg_ncol*(Xg_ncol-1)/2;
	} else {
		N_param += Xg_ncol;
	}
	if (Xc_nvar>0)
		for(int i = 0; i < Xc_nvar; i++)
			N_param += Xc_ncol[i];
	if (Xm_nvar>0)
		for(int i = 0; i < Xm_nvar; i++)
			N_param += Xm_ncol[i];
	N_param *= K_clust;
	// Goodness of fit statistics
	double AIC[2] = {0}, BIC[2] = {0};
	AIC[1] = std::numeric_limits<double>::infinity();
	BIC[1] = std::numeric_limits<double>::infinity();
	NumericVector GOF(2);
	
	while( iter_attempts < attempts )	{
		
		Rcout << "------------ Attempt ID: " << iter_attempts + 1 << " ------------" << std::endl;
		
		// initialize latent variables
		NumericVector alp(K_clust, 1.0);
		NumericMatrix Z_t(nrow,K_clust), Z_tp1 = cDirichlet2(nrow, alp);
		
		// initializing EM iteration
		double eps = tol + 1, eps_max = 0;
		std::vector<double> loglik_iter; loglik_iter.reserve(70);
		unsigned int iter = 0;
		bool err_break = false;
		
		while( eps > tol )	{
			
			loglik_iter.push_back(0);
			Z_t = Z_tp1;
			sum_zk = z_sum_i(Z_t, nrow, K_clust);
			
			// ----- [1] Maximization Step -------------------------------------
			// gaussian distribution parameter estimation
			if (Xg_ncol>0)	{
				mu_k = muEM(Xg, Z_t, sum_zk, nrow, Xg_ncol, K_clust);
				if ( est_cov == true ) {
					sig_k = covEM(Xg, Z_t, mu_k, sum_zk, nrow, Xg_ncol, K_clust);
				} else {
					sig_k = covEM_diag(Xg, Z_t, mu_k, sum_zk, nrow, Xg_ncol, K_clust);
				}
			}
			
			// categorical distribution parameter estimation
			if (Xc_nvar>0)
				pc_k = catEM(Xc, Z_t, sum_zk, Xc_nvar, nrow, Xc_ncol, K_clust);
			if (Xm_nvar>0)
				pm_k = multEM(Xm, multsum, Z_t, Xm_nvar, nrow, Xm_ncol, K_clust);
			
			// cluster probability
			pi_k = clust_prob(sum_zk, K_clust);
			
			// ----- [2] Expectation Step --------------------------------------
			if (Xg_ncol>0)
				try {
					lmvgauss_k = lmvgauss(Xg, mu_k, sig_k, nrow, Xg_ncol, K_clust);
				} catch (std::runtime_error& err) {
					Rcout << "Error: " << err.what() << "\n" << std::endl;
					err_break = true;
					break;
				}
			
			if (Xc_nvar>0)
				lcat_k = lcat(Xc, pc_k, Xc_nvar, nrow, Xc_ncol, K_clust);
			
			if (Xm_nvar>0)
				lmult_k = lmult(Xm, logmultconst, pm_k, Xm_nvar, nrow, Xm_ncol, K_clust);
			
			loglik_k = loglik( lmvgauss_k, lcat_k, lmult_k, pi_k, nrow, K_clust);
			
			Z_tp1 = data_prob(loglik_k, nrow, K_clust);
			
			// ----- [3] Full Likelihood and Convergence -----------------------
			for( int i = 0; i < nrow; i++ )	{
				for( int k = 0; k < K_clust; k++ )	{
					loglik_iter[iter] += Z_t(i,k)*loglik_k(i,k);
					if ( (Z_tp1(i,k) - Z_t(i,k)) > eps_max )
						eps_max = Z_tp1(i,k) - Z_t(i,k);
				}
			}
			eps = eps_max;
			eps_max = 0;
			
			iter++;
			if (iter % 1 == 0) {
				Rcout << "Current Iteration: " << iter << std::endl;
				Rcout << "      Convergence: " << eps << std::endl;
				Rcout << "    loglikelihood: " << loglik_iter[iter-1] << std::endl;
			}
		}
		
		//conduct output assignment only if no errors occurred.
		if (err_break == false)	{
			AIC[1] = 2*(N_param - loglik_iter[iter-1]);
			BIC[1] = -2*loglik_iter[iter-1]+N_param*std::log(nrow);
			if (criterion == "BIC" && iter_attempts > 0) {
				if( BIC[1] < BIC[0] ) {
					GOF[0] = AIC[1]; GOF[1] = BIC[1];
					AIC[0] = AIC[1]; BIC[0] = BIC[1];
					out = List::create( _["GOF"] = GOF,
										_["data_prob"] = Z_tp1, 
										_["loglikelihood"] = loglik_iter, 
										_["cluster_prob"] = pi_k, 
										_["gaussian_mean"] = mu_k, 
										_["gaussian_var"] = sig_k,
										_["categorical_prob"] = pc_k, 
										_["multinomial_prob"] = pm_k);
				}
			} else if (criterion == "AIC" && iter_attempts > 0){
				if( AIC[1] < AIC[0] ) {
					GOF[0] = AIC[1]; GOF[1] = BIC[1];
					AIC[0] = AIC[1]; BIC[0] = BIC[1];
					out = List::create( _["GOF"] = GOF,
										_["data_prob"] = Z_tp1, 
										_["loglikelihood"] = loglik_iter, 
										_["cluster_prob"] = pi_k, 
										_["gaussian_mean"] = mu_k, 
										_["gaussian_var"] = sig_k,
										_["categorical_prob"] = pc_k, 
										_["multinomial_prob"] = pm_k);
				}
			} else {
				GOF[0] = AIC[1]; GOF[1] = BIC[1];
				out = List::create( _["GOF"] = GOF,
									_["data_prob"] = Z_tp1, 
									_["loglikelihood"] = loglik_iter, 
									_["cluster_prob"] = pi_k, 
									_["gaussian_mean"] = mu_k, 
									_["gaussian_var"] = sig_k,
									_["categorical_prob"] = pc_k, 
									_["multinomial_prob"] = pm_k);
			}
			
			if (loglik_iter[iter-1] != loglik_iter[iter-1]) {
				Rcout << "\nNote: Loglikelihood is returning NaN. Please" << std::endl;
				Rcout << "\t[1]repeat the function call, \n\t[2]drop variables which may be causing the NaNs or \n\t[3]reduce the number of clusters if the problem persists." << std::endl;
			}
		}
		
		Rcout << "Attempt " << iter_attempts + 1 << " completed in " << iter + 1 << " iterations" << std::endl;
		Rcout << "---------------------------------------" << std::endl;
		iter_attempts++;
	}
	return out;
}
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


