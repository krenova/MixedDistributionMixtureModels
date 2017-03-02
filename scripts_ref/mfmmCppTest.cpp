#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// using namespace arma;
#define log2pi        1.8378770664093453

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Guassian Mean
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//[[Rcpp::export]]
NumericMatrix meanEM(NumericMatrix x, NumericMatrix z, 
				   const int nrow, const int ncol, const int K_clust)	{
	NumericMatrix out( K_clust, ncol);
	for(int k = 0; k < K_clust; k++)	{
		for(int j = 0; j < ncol; j++)	{
			for(int i = 0; i < nrow; i++)	{
				out(k,j) += x(i,j) * z(i,k);
			}
		}
		out.row(k) = out.row(k) / sum(z.column(k));
	}
	return out;
}



// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Elliptical Covariance
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

//[[Rcpp::export]]
List covEM1(NumericMatrix x, NumericMatrix z, NumericMatrix mu,
			 const int nrow, const int ncol, const int K_clust)	{
	List out(K_clust);
	// Precomputing the cluster counts total
	NumericVector sum_z(K_clust);
	for ( int k = 0; k < K_clust; k++ )	{
		for ( int h = 0; h < nrow; h++ )	{
			sum_z[k] += z(h,k);
		}
	}
	// iterating through each cluster
	for( int k = 0; k < K_clust; k++ ) {
		NumericMatrix out_k(ncol,ncol);
		// iterating through each data point
		for( int h = 0; h < nrow; h++ ) {
			// Covariance Matrix for each data point
			for( int i = 0; i < ncol; i++ ) {
				for( int j = i; j < ncol; j++ ) {
					out_k(i,j) += z(h,k)*(x(h,i)-mu(k,i))*(x(h,j)-mu(k,j))/sum_z[k];
					out_k(j,i) = out_k(i,j);
				}
			}
		}
		out(k) = out_k;
	}
	return out;
}

//[[Rcpp::export]]
List covEM2(NumericMatrix x, NumericMatrix z, NumericMatrix mu,
			 const int nrow, const int ncol, const int K_clust)	{
	List out(K_clust);
	// Precomputing the cluster counts total
	NumericVector sum_z(K_clust);
	for ( int k = 0; k < K_clust; k++ )	{
		for ( int h = 0; h < nrow; h++ )	{
			sum_z[k] += z(h,k);
		}
	}
	// iterating through each cluster
	for( int k = 0; k < K_clust; k++ ) {
		NumericMatrix out_k(ncol,ncol);
		// iterating through each data point
		for( int h = 0; h < nrow; h++ ) {
			// Covariance Matrix for each data point
			for( int i = 0; i < ncol; i++ ) {
				for( int j = i; j < ncol; j++ ) {
					out_k(i,j) += z(h,k)*(x(h,i)-mu(k,i))*(x(h,j)-mu(k,j))/sum_z[k];
				}
			}
		}
		out(k) = out_k;
	}
	return out;
}

//[[Rcpp::export]]
NumericVector covEM3(NumericMatrix x, NumericMatrix z, NumericMatrix mu,
					 const int nrow, const int ncol, const int K_clust)	{
	NumericVector out(ncol*ncol*K_clust);
	// Precomputing the cluster counts total
	NumericVector sum_z(K_clust);
	for ( int k = 0; k < K_clust; k++ )	{
		for ( int h = 0; h < nrow; h++ )	{
			sum_z[k] += z(h,k);
		}
	}
	// iterating through each cluster
	for( int k = 0; k < K_clust; k++ ) {
		// iterating through each data point
		for( int h = 0; h < nrow; h++ ) {
			// Covariance Matrix for each data point
			for(int i = 0; i < ncol; i++)	{
				for(int j = i; j < ncol; j++)	{
					out[i+(ncol*j)+(ncol*ncol*k)] +=z(h,k)*(x(h,i)-mu(k,i))*(x(h,j)-mu(k,j))/sum_z[k];
					out[j+(ncol*i)+(ncol*ncol*k)] = out[i+(ncol*j)+(ncol*ncol*k)];
				}
			}
		}
	}
	out.attr("dim") = Dimension(ncol,ncol,K_clust);
	return out;
}


//[[Rcpp::export]]
NumericVector covEM4(NumericMatrix x, NumericMatrix z, NumericMatrix mu,
					 const int nrow, const int ncol, const int K_clust)	{
	NumericVector out(ncol*ncol*K_clust);
	// Precomputing the cluster counts total
	NumericVector sum_z(K_clust);
	for ( int k = 0; k < K_clust; k++ )	{
		for ( int h = 0; h < nrow; h++ )	{
			sum_z[k] += z(h,k);
		}
	}
	// iterating through each cluster
	for( int k = 0; k < K_clust; k++ ) {
		// iterating through each data point
		for( int h = 0; h < nrow; h++ ) {
			// Covariance Matrix for each data point
			for(int i = 0; i < ncol; i++)	{
				for(int j = i; j < ncol; j++)	{
					out[i+(ncol*j)+(ncol*ncol*k)] +=z(h,k)*(x(h,i)-mu(k,i))*(x(h,j)-mu(k,j))/sum_z[k];
					if( i != j) out[j+(ncol*i)+(ncol*ncol*k)] = out[i+(ncol*j)+(ncol*ncol*k)];
				}
			}
		}
	}
	out.attr("dim") = Dimension(ncol,ncol,K_clust);
	return out;
}
//[[Rcpp::export]]
NumericVector covEM5(NumericMatrix x, NumericMatrix z, NumericMatrix mu,
					 const int nrow, const int ncol, const int K_clust)	{
	NumericVector out(ncol*ncol*K_clust);
	// Precomputing the cluster counts total
	NumericVector sum_z(K_clust);
	for ( int k = 0; k < K_clust; k++ )	{
		for ( int h = 0; h < nrow; h++ )	{
			sum_z[k] += z(h,k);
		}
	}
	// iterating through each cluster
	for( int k = 0; k < K_clust; k++ ) {
		// iterating through each data point
		for( int h = 0; h < nrow; h++ ) {
			// Covariance Matrix for each data point
			for(int i = 0; i < ncol; i++)	{
				for(int j = i; j < ncol; j++)	{
					out[i+(ncol*j)+(ncol*ncol*k)] +=z(h,k)*(x(h,i)-mu(k,i))*(x(h,j)-mu(k,j))/sum_z[k];
					out[j+(ncol*i)+(ncol*ncol*k)] = out[i+(ncol*j)+(ncol*ncol*k)];
				}
			}
		}
	}
	return out;
}
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Spherical Covariance
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

//[[Rcpp::export]]
List covEM_diag(NumericMatrix x, NumericMatrix z, NumericMatrix mu,
				const int nrow, const int ncol, const int K_clust)	{
	List out(K_clust);
	// Precomputing the cluster counts total
	NumericVector sum_z(K_clust);
	for ( int k = 0; k < K_clust; k++ )	{
		for ( int h = 0; h < nrow; h++ )	{
			sum_z[k] += z(h,k);
		}
	}
	// iterating through each cluster
	for( int k = 0; k < K_clust; k++ ) {
		NumericMatrix out_k(ncol,ncol);
		// iterating through each data point
		for( int h = 0; h < nrow; h++ ) {
			// Covariance Matrix for each data point
			for( int i = 0; i < ncol; i++ ) {
					out_k(i,i) += z(h,k)*(x(h,i)-mu(k,i))*(x(h,i)-mu(k,i))/sum_z[k];
			}
		}
		out(k) = out_k;
	}
	return out;
}

//[[Rcpp::export]]
NumericVector covEM_diag2(NumericMatrix x, NumericMatrix z, NumericMatrix mu,
				const int nrow, const int ncol, const int K_clust)	{
	NumericVector out(ncol*ncol*K_clust);
	// Precomputing the cluster counts total
	NumericVector sum_z(K_clust);
	for ( int k = 0; k < K_clust; k++ )	{
		for ( int h = 0; h < nrow; h++ )	{
			sum_z[k] += z(h,k);
		}
	}
	// iterating through each cluster
	for( int k = 0; k < K_clust; k++ ) {
		// iterating through each data point
		for( int h = 0; h < nrow; h++ ) {
			// Covariance Matrix for each data point
			for( int i = 0; i < ncol; i++ ) {
					out[i+(ncol*i)+(ncol*ncol*k)] +=z(h,k)*(x(h,i)-mu(k,i))*(x(h,i)-mu(k,i))/sum_z[k];
			}
		}
	}
	out.attr("dim") = Dimension(ncol,ncol,K_clust);
	return out;
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Categorical/Multinomial Probability
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//[[Rcpp::export]]
List multEM(List x, NumericMatrix z, 
			const int nvar, const int nrow,
			IntegerVector ncol, const int K_clust)	{
	List out(nvar);
	double sum_out_l = 0;
	for ( int l = 0; l < nvar; l++ )	{
		NumericMatrix x_l = x[l], out_l( K_clust, ncol[l]);
		for ( int k = 0; k < K_clust; k++ )	{
			for ( int j = 0; j < ncol[l]; j++ ) {
				for ( int i = 0; i < nrow; i++ ) {
					out_l(k,j) += z(i,k)*x_l(i,j);
				}
				sum_out_l += out_l(k,j);
			}
			for ( int j = 0; j < ncol[l]; j++ ) {
				out_l(k,j) = out_l(k,j) / sum_out_l;
			}
			sum_out_l = 0;
		}
		out(l) = out_l;
	}
	return out;
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Cluster Probability
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//[[Rcpp::export]]
NumericVector clust_prob(NumericMatrix z, const int nrow, const int K_clust)	{
	NumericVector out(K_clust);
	double sum_out = 0;
	for (int k = 0; k < K_clust; k++ ) {
		for (int i = 0; i < nrow; i++ ) {
			out[k] += z(i,k);
		}
		sum_out += out[k];
	}
	for (int k = 0; k < K_clust; k++ ) {
		out[k] = out[k]/sum_out;
	}
	return out;
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// log gaussian likelihood
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//[[Rcpp::export]]
NumericMatrix lmvgauss_k(NumericMatrix x, NumericMatrix mu, List sig,
					     const int nrow, const int ncol, const int K_clust)	{
	NumericMatrix out(nrow, K_clust);
	arma::mat sig_inv = arma::zeros(ncol,ncol);
	for (int k = 0; k < K_clust; k++ ) {
		double sign, sig_det;
		arma::mat sig_k = sig[k];
		sig_inv = arma::inv_sympd(sig_k);
		log_det(sig_det, sign, sig_k);
		for (int h = 0; h < nrow; h++ ) {
			for (int i = 0; i < ncol; i++ ) {
				for (int j = 0; j < ncol; j++ ) {
					out(h,k) += (x(h,i)-mu(k,i)) * sig_inv(i,j) * (x(h,j)-mu(k,j));
				}
			}
			out(h,k) = -0.5 * (ncol*log2pi + sig_det + out(h,k));
		}
	}
	return out;
}

//[[Rcpp::export]]
NumericMatrix lcat_k(List x, List p, const int nvar, 
				     const int nrow, IntegerVector ncol, const int K_clust)	{
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

//[[Rcpp::export]]
NumericMatrix lmultconst(List& x, const int& nvar, const int& nrow, 
						  IntegerVector& ncol) {
	NumericMatrix out(nrow,nvar);
	double sum_x_l = 0;
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

//[[Rcpp::export]]
NumericMatrix lmult_k2(List& x, NumericMatrix& lmultconst_x, List& p, const int& nvar, 
			           const int& nrow, IntegerVector& ncol, const int& K_clust)	{
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

//[[Rcpp::export]]
NumericMatrix lmult_k(List x, List p, const int nvar, 
			           const int nrow, IntegerVector ncol, const int K_clust)	{
	NumericMatrix out(nrow,K_clust);
	double sum_x_l = 0;
	for ( int l = 0; l < nvar; l++ )	{
		NumericMatrix x_l = x[l], p_l = p[l];
		for ( int k = 0; k < K_clust; k++ )	{
			for ( int i = 0; i < nrow; i++ ) {
				for ( int j = 0; j < ncol[l]; j++ ) {
					out(i,k) += -lgamma(x_l(i,j)+1) + x_l(i,j)*std::log(p_l(k,j));
					sum_x_l += x_l(i,j);
				}
				out(i,k) += lgamma(sum_x_l+1);
				sum_x_l = 0;
			}
		}
	}
	return out;
}

//[[Rcpp::export]]
NumericMatrix loglik_k( NumericMatrix lgauss_k, NumericMatrix lcat_k, 
						NumericMatrix lmult_k, NumericVector Pi,
						const int nrow, const int K_clust)	{
	NumericMatrix out(nrow, K_clust);
	for ( int k = 0; k < K_clust; k++ )	{
		for ( int i = 0; i < nrow; i++ ) {
			out(i,k) += lgauss_k(i,k) + lcat_k(i,k) + lmult_k(i,k) + std::log(Pi[k]);			
		}
	}
	return out;
}


//[[Rcpp::export]]
NumericMatrix data_prob(NumericMatrix& loglik_k, const int& nrow, const int& K_clust)	{
	NumericMatrix out(nrow, K_clust);
	NumericVector max_i(nrow);
	max_i = as<NumericVector>(wrap(max(as<arma::mat>(loglik_k),1)));
	long double logsumexp = 0;
	for ( int i = 0; i < nrow; i++ ) {
		logsumexp = 0;
		for ( int k = 0; k < K_clust; k++ )	{
			logsumexp += exp(loglik_k(i,k)-max_i[i]);
		}
		logsumexp = max_i[i] + log(logsumexp);
		for ( int k = 0; k < K_clust; k++ )	{
			out(i,k) = exp(loglik_k(i,k) - logsumexp) ;			
		}
	}
	return out;
}








/*
//[[Rcpp::export]]
List clmult(NumericMatrix Xg, NumericMatrix mu, List sig,
			List Xc_l, List Pc_l, IntegerVector ncolc, const int nvarc,
			List Xm_l, List Pm_l,  IntegerVector ncolm, const int nvarm,
			const int nrow, const int K_clust)	{
	List out(nvar);
	double sum_x_l = 0;
	for ( int k = 0; k < K_clust; k++ )	{
		N
		
		
		
	}
	return out;
}
*/







/*

List multEM(List x, NumericMatrix z, 
			const int nvar, const int nrow,
			IntegerVector ncol, const int K_clust)	{
	List out(nvar);
	double sum_out_l = 0;
	for ( int l = 0; l < nvar; l++ )	{
		NumericMatrix x_l = x[l], out_l( K_clust, ncol[l]);
		for ( int k = 0; k < K_clust; k++ )	{
			for ( int j = 0; j < ncol[l]; j++ ) {
				for ( int i = 0; i < nrow; i++ ) {
					out_l(k,j) += z(i,k)*x_l(i,j);
				}
				sum_out_l += out_l(k,j);
			}
			for ( int j = 0; j < ncol[l]; j++ ) {
				out_l(k,j) = out_l(k,j) / sum_out_l;
			}
			sum_out_l = 0;
		}
		out(l) = out_l;
	}
	return out;
}


/*


lapply(lapply(Pc_kl,function(x) lapply(x,log)), 
       function(p) mapply( function(p_2,x) x %*% p_2, p, Xc_l)	
		)
						 
						 

-0.5*( d*log(2*pi) + determinant(sig)$modulus + 
                  t(x-mu) %*% (chol2inv(chol(sig)) %*% (x-mu)) ))
				  */