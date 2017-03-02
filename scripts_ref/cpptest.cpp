#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// using namespace arma;

#include <math.h>
#include <valarray>
#include <vector>

//[[Rcpp::export]]
NumericVector rowSumsC(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow);

  for (int i = 0; i < nrow; i++) {
    double total = 0;
    for (int j = 0; j < ncol; j++) {
      total += x(i, j);
    }
    out[i] = total;
  }
  return out;
}

// flag = arma::inv_sympd(sig_inv,sig_k);
		// if(!flag) {
			// throw std::runtime_error("matrix inversion failed");
			// Rcout << "Please repeat the function call or reduce the number of clusters if the problem persists."<< std::endl;
		// }

//[[Rcpp::export]]
arma::mat TEST(arma::mat x)	{
	arma::mat out = arma::zeros(x.n_rows,x.n_cols);
	bool flag = arma::inv_sympd(out,x);
	if(!flag) {
		int k = 5;
		Rcout << "Note: Covariance matrix for cluster " << k << " could not be inverted. Please" << std::endl;
		Rcout << "\t[1]repeat the function call or \n\t[2]reduce the number of clusters \n      if the problem persists." << std::endl;
		throw std::runtime_error("matrix inversion failed");
	}
	return out;
}


// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Matrix Inverse
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//[[Rcpp::export]]
arma::mat Cinv(arma::mat& x)	{
	return arma::inv(x);
}

//[[Rcpp::export]]
arma::mat Cinv_sympd(arma::mat x)	{
	return arma::inv_sympd(x);
}

//[[Rcpp::export]]
arma::mat Cinv_chol(arma::mat& x)	{
	return arma::inv(arma::trimatu(arma::chol(x)));
}

// ***************
//[[Rcpp::export]]
NumericMatrix Rcppinv_chol(NumericMatrix x)	{
	return wrap(arma::inv(arma::trimatu(arma::chol(as<arma::mat>(x)))));
}

//[[Rcpp::export]]
NumericMatrix Rcppinv_chol2(NumericMatrix x)	{
	return wrap(arma::trimatu(arma::chol(as<arma::mat>(x))).i());
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Multinomial Constant
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//[[Rcpp::export]]
NumericVector Clmult(NumericMatrix x)	{
	int nrow = x.nrow(), ncol = x.ncol();
	NumericVector out(nrow);
	for(int i = 0; i < nrow; i++)	{
		double i_out1 = 0, i_out2 = 0;
		for(int j = 0; j < ncol; j++)	{
			i_out1 += x(i,j);
			i_out2 += lgamma(x(i,j)+1);
		}
		out[i] = lgamma(i_out1+1) - i_out2;
	}
	return out;
}

// ***************
//[[Rcpp::export]]
arma::vec Clmult2(arma::mat x)	{
	int nrow = x.n_rows, ncol = x.n_cols;
	arma::vec out = arma::zeros(nrow);
	for(int i = 0; i < nrow; i++)	{
		double i_out1 = 0, i_out2 = 0;
		for(int j = 0; j < ncol; j++)	{
			i_out1 += x(i,j);
			i_out2 += lgamma(x(i,j)+1);
		}
		out[i] = lgamma((i_out1)+1) - i_out2;
	}
	return out;
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// element wise vector product
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//[[Rcpp::export]]
arma::vec multiply_arma(arma::vec x, arma::vec y)	{
	int nrow = x.n_rows;
	arma::vec out = arma::zeros(nrow);
	for(int i = 0; i < nrow; i++)	{
		out[i] = x[i]*y[i];
	}
	return out;
}

//[[Rcpp::export]]
NumericVector multiply_Rcpp(NumericVector x, NumericVector y)	{
	const int nrow = x.size();
	NumericVector out(nrow);
	for(int i = 0; i < nrow; i++)	{
		out[i] = x[i]*y[i];
	}
	return out;
}

//[[Rcpp::export]]
arma::vec multiply_arma2(arma::vec x, arma::vec y)	{
	return x % y;
}

// ***************
//[[Rcpp::export]]
NumericVector multiply_Rcpp2(NumericVector x, NumericVector y)	{
	return x * y;
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// element wise vector to matrix column product
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

//[[Rcpp::export]]
arma::mat multiply_matArma(arma::mat x, arma::vec y)	{
	int nrow = x.n_rows, ncol = x.n_cols;
	arma::mat out = arma::zeros(nrow,ncol);
	for(int j = 0; j < ncol; j++)	{
		out.col(j) = x.col(j) % y;
	}
	return out;
}

//[[Rcpp::export]]
NumericMatrix multiply_matRcpp(NumericMatrix x, NumericVector y)	{
	int nrow = x.nrow(), ncol = x.ncol();
	NumericMatrix out( nrow, ncol);
	for(int j = 0; j < ncol; j++)	{
		out.column(j) = x.column(j) * y;
	}
	return out;
}

// ***************
//[[Rcpp::export]]
NumericMatrix multiply_matMan(NumericMatrix x, NumericVector y)	{
	int nrow = x.nrow(), ncol = x.ncol();
	NumericMatrix out( nrow, ncol);
	for(int j = 0; j < ncol; j++)	{
		for(int i = 0; i < ncol; i++)	{
			out(i,j) = x(i,j) * y(i);
		}
	}
	return out;
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// element wise vector to matrix column product
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

//[[Rcpp::export]]
double sumSugar(NumericVector x)	{
	return sum(x);
}

//[[Rcpp::export]]
double sumArma(NumericVector x)	{
	return arma::accu(as<arma::vec>(x));
}

//[[Rcpp::export]]
double sumArma2(arma::vec x)	{
	return arma::accu(x);
}

//[[Rcpp::export]]
double sumRcpp(NumericVector x)	{
	const int size = x.size();
	double out = 0;
	for(int i; i < size; i++)	{
		out += x[i];
	}
	return out;
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// estimation of mean
// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
NumericMatrix mu_k(NumericMatrix x, NumericMatrix z, 
				   int nrow, int ncol, int K_clust)	{
	NumericMatrix out( K_clust, ncol);
	for(int k; k < K_clust; k++)	{
		for(int j = 0; j < ncol; j++)	{
			for(int i = 0; i < nrow; i++)	{
				out(k,j) += x(i,j) * z(i,k);
			}
		}
	}
	return out;
}


double * rDirichlet(int& N, std::vector<double>& alp)	{
	
	int K = alp.size();
	double* out = new double[N*K];
	double sum_k = 0;
	
	for (int i = 0; i < N; i++)	{
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

std::vector<double> rDirichlet2(int& N, std::vector<double>& alp)	{
	
	int K = alp.size();
	std::vector<double> out (N*K,1);
	double sum_k = 0;
	
	for (int i = 0; i < N; i++)	{
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


//[[Rcpp::export]]
NumericMatrix cDirichlet(int& N, NumericVector& alp) {
	
	int K = alp.size();
	std::vector<double> alpvec = as< std::vector<double> >(alp);
	double *in;
	NumericMatrix out(N,K);
	
	in = rDirichlet(N, alpvec);
	for(int i = 0; i < N; i++)  {
        for(int j = 0; j < K; j++)  {
             out(i,j) = in[j + i*K];
        }
	}
	delete[] in;
	return out;
}

//[[Rcpp::export]]
NumericMatrix cDirichlet2(int& N, NumericVector& alp) {
	
	int K = alp.size();
	std::vector<double> alpvec = as< std::vector<double> >(alp);
	std::vector<double> in = rDirichlet2(N, alpvec);
	NumericMatrix out(N,K);
		
	for(int i = 0; i < N; i++)  {
        for(int j = 0; j < K; j++)  {
             out(i,j) = in[j + i*K];
        }
	}
	return out;
}
  
//[[Rcpp::export]]
NumericMatrix cDirichlet22(int& N, NumericVector& alp) {
	
	NumericMatrix out = cDirichlet2(N, alp);
	return out;
}
// Note: can embed R code in the following manner
// /*** R
// library(microbenchmark)
// x <- runif(1e5)
// microbenchmark(
//   mean(x),
//   meanC(x)
// )
// */