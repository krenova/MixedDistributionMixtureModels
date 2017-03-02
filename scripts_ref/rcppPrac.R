library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
sourceCpp("E:/mfmm/Algorithms - 20160714/cpptest.cpp")

library(gtools)

RRdiri <- function(N,K)	{
	out <- matrix(0,ncol=K,nrow = N)
	for( i in 1:N)	{
		temp <- rgamma(K,1,1);
		out[i,] <- temp/sum(temp);
	}
	return(out);
}

microbenchmark(
	cDirichlet(N,rep(1,K)),
	cDirichlet2(N,rep(1,K)),
	cDirichlet22(N,rep(1,K)),
	rdirichlet(N, rep(1,K)),
	RRdiri(N,K)
)

Posdef <- function (n, ev = runif(n, 0, 10)) {
	Z <- matrix(ncol=n, rnorm(n^2))
	decomp <- qr(Z)
	Q <- qr.Q(decomp) 
	R <- qr.R(decomp)
	d <- diag(R)
	ph <- d / abs(d)
	O <- Q %*% diag(ph)
	Z <- t(O) %*% diag(ev) %*% O
	return(Z)
}

A <- Posdef(n=10)

x <- y <- runif(100000*10,-100,100)
X <- matrix(round(abs(x)), nrow = 100000, ncol=10)

y_100000 <- runif(100000,-100,100)













microbenchmark(
	multiply_arma(x, y),
	multiply_Rcpp(x, y),	#3rd fastest (not too far off from top 2)
	multiply_arma2(x, y),
	multiply_Rcpp2(x, y), 	#fastest
	x*y						#2nd fastest (v close)
)

microbenchmark(
	multiply_matArma(X, y_100000),
	multiply_matRcpp(X, y_100000),
	multiply_matMan(X, y_100000), #fastest
	X*y_100000
)

microbenchmark(
	apply(X,1,function(x) lfactorial(sum(x)) - sum(lfactorial(x))),
	Clmult(X), #fastest
	Clmult2(X)
)

microbenchmark(
	solve(A),
	chol2inv(chol(A)),
	Cinv(A),
	Rcppinv_chol(A), #to be used
	Rcppinv_chol2(A), #uses .i() instead of inv(), very slow
	Cinv_chol(A), #fastest but in arma
	Cinv_sympd(A)
)

microbenchmark(
	sumSugar(y_100000),
	sumArma(y_100000),
	sumArma2(y_100000),
	sumRcpp(y_100000),
	times = 5000
)

