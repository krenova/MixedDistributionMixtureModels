library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
sourceCpp("mfmmCppTest.cpp")

# log the 2 mistakes i made.

library(gtools)

N <- 10000
Z <- rdirichlet(N,rep(1,5))
X <- matrix(cbind(rnorm(N,5,2),rnorm(N,1,2)), nrow = N, ncol=2)


#-----------------------------------------------------
# Benchmarking the mean function
R_mu_k <- function(Z,X,M,K_clust)	{
	t(matrix(unlist(lapply(data.frame(Z), function(z) apply(z*X,2,sum)/sum(z))),nrow=M,ncol=K_clust))
}

microbenchmark(
	meanEM(X,Z,N,2,5),
	R_mu_k(Z,X,2,5)
)
all.equal(meanEM(X,Z,N,2,5), R_mu_k(Z,X,2,5))
#-----------------------------------------------------


#-----------------------------------------------------
# Benchmarking the covariance function
R_sig_k <- function(mu_k, X_n,Z)	{
    #est - gaussian covariance
	XmMU_ik <- lapply(data.frame(t(mu_k)), function(x) t(t(X_n) - x) )
	XmMU_sq_ik <- lapply(XmMU_ik, function(x) apply(x, 1, function(y) y %*% t(y)))
	sig_sq_k <- mapply(function(y,z) t(y)*z / sum(z), XmMU_sq_ik, data.frame(Z) ,SIMPLIFY = FALSE)
	sig_sq_k <- lapply(sig_sq_k,function(x) apply(x,2,sum))
	sig_sq_k <- lapply(sig_sq_k, function(x) matrix(x, nrow=2,ncol=2))
	return(sig_sq_k)
}
mu <- meanEM(X,Z,N,2,5)

microbenchmark(
	covEM1(X,Z,mu,N,2,5),
	covEM2(X,Z,mu,N,2,5),
	covEM3(X,Z,mu,N,2,5),
	covEM32(X,Z,mu,N,2,5),
	covEM4(X,Z,mu,N,2,5),
	covEM5(X,Z,mu,N,2,5),
	#R_sig_k(mu,X,Z),
	times = 500
)
microbenchmark(
	covEM_diag(X,Z,mu,N,2,5),
	covEM_diag2(X,Z,mu,N,2,5),
	times = 100
)


#-----------------------------------------------------
# Benchmarking the multinomial probabilities
idx_m <- list(4:6,7:8)
Xm_l <- lapply(idx_m,function(x) dat_sim[,x])


RmultEM <- function(Xc_l,Y)	{
	lapply( data.frame(Y), function(z) lapply( Xc_l, function(x) {
      P <- apply(z*x,2,sum)
      return(P/sum(P))
	}))
}
microbenchmark(
	multEM(Xm_l,Z[1:nrow(dat_sim),], length(Xm_l), nrow(dat_sim), sapply(Xm_l,ncol), 5),
	RmultEM(Xm_l,Z[1:nrow(dat_sim),])
)

#-----------------------------------------------------
# Benchmarking the cluster probabilities
Rclust_prob <-  function(input)	{
	temp <- apply(input,2,sum) 
	return(temp/sum(temp))
}
microbenchmark(
	clust_prob(Z,N,5),
	Rclust_prob(Z)
)

#-----------------------------------------------------
# Benchmarking multivariate gaussian log likelihood
N <- 1000
Z <- rdirichlet(N,rep(1,5))
X <- matrix(cbind(rnorm(N,5,2),rnorm(N,1,2)), nrow = N, ncol=2)

mu <- meanEM(X,Z,N,2,5)
sig <- covEM1(X, Z, mu, N, 2, 5)
mu_k <- data.frame(t(mu))
sig_k <- sig

lgaussian <- function(X_n,mu,sig)	{
	lmvnorm <- function( x, mu, sig, d)  {
		return( -0.5*( d*log(2*pi) + determinant(sig)$modulus + 
                t(x-mu) %*% (chol2inv(chol(sig)) %*% (x-mu)) ))
	}
	return(data.frame(mapply( function(p,r) { 
							d <- length(p)
							apply( X_n, 1, function(x) lmvnorm(x,p,r,d))
						},
						mu, sig),
				stringsAsFactors = F
				))
}

microbenchmark(
	lgaussian(X, mu_k,sig_k),
	lmvgauss( X, mu, sig, N, 2, 5)
)

#-----------------------------------------------------
# Benchmarking categorical and multinomial log likelihood

P <- multEM(Xm_l,Z[1:nrow(dat_sim),], length(Xm_l), nrow(dat_sim), sapply(Xm_l,ncol), 5)
Pc_kl <- RmultEM(Xm_l,Z[1:nrow(dat_sim),])

clcat <- function(x) lcat(Xm_l, P, 2, nrow(dat_sim), c(3,2), 5)

rlcat <- function() { lapply(lapply(Pc_kl,function(x) lapply(x,log)), 
       function(p) mapply( function(p_2,x) x %*% p_2, p, Xm_l)	
		)}

microbenchmark(
	clcat(),
	rlcat()
)

clmult_out <- function() clmult(Xm_l, P, 2, nrow(dat_sim), c(3,2), 5)
rlmult <- function() {
			 lapply(lapply(Pc_kl,function(x) lapply(x,log)), 
                          function(p) mapply( function(p_2,r) apply(r,1,function(x) lmult(x,p_2)),
                                              p, Xm_l)
						)}
microbenchmark(
	clmult_out(),
	rlmult(),
	times = 50
)

#-----------------------------------------------------
# Benchmarking full likelihood calculations log likelihood
mu <- meanEM(as.matrix(X_n),as.matrix(Z_t),N,len_num,5)
sig <- covEM1(as.matrix(X_n),as.matrix(Z_t),mu,N,len_num,5)
pm <- multEM(lapply(Xm_l,as.matrix),as.matrix(Z_t), 2, N, sapply(Xm_l,ncol), 5)
pc <- multEM(lapply(Xc_l,as.matrix),as.matrix(Z_t), 1, N, sapply(Xc_l,ncol), 5)

loggauss <- lmvgauss_k( as.matrix(X_n), mu, sig, N, len_num, 5)
logcat <- lcat_k(lapply(Xc_l,as.matrix), pc, 1, N, 5, 5)
logmult <- lmult_k(lapply(Xm_l,as.matrix), pm, 2, N, c(3,2), 5)


lmulconst_x <- lmultconst(lapply(Xm_l,as.matrix),2,N,c(3,2))
logmult2 <- lmult_k2(lapply(Xm_l,as.matrix), lmulconst_x, pm, 2, N, c(3,2), 5)


microbenchmark(
lmult_k(lapply(Xm_l,as.matrix), pm, 2, N, c(3,2), 5),
lmult_k2(lapply(Xm_l,as.matrix), lmulconst_x, pm, 2, N, c(3,2), 5)
)

test <- loglik_k(loggauss, lcat_ikl, logmult, 1, 2, N, 5)
test2 <- mapply( function(a1,a2,a3,a4) rep(a1,N) + a2 + apply(a3,1,sum) + apply(a4,1,sum), 
                        as.list(log(PI)), lmvnorm_ik, lcat_ikl, lmult_ikl)

						
test2 <- mapply( function(a1,a2,a3,a4) rep(a1,N) + a2 + apply(a3,1,sum) + apply(a4,1,sum), 
                        as.list(log(PI)), lmvnorm_ik, lcat_ikl, lmult_ikl)
						
a <- apply(loggauss,2,sum)
b <- apply(lmvnorm_ik,2,sum)		
a-b				

						
a <- unlist(lapply(lapply(lmult_ikl,function(x) apply(x,1,sum))	,sum))				
b <- apply(logmult,2,sum)
a-b

a <- unlist(lapply(lapply(lcat_ikl,function(x) apply(x,1,sum))	,sum))	
b <- apply(logcat,2,sum)
a-b

a <- clust_prob(as.matrix(Z_t), N, 5)
b <- apply(Z,2,sum) / sum(apply(Z,2,sum))
a-b

a <- loglik_k(loggauss, logcat, logmult, clust_prob(as.matrix(Z_t), N, 5) , N, 5)
b <- mapply( function(a1,a2,a3,a4) rep(a1,N) + a2 + apply(a3,1,sum) + apply(a4,1,sum), 
                        as.list(log(PI)), lmvnorm_ik, lcat_ikl, lmult_ikl)
apply(a-b,2,sum)
microbenchmark(
loglik_k(loggauss, logcat, logmult, clust_prob(as.matrix(Z_t), N, 5) , N, 5),
mapply( function(a1,a2,a3,a4) rep(a1,N) + a2 + apply(a3,1,sum) + apply(a4,1,sum), 
                        as.list(log(PI)), lmvnorm_ik, lcat_ikl, lmult_ikl)
)

prob <-function()	{ 
	Z_t <- a - apply(a,1,function(x) {
						x_max <- max(x)
						return(x_max + log(sum(exp(x-x_max))))
					   })
	return(exp(Z_t))
}

microbenchmark(
data_prob(a,N,5),
prob()
)
apply(data_prob(a,N,5) - prob(),2,sum)
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# simulating 3 distinct mix multivariate normal, categorical and 
# multinomial sets of data
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
library(gtools)

# N <- 10000
# Z <- rdirichlet(N,rep(1,5))
# X <- matrix(cbind(rnorm(N,5,2),rnorm(N,1,2)), nrow = N, ncol=2)
library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
sourceCpp("C:/Users/ACER/Documents/R Scripts/MFMM/mfmmCppTest.cpp")
setwd('C:/Users/ACER/Documents/R Scripts/MFMM')
source('MFMM.R')
mu <- list(c(1,2),c(4,6),c(0,2))
varc <- list( cor2cov(matrix(c(1,0.8,0.8,1),ncol=2,nrow=2),c(1,2.4)),
              cor2cov(matrix(c(1,0.5,0.5,1),ncol=2,nrow=2),c(1,1)),
              cor2cov(matrix(c(1,0.2,0.2,1),ncol=2,nrow=2),c(3,2))
)

P_c <- list(c(0.2,0.1,0.5,0.1,0.1),c(0.4,0.4,0.05,0.1,0.05),c(0.1,0.1,0.1,0.1,0.6))

P_m <- list(list(c(0.9,0.1,0.1),c(0.5,0.5)),list(c(0.9,0.1,0.1),c(0.5,0.5)),list(c(0.2,0.4,0.4),c(0.2,0.8)))

dat_sim <- list()
for(i in 1:3) {
  N_g <- round(runif(1,1000,3000))
  dat_sim <- c(dat_sim, list( cbind(rep(i,N_g),
                                    mvrnorm(N_g, mu[[i]], varc[[i]]),
                                    apply( t(rmultinom(N_g, 1, P_c[[i]])),1,function(x) which(x>0)),
                                    t(sapply(1:N_g, function(x) rmultinom(1, round(runif(1, 30, 300)), P_m[[i]][[1]]))),
                                    t(sapply(1:N_g, function(x) rmultinom(1, round(runif(1, 5, 30)), P_m[[i]][[2]])))
  )
  )
  )
}
dat_sim <- do.call('rbind',dat_sim)
original_clusters <- dat_sim[,1]
dat <- data.frame( dat_sim[,2:3], as.factor(dat_sim[,4]), apply(dat_sim[,5:9],2,as.integer) )
colnames(dat) <- paste0('x',1:8)
idx_m <- list(4:6,7:8)


K = 5
multinomial_indices <- idx_m
est_cov = TRUE
N <- nrow(dat)
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo