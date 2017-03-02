
library(MASS)
library(gtools)


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# predefined functions
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
cor2cov <- function( corr, std)  {
  return( diag(std) %*% corr %*% diag(std) )
}

lmvnorm <- function( x, mu, sig, d)	{
	return(-0.5*( d*log(2*pi) + determinant(sig)$modulus + 
			t(x-mu) %*% (chol2inv(chol(sig)) %*% (x-mu)) ))
}

lmult <- function( x, p)	{
	return( lfactorial(sum(x))-sum(lfactorial(x)) + x %*% p )
}
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# simulating 3 distinct mix multivariate normal, categorical and 
# multinomial sets of data
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

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
dim(dat_sim)
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo



original_clusters <- dat_sim[,1]
dat <- data.frame( dat_sim[,2:3], as.factor(dat_sim[,4]), apply(dat_sim[,5:9],2,as.integer) )
colnames(dat) <- paste0('x',1:8)



N <- nrow(dat)


#number of clusters
K <- 3

#identify classes of each variable
ipt_class <- unlist(lapply(dat,class))
idx_ftr <- which(ipt_class %in% 'factor')
idx_mult <- list(4:6,7:8)
idx_num <- setdiff(1:ncol(dat),c(idx_ftr,unlist(idx_mult)))

#continuous variable data
X_n <- dat[,idx_num]
plot(X_n)

#1 hot encoding of factor variables
if(length(idx_ftr) > 1){
	ipt_ftr_lvls <- sapply(dat[,idx_ftr], levels)
	Xc_l <- mapply(	function(x,y) t(do.call('rbind', lapply(y, function(z) x == z))), 
					dat[,idx_ftr], ipt_ftr_lvls 
					)  
	Xc_l <- mapply( function(x,y,z) {colnames(y) <- paste0( x, '_', z); return(y)} , 
					colnames(dat[,idx_ftr]), Xc_l, ipt_ftr_lvls)
} else {	
	ipt_ftr_lvls <- levels(dat[,idx_ftr])
	Xc_l <- list(sapply(ipt_ftr_lvls, function(x) dat[,idx_ftr] == x))
}

#multinomial variable data
if(length(idx_mult)>1){
	Xm_l <- lapply(idx_mult,function(x) dat[,x])
} else	{
	Xm_l <- list(dat[,idx_mult])
}


#initialize latent class
Z_t <- data.frame( rdirichlet(N, rep(1,K)), stringsAsFactors = FALSE)
Z <- matrix(0, nrow = N, ncol = K)
loglik_tm1 <- 0
ll=as.numeric()
i=0
t_0 <- proc.time()
while( max(apply(abs(Z_t - Z),2,max)) > 0.01)	{
	
	i = i + 1
	
	Z <- data.frame(Z_t, stringsAsFactors = FALSE)
	#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
	# Maximization
	#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
	#est - gaussian mean
	mu_k <- lapply(Z, function(z) apply(z*X_n,2,sum)/sum(z))

	#est - gaussian covariance
	XmMU_ik <- lapply(mu, function(x) data.frame(t(X_n) - x,stringsAsFactors = F) )
	XmMU_sq_ik <- lapply(XmMU_ik, function(x) lapply(x, function(y) y %*% t(y)))
	sig_sq_k <- mapply(function(q,r) apply(mapply( function(y,z) z*y, q, as.list(r)),1,sum) / sum(r), 
					   XmMU_sq_ik, Z )
	sig_sq_k <- lapply(data.frame(sig_sq_k), function(x) matrix(x, nrow=2,ncol=2))

	#est - categorical probabilities
	Pc_kl <- lapply( Z, function(z) lapply( Xc_l, function(x) apply(z*x,2,sum)/sum(apply(z*x,2,sum))) )
	
	#est - categorical probabilities
	Pm_kl <- lapply( Z, function(z) lapply( Xm_l, function(x) apply(z*x,2,sum)/sum(apply(z*x,2,sum))) )
	
	#est - mixture probabilities
	PI <- apply(Z,2,sum) / sum(apply(Z,2,sum))
	#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

	#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
	# Expectations
	#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
	lmvnorm_ik <- data.frame(mapply( function(p,r) { 
										d <- length(p)
										apply( X_n, 1, function(x) lmvnorm(x,p,r,d))
									 },
								    mu_k, sig_sq_k),
							 stringsAsFactors = F)	
	
	lcat_ikl <- lapply(lapply(Pc_kl,function(x) lapply(x,log)), 
					function(p) mapply( function(p_2,r) apply(apply(r,1,function(x) x*p_2),2,sum), 
										p, Xc_l)
					)
								
	lmult_ikl <- lapply(lapply(Pm_kl,function(x) lapply(x,log)), 
					function(p) mapply( function(p_2,r) apply(r,1,function(x) lmult(x,p_2)),
									    p, Xm_l)
					)
					
	Z_t <- mapply( function(a1,a2,a3,a4) rep(a1,N) + a2 + apply(a3,1,sum) + apply(a4,1,sum), 
				   as.list(log(PI)), lmvnorm_ik, lcat_ikl, lmult_ikl)
	expZ_t <- exp(Z_t)
	Z_t <- expZ_t / apply(expZ_t,1,sum)
	#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
	
	loglik <- mapply( function(a1,a2,a3,a4,a5) a1*(rep(a2,N) + a3 + apply(a4,1,sum) + apply(a5,1,sum)), 
						data.frame(Z), as.list(log(PI)), lmvnorm_ik, lcat_ikl, lmult_ikl	)
	loglik <- sum(apply(loglik,1,sum))
	
	ll <- c(ll,loglik)
	print(i)
	print(loglik,digits=20)
	print(loglik>loglik_tm1,digits=20)
	flush.console()
	loglik_tm1 <- loglik
}
t_N <- proc.time()
t_N-t_0


clust_output <- apply(Z_t, 1, function(x) which(x == max(x)))
table(clust_output)
table(original_clusters)
sum(clust_output == original_clusters)/length(clust_output)


mu_k;mu
sig_sq_k;varc
Pc_kl;P_c
Pm_kl;P_m


plot(X_n[clust_output==1,], 
	 xlim = c(min(X_n[,1]),max(X_n[,1])), 
	 ylim = c(min(X_n[,2]),max(X_n[,2])), col = 'blue' )
points(X_n[clust_output==2,],col='red')
points(X_n[clust_output==3,],col='lightblue')
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# system.time({
    # # Do something that takes time
    # lapply(sig_sq_k , function(x) chol2inv(chol(x)))
	# lapply(sig_sq_k , function(x) solve(x))
	# lapply(sig_sq_k , function(x) log(det(x)))
	# lapply(sig_sq_k , function(x) determinant(x)$modulus)
# })

# sum(eigen(siq_sq_k$X1)$values<0)
