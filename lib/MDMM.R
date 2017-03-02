
#**********************************************************************************
#**********************************************************************************
#**********************************************************************************
# Original pure R prototype of the MDMM algorithm. 
#	- used to illustrate the difference in speed between the R and Rcpp 
#	  implementation
#	- also contains thenon-monotonic likelihood bug which needs to be fixed 
#**********************************************************************************
#**********************************************************************************
#**********************************************************************************


library(MASS)
library(gtools)


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# predefined functions
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
cor2cov <- function( corr, std)  {
  return( diag(std) %*% corr %*% diag(std) )
}

# log multivariate normal
lmvnorm <- function( x, mu, sig, d)  {
  return(-0.5*( d*log(2*pi) + determinant(sig)$modulus + 
                  t(x-mu) %*% (chol2inv(chol(sig)) %*% (x-mu)) ))
}

# log multinomial function 
# note that numerical overflow issues have been observed in the usage of lfactorial
# which results in the occasional reduction of MDMM log-likelihood calculations.
lmult <- function( x, p)  {
  return( lfactorial(sum(x))-sum(lfactorial(x)) + x %*% p )
}
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# Mixed Finite Mixture Model Function (Via EM Algorithm)
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
MDMM <- function(dat, K, multinomial_indices = integer(), tol = 0.01, est_cov = TRUE)  {
  
  #number of data points
  N <- nrow(dat)
  
  #identify classes of each variable
  ipt_class <- character()
  for(i in 1:ncol(dat)) ipt_class[i] <- class(dat[,i])
  idx_ftr <- which(ipt_class %in% c('factor','logical'))
  idx_mult <- multinomial_indices
  idx_num <- setdiff(1:ncol(dat),c(idx_ftr,unlist(idx_mult)))
  
  #coerce logical data fields into factor
  if ( sum(ipt_class %in% 'logical') > 0 )	{
	  dat[,which(ipt_class %in% 'logical')] <- lapply(data.frame(dat[,which(ipt_class %in% 'logical')],stringsAsFactors=F),
													  as.factor)
	}
  
  #continuous variable data
  len_num <- length(idx_num)
  if(len_num>0) {
    X_n <- dat[,idx_num]
  }  else  {
    lmvnorm_ik <- list(t(t(rep(0,N))))
  }
  
  #1 hot encoding of factor variables
  len_ftr <- length(idx_ftr)
  if(len_ftr > 1){
    ipt_ftr_lvls <- sapply(dat[,idx_ftr], levels)
    Xc_l <- mapply(	function(x,y) t(do.call('rbind', lapply(y, function(z) x == z))), 
                    dat[,idx_ftr], ipt_ftr_lvls 
    )  
    Xc_l <- mapply( function(x,y,z) {colnames(y) <- paste0( x, '_', z); return(y)} , 
                    colnames(dat[,idx_ftr]), Xc_l, ipt_ftr_lvls)
  } else if (len_ftr > 0) {	
    ipt_ftr_lvls <- levels(dat[,idx_ftr])
    Xc_l <- list(sapply(ipt_ftr_lvls, function(x) dat[,idx_ftr] == x))
  } else {
    lcat_ikl <- list(t(t(rep(0,N))))
  }
  Xc_l <- lapply(Xc_l,"*",1)
  
  #multinomial variable data
  len_mult <- length(idx_mult)
  if(is.list(idx_mult) & len_mult>1){
    Xm_l <- lapply(idx_mult,function(x) dat[,x])
  } else if(!is.list(idx_mult))	{
    Xm_l <- list(dat[,idx_mult])
  } else {
    lmult_ikl <- list(t(t(rep(0,N))))
  }
  
  #initialize latent class
  Z_t <- data.frame( rdirichlet(N, rep(1,K)), stringsAsFactors = FALSE)
  Z <- matrix(0, nrow = N, ncol = K)
#**********************************************************
#  min_z <- max_z <- avg_z <- matrix(0, nrow = 0, ncol = K)
#**********************************************************
  #while loop initialization
  eps <- 1e99
  loglik <- numeric()
  i=0
  while( eps > 0.01)	{
    
    i = i + 1
    Z <- data.frame(Z_t, stringsAsFactors = FALSE)
    sum_zk <- apply(Z,2,sum)
    #--------------------------- Maximization (M-Step) --------------------------
    #est - gaussian parameters
    if( len_num>0  )	{
      #est - gaussian mean
      mu_k <- lapply(Z, function(z) apply(z*X_n,2,sum)/sum(z))
	  
      if( est_cov )	{
        #est - gaussian covariance
        XmMU_ik <- lapply(mu_k, function(x) t(t(X_n) - x) )
        XmMU_sq_ik <- lapply(XmMU_ik, function(x) apply(x, 1, function(y) y %*% t(y)))
        sig_sq_k <- mapply(function(y,z) t(y)*z, XmMU_sq_ik, Z ,SIMPLIFY = FALSE)
        sig_sq_k <- mapply(function(x,z) apply(x,2,sum)/z,sig_sq_k,sum_zk)
        sig_sq_k <- lapply(data.frame(sig_sq_k), function(x) matrix(x, nrow=len_num,ncol=len_num))
      } else {
        #est - gaussian 0 correlation covariance matrix
        XmMU_ik <- lapply(mu_k, function(x) t(t(X_n) - x) )
        XmMU_sq_ik <- lapply(XmMU_ik, '^', 2)
        sig_sq_k <- mapply(function(y,z) y*z , XmMU_sq_ik, Z , SIMPLIFY = FALSE)
        sig_sq_k <- mapply(function(x,z) apply(x,2,sum)/z,sig_sq_k,sum_zk)
        sig_sq_k <- lapply(data.frame(sig_sq_k), diag)
      }
    }
    
    #est - categorical probabilities
    if( len_ftr>0 ) Pc_kl <- lapply( Z, function(z) lapply( Xc_l, function(x) {
      P <- apply(z*x,2,sum)
      P <- P/sum(z)
    }))
    
    #est - multinomial probabilities
    if( len_mult>0 ) Pm_kl <- lapply( Z, function(z) lapply( Xm_l, function(x) {
	  den <- z*apply(x,1,sum)
      P <- apply(z*x,2,sum)
      P <- P/sum(den)
    }) )
    
    #est - mixture probabilities
    PI <- apply(Z,2,sum) / sum(apply(Z,2,sum))
    
    #--------------------------- calculating the log-likelihood -----------------
    #gaussian log-likelihood
    if( len_num>0 )	{
		lmvnorm_ik <- data.frame(mapply( function(p,r) { 
											d <- length(p)
											apply( X_n, 1, function(x) lmvnorm(x,p,r,d))
											},
									mu_k, sig_sq_k),
								stringsAsFactors = F
								)
    }
    
    #categorical log-likelihood
    if( len_ftr>0 ) {
      lcat_ikl <- lapply(lapply(Pc_kl,function(x) lapply(x,log)), 
                         function(p) mapply( function(p_2,x) x %*% p_2, p, Xc_l)									  
      )
    }
    
    #multinomial log-likelihood
    if( len_mult>0 ) {
      lmult_ikl <- lapply(lapply(Pm_kl,function(x) lapply(x,log)), 
                          function(p) mapply( function(p_2,r) apply(r,1,function(x) lmult(x,p_2)),
                                              p, Xm_l)
      )
    }
    
    #--------------------------- Expectations (E-Step) --------------------------
    loglik_k <- mapply( function(a1,a2,a3,a4) rep(a1,N) + a2 + apply(a3,1,sum) + apply(a4,1,sum), 
                        as.list(log(PI)), lmvnorm_ik, lcat_ikl, lmult_ikl )
    Z_t <- loglik_k - apply(loglik_k,1,function(x) {
						x_max <- max(x)
						return(x_max + log(sum(exp(x-x_max))))
					   })
    Z_t <- exp(Z_t)
    
    #--------------------------- full log-likelihood ----------------------------
#********************************************************************
# min_z <- rbind(min_z,apply(loglik_k,2,min))
# max_z <- rbind(max_z,apply(loglik_k,2,max))
# avg_z <- rbind(avg_z,apply(loglik_k,2,mean))
#********************************************************************
    loglik <- c(loglik, sum(apply(Z*loglik_k,1,sum)))
    
    
    #--------------------------- end of iteration -------------------------------
    eps <- max(apply(abs(Z_t - Z),2,max))
    print(paste0('        Current iteration: ', i))
    print(paste0('|theta_tp1 - theta_t|_inf: ', eps))
    print(paste0('            loglikelihood: ', loglik[i]))
    flush.console()
	diff(loglik)
  }
  
  # identify cluster output (using mode)
  clust_output <- apply(Z_t, 1, function(x) which(x == max(x)))
  
  # rename estimated parameters
  names(mu_k) <- paste0('mu_', 1:K)
  names(sig_sq_k) <- paste0('sig_', 1:K)
  names(Pc_kl) <- paste0('Pc_', 1:K)
  names(Pm_kl) <- paste0('Pm_', 1:K)
  
  N_param <- (len_num + len_num*len_num +
                len_ftr*length(unlist(ipt_ftr_lvls)) +
                length(unlist(idx_mult)))*K
  
  AIC <- 2*(N_param -	loglik[i])
  AICc <- AIC + 2*N_param*(N_param+1)/(N-N_param-1)
  BIC <- -2*loglik[i]+N_param*log(N)
  
  output_summary <- list( clusters = clust_output,
                          cluster_assignment_probability = Z_t,
                          GOF = c(AIC = AIC, AICc = AICc, BIC = BIC),
                          Data_Attributes = list(N_datapoints = N,
                                                 N_parameters = N_param,
                                                 estimated_parameters = list( gaussian_mu = mu_k, 
                                                                              gaussian_sigma = sig_sq_k, 
                                                                              categorical_p = Pc_kl,
                                                                              multinomial_p = Pm_kl
																			)
												),
                          loglikelihood = loglik[length(loglik)],
                          loglike_profile = loglik
                          )
  class(output_summary) <- 'mdmm' 
  attr(output_summary, "hidden") <-c('clusters',
                                     'cluster_assignment_probability',
                                     'Data_Attributes')
  return(output_summary)
}

# To prevent printing of all data when MDMM output is called
print.mdmm <- function (x) {
  hid <- attr(x, "hidden")
  print(x[!names(x) %in% hid])
}
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# To enable use of summary function on output_summary
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
summary <- function (x, ...) {
  UseMethod("summary", x)
}

summary.mdmm <- function(x, ...) {
  cat('\n\nGoodness of Fit Statistics:\n')
  print(x$GOF, digits = max(5, getOption("digits") - 3))
  cat('\nCluster Membership Counts:')
  print(table(x$clusters))
}
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
