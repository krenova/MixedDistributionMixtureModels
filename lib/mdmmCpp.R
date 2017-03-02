
library(Rcpp)
library(RcppArmadillo)
sourceCpp('lib/mdmmCore.cpp')

#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# Mixed Finite Mixture Model Function (Via EM Algorithm)
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
MDMM <- function(dat, K, multinomial_indices = integer(), tol = 0.01, est_cov = TRUE, attempts = 1, criterion = "BIC")  {
  
  #----------- [1] Algorithm Initialization -----------
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
    X_n <- as.matrix(dat[,idx_num])
  }  else  {
    X_n <- matrix()
  }
  
  #1 hot encoding of factor variables
  len_ftr <- length(idx_ftr)
  if(len_ftr > 1){
    ipt_ftr_lvls <- sapply(dat[,idx_ftr], levels)
    Xc_l <- mapply(	function(x,y) t(do.call('rbind', lapply(y, function(z) x == z))), 
                    dat[,idx_ftr], ipt_ftr_lvls )  
    Xc_l <- mapply( function(x,y,z) {colnames(y) <- paste0( x, '_', z); return(y)} , 
                    colnames(dat[,idx_ftr]), Xc_l, ipt_ftr_lvls )
  } else if (len_ftr > 0) {	
    ipt_ftr_lvls <- levels(dat[,idx_ftr])
    Xc_l <- list(sapply(ipt_ftr_lvls, function(x) dat[,idx_ftr] == x))
  } else {
    Xc_l <- list()
  }
  Xc_l <- lapply(Xc_l,"*",1)
  if(len_ftr>0) Xc_ncol <- sapply(Xc_l,ncol) else Xc_ncol <- 0
  
  #multinomial variable data
  len_mult <- length(idx_mult)
  if(is.list(idx_mult) & len_mult>1){
    Xm_l <- lapply(idx_mult,function(x) as.matrix(dat[,x]))
  } else if(len_mult>0)	{
    Xm_l <- list(as.matrix(dat[,idx_mult]))
  } else {
    Xm_l <- list()
  }
  if(len_ftr>0) Xm_ncol <- sapply(Xm_l,ncol) else Xm_ncol <- 0
  
  #----------- [2] Algorithm Execution -----------
  output = mdmmCore(X_n, len_num,
					Xc_l, Xc_ncol, len_ftr,
					Xm_l, Xm_ncol, len_mult,
					N, K, tol, est_cov, attempts, criterion)
  
  if( length(output)>0 ) {
	  #----------- [3] Output Preparation -----------
	  #name GOF fields
	  names(output$GOF) <- c("AIC", "BIC")
	  
	  #name data_prob output cluster groups
	  colnames(output$data_prob) <- paste0('clust', 1:K)
	  
	  #name cluster probabilities
	  names(output$cluster_prob) <- paste0('clust', 1:K)
	  
	  #name gaussian params or delete output if not used
	  if( len_num > 0 ) {
		colnames(output$gaussian_mean) <- colnames(X_n)
		rownames(output$gaussian_mean) <- paste0('clust', 1:K)
		names(output$gaussian_var) <- paste0('clust', 1:K)
		output$gaussian_var <- lapply(output$gaussian_var,
									function(x) {
										colnames(x) <- colnames(X_n)
										rownames(x) <- colnames(x)
										return(x)
									})
	  } else {
		output <- output[!( names(output) %in% c("gaussian_mean", "gaussian_var") )]
	  }
	  
	  #name catagorical params or delete output if not used
	  if( len_ftr > 1 ) {
		names(output$categorical_prob) <- paste0("Xc_", names(Xc_l))
		output$categorical_prob <- mapply(function(x,y) {
										colnames(x) <- colnames(y)
										rownames(x) <- paste0('clust', 1:K)
										return(x)
									 },
									output$categorical_prob, 
									Xc_l )
	  } else if( len_ftr == 1 ) {
		names(output$categorical_prob) <- names(Xc_l)
		colnames(output$categorical_prob[[1]]) <- colnames(Xc_l[[1]])
		rownames(output$categorical_prob[[1]]) <- paste0('clust', 1:K)
	  } else {
		output <- output[!( names(output) %in% c("categorical_prob") )]
	  }
	  
	  #name multinomial params or delete output if not used
	  if( len_mult > 1 ) {
		names(output$multinomial_prob) <- paste0("Xm_", 1:len_mult)
		output$multinomial_prob <- mapply(function(x,y) {
										colnames(x) <- colnames(y)
										rownames(x) <- paste0('clust', 1:K)
										return(x)
									 },
									output$multinomial_prob, 
									Xm_l )
	  } else if ( len_mult == 1) {
		colnames(output$multinomial_prob[[1]]) <- colnames(Xm_l[[1]])
		rownames(output$multinomial_prob[[1]]) <- paste0('clust', 1:K)
	  } else{
		output <- output[!( names(output) %in% c("multinomial_prob") )]
	  }
	  
	  #----------- [4] Finalize Output -----------
	  attr(output, "class") <- c("mdmm","list")
	  attr(output, "truncated") <- 'data_prob'
	  return(output)
	} else{
		message('None of the evaluation attempts have been successful.')
		message('Please reconfigure the function parameters.')
	}
}

#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# To enable use of summary function on output_summary
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# To prevent printing of all data when MDMM output is called
# print.mdmm <- function (x) {
  # hid <- attr(x, "hidden")
  # print(x[!names(x) %in% hid])
# }
print.mdmm <- function (x) {
  trnctd <- attr(x, "truncated")
  x[[trnctd]] <-  rbind( x[[trnctd]][1:20,], rep("...truncated...",ncol(x[[trnctd]])) )
  print(x[1:length(x)], quote=FALSE)
}

summary <- function (x, ...) {
  UseMethod("summary", x)
}

summary.mdmm <- function(x, ...) {
  cat('\n\nGoodness of Fit Statistics:\n')
  print(x$GOF, digits = max(5, getOption("digits") - 3))
  cat('\nCluster Membership Counts:')
  print(table(apply(x$data_prob,1,function(x) which(x==max(x)))))
}
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
