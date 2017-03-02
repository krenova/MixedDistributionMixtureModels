
setwd('D:/Project - Research')


source('MFMM.R')


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


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# test run the MFMM function
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

original_clusters <- dat_sim[,1]
dat <- data.frame( dat_sim[,2:3], as.factor(dat_sim[,4]), apply(dat_sim[,5:9],2,as.integer) )
colnames(dat) <- paste0('x',1:8)

idx_m <- list(4:6,7:8)

clust_output <- MFMM( dat, K = 3, multinomial_indices = idx_m, est_cov = FALSE)
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# evaluating the MFMM results
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
mu
varc
P_c
P_m

sort(table(original_clusters),decreasing = T)
summary(clust_output)

#identify classes of each variable
ipt_class <- unlist(lapply(dat,class))
idx_ftr <- which(ipt_class %in% c('factor','logical'))
idx_mult <-list(4:6,7:8)
idx_num <- setdiff(1:ncol(dat),c(idx_ftr,unlist(idx_mult)))

#identify classes of each variable	
plot(dat[,idx_num][clust_output$clusters==1,], 
     xlim = c(min(dat[,idx_num][,1]),max(dat[,idx_num][,1])), 
     ylim = c(min(dat[,idx_num][,2]),max(dat[,idx_num][,2])), col = 'blue' )
points(dat[,idx_num][clust_output$clusters==2,],col='red')
points(dat[,idx_num][clust_output$clusters==3,],col='lightblue')
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo



source('MFMM_plot_tools.R')

# plot gaussian data
for ( i in 1:2)	{
  for (k in 1:3)	{
    hist_mfmm(dat[,i], dat[clust_output$clusters == k,i])
    readline()
  }
}
#plot categorical data
for (k in 1:3)	{
  barplot_mfmm(dat[,3], dat[clust_output$clusters == k,3], ylim = c(0,0.7))
  readline()
}
#plot multinomial data
for ( i in 1:2)	{
  for (k in 1:3)	{
    barplot_mfmm(dat[,idx_m[[i]]], dat[clust_output$clusters == k,idx_m[[i]]], ylim = c(0,1))
    readline()
  }
}


