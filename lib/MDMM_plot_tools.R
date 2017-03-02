
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# general purpose functions used.
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
numsplit <- function(x) {
  if(x>=1)  {
    return(as.integer(unlist(strsplit(strsplit(as.character(x), '\\.')[[1]][1], ''))))
  } else if (x<1 & nchar(x)>1)	{
    return(strsplit(as.character(x),'')[[1]])
  } else {
    stop('input must be a number and must be positve!')
  }
}
round_D <- function(x) {
  options(scipen=999)
  num <- numsplit(x)
  N_num <- length(num)
  if(x>=1)	{	
    if( num[2]>=5 & N_num>2 ) {
      num[1] <- num[1] + 1
      num[2:length(num)] <- 0
    } else if (num[2]<5 & N_num>2) {
      num[2] <- 5
      num[3:length(num)] <- 0
    }
    return(as.integer(paste0(num,collapse='')))
  } else {
    idx_sig_fig <- which(!(num[-(1:2)] %in% '0'))[1]
    if(!is.na(num[5+3]))	{
      return(eval(parse(text = paste0('1e-',idx_sig_fig)))*(as.integer(num[idx_sig_fig+2])+1))
    } else	{
      return(x)
    }
  }
}

#wrap.labels(names(x_plot),20)
wrap.it <- function(x, len)
{ 
  sapply(x, function(y) paste(strwrap(y, len), 
                              collapse = "\n"), 
         USE.NAMES = FALSE)
}


# Call this function with a list or vector
wrap.labels <- function(x, len)
{
  if (is.list(x))
  {
    lapply(x, wrap.it, len)
  } else {
    wrap.it(x, len)
  }
}

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo



# #ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# # Function to prep data for plotting
# #ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# mdmm_plotdata <- function(dat, idx_m)  {

# #intitialization
# N <- nrow(dat)
# ipt_class <- unlist(lapply(dat,class))
# idx_ftr <- which(ipt_class %in% c('factor','logical'))
# if (is.list(idx_m)) idx_mult <- idx_m else idx_mult <- list(idx_m)
# idx_num <- setdiff(1:ncol(dat),c(idx_ftr,unlist(idx_mult)))

# N_ftr <- ifelse( !is.null(idx_ftr), length(idx_ftr), 0)
# N_mult <- ifelse(!is.null(idx_mult), length(idx_mult), 0)
# N_num <- ifelse( !is.null(idx_num), length(idx_num), 0)
# N_fields <- N_ftr + N_mult + N_num


# output <- list()
# #numeric plot handles
# if(N_num>0)	{
# plot_num <- list()
# for (i in 1:N_num)	{

# x_rge <- quantile(dat[,idx_num[i]], c(0.1,0.9))
# idx_rge <- (dat[,idx_num[i]] >= x_rge[1]) & (dat[,idx_num[i]] <= x_rge[2])
# x_plt <- dat[idx_rge,idx_num[i]]
# N_x_uqe <- length(unique(dat[idx_rge,idx_num[i]]))

# if(N_x_uqe < 20)	{
# x_table <- table(x_plt) / sum(table(x_plt))
# x_plot_bar <- list(data = x_table)
# x_plot_bar <- c(x_plot_bar, list( ylim = c(0,round_D(max(x_table)))) )
# class(x_plot_bar) <- 'mdmm_bar'
# plot_num <- c( plot_num, list(x_plot_bar) )
# }	else 	{
# x_plot_hist <- list( data = x_plt)
# x_plot_hist <- c(x_plot_hist, list( breaks = ifelse( N_x_uqe/50 < 50, 20, 50)) )
# class(x_plot_hist) <- 'mdmm_hist'
# plot_num <- c( plot_num, list(x_plot_hist) )
# }
# }
# output <- c(output, Num_Data = list(plot_num))
# }

# #categorical plot handles
# if(N_ftr>0)	{
# plot_ftr <- list()
# for (i in 1:N_ftr)	{
# x_table <- table(dat[,idx_ftr[i]]) / sum(table(dat[,idx_ftr[i]]))
# x_plot_bar <- list(data = x_table)
# x_plot_bar <- c(x_plot_bar, list( ylim = c(0,round_D(max(x_table)))) )
# class(x_plot_bar) <- 'mdmm_bar'
# plot_ftr <- c( plot_ftr, list(x_plot_bar) )
# }
# output <- c(output, Cat_Data = list(plot_ftr))
# }

# #multinomial plot handles
# if(N_mult>0)	{
# plot_mult <- list()
# for (i in 1:N_mult)	{
# x_table <- apply(dat[,idx_mult[[i]]],2,sum)
# x_table <- x_table/sum(x_table)
# x_plot_bar <- list(data = x_table)
# x_plot_bar <- c(x_plot_bar, list( ylim = c(0,round_D(max(x_table)))) )
# class(x_plot_bar) <- 'mdmm_bar'
# plot_mult <- c( plot_mult, list(x_plot_bar) )
# }
# output <- c(output, Mult_Data = list(plot_mult))
# }


# }
# #ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


# #ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# # Functions to plot the data created by dat_plot()
# #ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# plot <- function (x, ...) {
# UseMethod("plot", x)
# }

# plot.mdmm_bar <- function(x, ...) {
# barplot(x$data, ylim = x$ylim, ...)
# }

# plot.mdmm_hist <- function( x_base, x_cmpre = NULL, breaks = NULL, ...) {
# if( is.null(x_cmpre) )	{
# hist(x_base$data, 
# breaks = ifelse(is.null(breaks), x_base$breaks, breaks), 
# ...)
# } else {
# x_base_rg <- diff(range(x_base$data))
# n_base <- length(x_base$data)
# x_cmpre_rg <- diff(range(x_cmpre$data))
# n_cmpre <- length(x_cmpre$data)
# x_ttl <- c(x_base$data,x_cmpre$data)
# if ( (x_base_rg/n_base) > (x_cmpre_rg/n_cmpre) )	{
# plt1 <- plot( x_base , plot = F)
# breaks <- x_cmpre_rg/(x_base_rg/length(plt1$breaks))
# breaks <- ifelse(breaks>1, round(breaks), 1)
# plt2 <- plot( x_cmpre, breaks = breaks, plot = F)
# ymax <- round_D( max(c(plt1$density, plt2$density)) )
# plot( x_base, 
# col = rgb(0.8, 0.8, 0.8, 0.5), 
# xlim=c(min(x_ttl),max(x_ttl)), 
# ylim = c(0,ymax), 
# cex.axis = 0.8,
# prob = T)
# plot( x_cmpre, breaks = breaks, col = rgb(1, 0, 0, 0.2), add = T, prob = T )
# } else if( (x_base_rg/n_base) < (x_cmpre_rg/n_cmpre) ) {
# plt1 <- plot( x_cmpre , plot = F)
# breaks <- x_base_rg/(x_cmpre_rg/length(plt1$breaks))
# breaks <- ifelse(breaks>1, round(breaks), 1)
# plt2 <- plot( x_base, breaks = breaks, plot = F)
# ymax <- round_D( max(c(plt1$density, plt2$density)) )
# plot( x_base, breaks = breaks, 
# col = rgb(0.8, 0.8, 0.8, 0.5), 
# xlim=c(min(x_ttl),max(x_ttl)), 
# ylim = c(0,ymax), 
# cex.axis = 0.8,
# prob = T)
# plot( x_cmpre, col = rgb(1, 0, 0, 0.2), add = T, prob = T )
# } else	{
# plot( x_base, 
# xlab = 'adfafdasa',
# cex.axis = 0.8,
# col = rgb(0.8, 0.8, 0.8, 0.5), 
# xlim=c(min(x_ttl),max(x_ttl)), 
# prob = T)
# plot( x_cmpre, col = rgb(1, 0, 0, 0.2), add = T, prob = T )
# }
# rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(0, 0, 0, 0.01) )
# grid(5,5)
# box()
# }
# }
# #ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


# #ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# # custom histogram plot (previous not in use)
# #ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# test_hist <- function(x, breaks = 50, min_Nbin = 10, 
# outlier_precentiles = 0.05,
# col = NULL, prob = FALSE, ...)	{

# N <- length(x)
# avg_N <- floor(N/breaks)

# # if there are too many bins
# if( avg_N < min_Nbin )	{
# #at least have an average of 5 or 10 items a break
# breaks <- floor(N/min_Nbin)
# if(breaks<1) breaks <- 1
# avg_N <- floor(N/breaks)
# }

# # if there are lots of bins
# #	[1] represent bottom and top outliers in seperate bars
# #	[2] else; as per normal.
# if ( avg_N/N <= outlier_precentiles  )	{
# x_order <- order(x,decreasing = F)
# x_lrg <- x[head(x_order, avg_N)]
# x_mrg <- x[-c(head(x_order, avg_N),tail(x_order, avg_N))]
# x_urg <- x[tail(x_order, avg_N)]
# x_bar_mid <- table(cut(x_mrg, seq(min(x_mrg),max(x_mrg), length.out = breaks-1)))
# x_bar_all <- c(avg_N, x_bar_mid, avg_N)
# x_names <- as.numeric(gsub("]", "", sapply(lapply(names(x_bar_mid), function(x) strsplit(x,',')), "[[", 1)[2,]))
# x_names <- c(min(x_lrg),max(x_lrg),x_names,max(x_urg))
# x_names <- round(x_names,2)
# if(is.null(col)) col <- c(rgb(0.8, 0.8, 0.8, 0.1),rgb(0.8, 0.8, 0.8, 0.5))[(1:breaks %in% c(1,breaks)) + 1]
# } else	{
# x_bar_all <- table(cut(x, seq(min(x)-1e-15,max(x), length.out = breaks+1)))
# x_names <- as.numeric(gsub("]", "", sapply(lapply(names(x_bar_all), function(x) strsplit(x,',')), "[[", 1)[2,]))
# x_names <- round(x_names,2)
# if(is.null(col)) col <- rgb(0.8, 0.8, 0.8, 0.1)
# }

# barplot(x_bar_all / (1 + prob*sum(x_bar_all)),
# space = 0,
# col = col,
# axisnames = FALSE, ...)
# axis(1, at=seq(0, length(x_bar_all), by=1), labels = FALSE)
# text(seq(0, length(x_bar_all), by=1), par("usr")[3] - 30, labels = x_names, cex = 0.8, srt = 60, xpd = TRUE)
# }
# #ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo



#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# customized hist plot function for data and cluster comparison
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
hist_mdmm <- function(x, y = NULL, 
                      ylim = NULL,
                      breaks = 50, min_Nbin = 10, 
                      outlier_show = TRUE,
                      outlier_q = NULL, 
                      x_col = NULL, y_col = NULL,
                      nlabs = 5,
                      prob = TRUE, ...)	{
  
  N <- length(x)
  avg_N <- floor(N/breaks)
  
  # if there are too many bins
  if( avg_N < min_Nbin )	{
    #at least have an average of min_Nbin items a break
    breaks <- floor(N/min_Nbin)
    if(breaks<1) breaks <- 1
    avg_N <- floor(N/breaks)
  }
  
  if( is.null(outlier_q) ) {
    IQR <- quantile(x,0.75)-quantile(x,0.25)
    lbound <- quantile(x,0.25) - 1.5*IQR
    ubound <- quantile(x,0.75) + 1.5*IQR
  } else	{
    lbound <- quantile(x,outlier_q)
    ubound <- quantile(x,1-outlier_q)
  }
  
  
  if ( outlier_show  )	{
    
    #define x, base outlier, and full data bins
    x_mrg <- x[ x > lbound & x <= ubound]
    bin_mid_cut <- seq(min(x_mrg),max(x_mrg), length.out = breaks-1)
    x_bar_mid <- table(cut(x_mrg, bin_mid_cut))
    x_lrg <- x[ x <= lbound ]
    x_urg <- x[	x > ubound ]
    x_bar_all <- c(length(x_lrg), x_bar_mid, length(x_lrg))
    
    #define y, comparison outlier, and full data bins
    if(!is.null(y))	{
      y_mrg <- y[ y >= lbound & y <= ubound ]
      y_bar_mid <- table(cut(y_mrg, bin_mid_cut))
      y_bar_lrg <- sum(y < lbound)
      y_bar_urg <- sum(y > ubound)
      y_bar_all <- c(y_bar_lrg, y_bar_mid, y_bar_urg)
    }
    
    #histogram bar colours and horizontal label definition
    if(is.null(x_col)) x_col <- c(rgb(0.8, 0.8, 0.8, 0.1),rgb(0.8, 0.8, 0.8, 0.5))[(1:breaks %in% c(1,breaks)) + 1]
    if(is.null(y_col) & !is.null(y)) y_col <- c(rgb(1, 0, 0, 0.1),rgb(1, 0, 0, 0.5))[(1:breaks %in% c(1,breaks)) + 1]
    horizontal_labs <- as.numeric(gsub("]", "", sapply(lapply(names(x_bar_mid), function(x) strsplit(x,',')), "[[", 1)[2,]))
    if( length(x_lrg) != 0) llabs <- c(min(x_lrg),max(x_lrg)) 	else llabs <- c(min(x_mrg),min(x_mrg))
    if( length(x_urg) != 0) ulabs <- max(x_urg)					else ulabs <- 0 
    horizontal_labs <- c(llabs, horizontal_labs, ulabs)
    horizontal_labs[2] <- paste0( '[', horizontal_labs[1], ', ', horizontal_labs[2], ']')
    horizontal_labs <- horizontal_labs[-1]
    horizontal_labs[breaks] <- paste0( '[', horizontal_labs[breaks-1], ', ', horizontal_labs[breaks], ']')
    horizontal_labs[-c(1,breaks)] <- signif(as.numeric(horizontal_labs[-c(1,breaks)]) ,3)
    horizontal_labs[!1:(breaks) %in% c(1, round(seq(1,breaks,length.out = nlabs))[-c(1,nlabs)], breaks)] <- ""
    
  } else	{
    
    #define x, base data bins, and leave out outliers
    x_mrg <- x[ x >= lbound & x <= ubound]
    bin_mid_cut <- seq(min(x_mrg),max(x_mrg), length.out = breaks+1)
    x_bar_all <- table(cut(x_mrg, bin_mid_cut))
    
    #define y, comparison data bins, and leave out outliers
    if(!is.null(y))	{
      y_mrg <- y[ y >= lbound & y <= ubound ]
      y_bar_all <- table(cut(y_mrg, bin_mid_cut))
    }
    
    #histogram bar colours and horizontal label definition
    if(is.null(x_col)) x_col <- rgb(0.8, 0.8, 0.8, 0.1)
    if(is.null(y_col) & !is.null(y)) y_col <- rgb(1, 0, 0, 0.1)
    horizontal_labs <- as.numeric(gsub("]", "", sapply(lapply(names(x_bar_all), function(x) strsplit(x,',')), "[[", 1)[2,]))
    horizontal_labs[1] <- paste0( '[', horizontal_labs[1], ', ', horizontal_labs[2], ']')
    horizontal_labs[breaks] <- paste0( '[', horizontal_labs[breaks-1], ', ', horizontal_labs[breaks], ']')
    horizontal_labs[-c(1,breaks)] <- signif(as.numeric(horizontal_labs[-c(1,breaks)]) ,3)
    horizontal_labs[!1:(breaks) %in% c(1, round(seq(1,breaks,length.out = nlabs))[-c(1,nlabs)], breaks)] <- ""
    
  }
  
  x_plot <- x_bar_all / (1 + prob*sum(x_bar_all))
  if(!is.null(y))	{
    y_plot <- y_bar_all / (1 + prob*sum(y_bar_all))
    if(is.null(ylim)) ylim <- c(0, round_D(max(c(x_plot,y_plot))))
  } else	{
    if(is.null(ylim)) ylim <- c(0, round_D(max(x_plot)))
  }
  
  #--- prep plot background
  barplot(x_plot,
          ylim = ylim,
          space = 0,
          col = rgb(0, 0, 0,0),
          border = rgb(0, 0, 0, 0),
          yaxt = "n", xaxt = "n",	axisnames = FALSE	)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(0, 0, 0, 0.05) )
  barplot(x_plot,
          ylim = ylim,
          space = 0,
          col = rgb(0, 0, 0,0),
          border = rgb(0, 0, 0, 0),
          panel.first = grid(5,5, col = rgb(1, 1, 1, 1),lty = "solid"),
          yaxt = "n", xaxt = "n",	axisnames = FALSE, add = TRUE)
  #------------------------
  xhist <- barplot(x_plot,
                   ylim = ylim,
                   space = 0,
                   col = x_col,
                   yaxt = "n", xaxt = "n",	axisnames = FALSE, add = TRUE	)
  if(!is.null(y))	{
    barplot(y_plot,
            space = 0,
            col = y_col,
            yaxt = "n", xaxt = "n", axisnames = FALSE, add = TRUE, ...)
  }
  axis(2,cex.axis=0.8)
  axis(1, at=seq(0, breaks, by=1)[which(!horizontal_labs %in% "")], labels = FALSE)
  text(xhist, rep(par("usr")[3] - par("usr")[3]/3,length(xhist)),
       labels = horizontal_labs, adj = c(0.5,2), cex = 0.8, xpd = TRUE)
}
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo



#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# customized bar plot function for data and cluster comparison
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
barplot_mdmm <- function(x, y = NULL,
                         ylim = NULL,
                         labels = NULL,
                         x_col = NULL, y_col = NULL,
                         prob = TRUE, ...)	{
  
  if(is.factor(x))	{
    x_plot <- table(x)
  }	else	{
    x_plot <- apply(x,2,sum)
  }
  x_plot <- x_plot / (1 + prob*sum(x_plot))
  if(is.null(x_col)) x_col <- rgb(0, 0, 0.8, 0.2)
  
  if (!is.null(y))	{
    if(is.factor(y))	{
      y_plot <- table(y)
    } else	{
      y_plot <- apply(y,2,sum)
    }
    y_plot <- y_plot / (1 + prob*sum(y_plot))
    if(is.null(ylim)) ylim <- c(0, round_D(max(c(x_plot,y_plot))))
    if(is.null(y_col)) y_col <- rgb(1, 0, 0, 0.4)
  } else	{
    if(is.null(ylim)) ylim <- c(0, round_D(max(x_plot)))
  }
  
  #--- prep plot background
  barplot(x_plot,
          ylim = ylim,
          space = 0,,
          col = rgb(0, 0, 0, 0),
          border = rgb(0, 0, 0, 1),
          yaxt = "n", xaxt = "n", axisnames = FALSE )
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(0, 0, 0, 0.05) )
  barplot(x_plot,
          ylim = ylim,
          space = 0,
          col = rgb(0, 0, 0, 0),
          border = rgb(0, 0, 0, 0),
          panel.first = grid(5,5, col = rgb(1, 1, 1, 1),lty = "solid"),
          yaxt = "n", xaxt = "n", axisnames = FALSE, add = TRUE)
  #------------------------
  xbp <- barplot( x_plot,
                  ylim = ylim,
                  space = 0,
                  col = x_col,
                  yaxt = "n", xaxt = "n", axisnames = FALSE, add = TRUE )
  if(!is.null(y))	{
    barplot(y_plot,
            space = 0,
            col = y_col,
            yaxt = "n", xaxt = "n", axisnames = FALSE, add = TRUE, ...)
  }
  if(is.null(labels)) labels <- names(x_plot)
  axis(2,cex.axis=0.8)
  text(xbp, rep(par("usr")[3] - par("usr")[3]/3,length(xbp)),
       labels = names(x_plot), adj = c(1.05,0.5), cex = 0.8, srt = 60, xpd = TRUE)
}
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

