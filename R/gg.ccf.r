
## other.vars and var are both data frames
## vars has the variable of interest, against which other.vars are ccf'd, and some POSIXct abcissa
## other.vars has multiple timeseries against another POSIXct abcissa which they all share. 
## both data frames may contain NAs. Note that no holes are allowed when calculating the ccf of the timeseries
## so the biggest chunk of the timeseries with no holes in common to both vars and each (separately)
## other.vars is used. 

## the return value is a list with 2 parts, ccf and gg. 
## ccf is a list with the names of the other.vars, each contains a $frame with their timeseries and a POSIXct
## that was used in the calculation. 
## gg is the ggplot object formed by default. if you dont like the default you can get out the data frame at its
## core via return.value$gg$data and build a new ggplot object from there.

## function to find the largest continuous part of a timeseries with no NAs.
take.biggest.continuous <- function(ts) {
  wh.na <- which(is.na(ts))  
  n.wh.na <- length(wh.na)
  if (n.wh.na==0) return( list( inds=1:length(ts), ts=ts) )
  extra.length=0
  if (wh.na[1]!=1) { wh.na <- c(0,wh.na); extra.length=extra.length+1}
  if (wh.na[n.wh.na+extra.length]!=length(ts)) { wh.na <- c(wh.na,length(ts)+1); extra.length=extra.length+1}
  sizes <-  wh.na[2:(n.wh.na+extra.length)]-wh.na[1:(n.wh.na+extra.length-1)]-1
  wh.biggest <- which( sizes == max(sizes) )
  inds=(wh.na[wh.biggest]+1):(wh.na[wh.biggest+1]-1)
  invisible( list(  inds=inds , ts=ts[ inds ] ) )
}

gg.ccf <- function( other.vars, var, plot=TRUE, ... ) {

  apply.ccf <- function( ts, ... ) {

    ts.cont <- take.biggest.continuous( ts )
    ts <- data.frame( POSIXct=other.vars$POSIXct[ts.cont$inds], data=ts.cont$ts )    

    wh.ts <- which( ts$POSIXct %in% var$POSIXct )
    wh.var <- which( var$POSIXct %in% ts$POSIXct )
    ts <- ts[wh.ts,]
    var <- var[wh.var,]

    ## sort the data sets to be sure
    ts <- ts[sort(as.numeric(ts$POSIXct),index.return=TRUE)$ix,]
    var <- var[sort(as.numeric(var$POSIXct),index.return=TRUE)$ix,]

    ## now find the largest complet section of with out missing data in each data set.
    wh.POSIXct.ts=which(names(ts)=='POSIXct')
    wh.POSIXct.var=which(names(var)=='POSIXct')
    ccf.out <- ccf( ts[,-wh.POSIXct.ts], var[,-wh.POSIXct.var], plot=FALSE, ...)
    return( list( frame=ts, ccf=ccf.out ) )
  }

  wh.POSIXct.other <- which(names(other.vars)=='POSIXct')
  ccf <- apply( other.vars[-wh.POSIXct.other], 2, apply.ccf, ...)

  ## gather for ggplotting
  gather<- unlist(lapply( ccf, function(z){ return( list( ccf=as.vector(z$ccf$acf), lag=as.vector(z$ccf$lag) ) ) } ) )
  names <- strsplit( names(gather), '\\.' )
  
  variable <- unlist(lapply( names, function(s){ len=length(s); paste(s[1:(len-1)],collapse='.') } ))
  quantity <- unlist(lapply( names, function(s){ len=length(s); substr(s[len],1,3) } ))

  if (!all( variable[which(quantity=='ccf')] == variable[which(quantity=='lag')] )) 
     warning('Somthing has gone wrong, check code.', immediate.=TRUE )
  plot.frame <- data.frame( ccf=gather[which(quantity=='ccf')], lag=gather[which(quantity=='lag')], 
                            variable=variable[which(quantity=='ccf')] )
  breaks=((max(plot.frame$lag)-min(plot.frame$lag)) %/% 5) %/% 2 
  breaks=5*((-1*breaks):breaks)
  plot.frame$lag <- factor(plot.frame$lag)
  print("you can ignore these warnings about ymin:")
  gg.ccf <-ggplot( plot.frame, aes(x=lag,y=ccf) ) + geom_bar( size=.3 ) + facet_wrap( ~variable, ncol=1 ) + 
                   scale_x_discrete(breaks=breaks)

  if (plot) print(gg.ccf)

  invisible( list( ccf=ccf, gg=gg.ccf ) )
 
  
}


summary.gg.ccf <- function( ggCCF ) {

  ## done in 2 steps.
  ## 1) pull the min/max lag.min/lag.max data from each ggCCF$ccf individually,
  ## 2) pull the individual stats from that list and list them in order.

  ## step 1
  get.max.corr.lag <- function( ccf.list ) {
    ccf <- ccf.list$ccf
    wh.max=which(ccf$acf == max(ccf$acf))
    wh.min=which(ccf$acf == min(ccf$acf))
    list( max=ccf$acf[wh.max],
         lag.max=ccf$lag[wh.max],
         min=ccf$acf[wh.min],
         lag.min=ccf$lag[wh.min]
         )
  }

  ## list of min max with lags of each= lmm.ccf
  lmm.ccf <- lapply( ggCCF$ccf, get.max.corr.lag )

  ## step 2  
  max.sort <- sort(unlist(lapply( lmm.ccf, '[[', 'max')), index.return=TRUE)
  max.corr <- rev(max.sort$x)
  max.lag <- unlist(lapply( lmm.ccf, '[[', 'lag.max'))[rev(max.sort$ix)]
  max.frame <- data.frame(name=attributes(max.corr)$names,
                          max.corr=as.numeric(round(100*max.corr)/100),
                          lag.max=as.numeric(max.lag))
  
  min.sort <- sort(unlist(lapply( lmm.ccf, '[[', 'min')), index.return=TRUE)
  min.corr <- min.sort$x
  min.lag <- unlist(lapply( lmm.ccf, '[[', 'lag.min'))[min.sort$ix]
  min.frame <- data.frame(name=attributes(min.corr)$names,
                          min.corr=as.numeric(round(100*min.corr)/100),
                          lag.min=as.numeric(min.lag))
  
  cbind( max.frame, min.frame )  
}
