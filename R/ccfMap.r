## this function computes the ccf (using ccf() in R core package stats) of a timeseries
## against all the timeseries in a spatial field (a spaceTime object). It is assumed that the
## spaceTime object contains no NAs. The returned quantities are spatial maps of
## maximum/minimum correlation,  the lags of maximum/minimum correlation and the significances of each. 
## Since this is a massively parrallell computation, mclapply (from the multicore package)
## is invoked whenever nprocessors>1.
## 

## I assume that ccf uses pearson correlation, though it's not specific.
## Significances are computed using pearson under corr.test() under the lags found by ccf.
## TODO: I should probably compute the correlation under corr.test() too and allow the method to be specified. 

## Class definition of a correlation map
setClass("ccfMap", 
         representation( max.corr='vector',
                         min.corr='vector',
                         lag.max='vector',
                         lag.min='vector',
                         p.max='vector',
                         p.min='vector',
                         is.signif.max='vector',
                         is.signif.min='vector',
                         confidence.level='numeric',
                         lon='vector',
                         lat='vector',
                         POSIXct='POSIXct',
                         time.range='character',
                         timeSeries.name='character',
                         spaceTime.name='character',
                         type='character'
                        ) 
         )
set.accessors( 'ccfMap' )

## Generic
if (!isGeneric("ccfMap")) {  ## creates a generic function
  fun <- if (is.function("ccfMap")) ccfMap else function( timeSeries, spaceTime, ...) standardGeneric("ccfMap")
  setGeneric("ccfMap", fun)
}

setMethod("ccfMap", c("timeSeries", "spaceTime"),
          function( timeSeries, spaceTime, confidence.level=.95, lag.max=12, nprocessors=1, ... )
{  

  if (nprocessors >1) {
    have.multicore <- require(multicore)      
    if (!have.multicore) {
      print(paste("multicore package is NOT present, will run",nsimulations,"simulations on 1 processor."))
      print("Do you wish to continue? (y/n)")
      if (tolower(substr(readLines(n=1),1,1))=='n') return(NULL) else nprocessors <- 1
    }
  }    

  wh.ts <- which( a.POSIXct(timeSeries) %in% a.POSIXct(spaceTime) )
  wh.st <- which( a.POSIXct(spaceTime) %in% a.POSIXct(timeSeries) )
  
  a.data(timeSeries) <- a.data(timeSeries)[wh.ts]
  a.POSIXct(timeSeries) <- a.POSIXct(timeSeries)[wh.ts]
  
  spaceTime <- subset( spaceTime, wh.st, dim='time' )

  njoint.time=length(wh.ts)
  
  cor.func <- function (st.time) {

    get.p <- function(lag) { ## must be defined here to find st.time[
      p <- if (lag>0) cor.test( a.data(timeSeries)[(1+lag):njoint.time], st.time[1:(njoint.time-lag)] )$p.value else
      cor.test( a.data(timeSeries)[1:(njoint.time+lag)], st.time[(1-lag):njoint.time] )$p.value
      p      
    }
    
    ccf <- ccf( a.data(timeSeries), st.time, lag.max=lag.max, plot=FALSE, ...)
    if (any(is.na(ccf$acf))) return( list( max=NA, lag.max=NA, p.max=NA, min=NA, lag.min=NA, p.min=NA ) )
    wh.max=which(ccf$acf == max(ccf$acf)); wh.min=which(ccf$acf == min(ccf$acf))    
    max <- ccf$acf[wh.max]     ; min <- ccf$acf[wh.min]
    lag.max <- ccf$lag[wh.max] ; lag.min <- ccf$lag[wh.min]
    p.max <- get.p(lag.max)    ; p.min <- get.p(lag.min)
    list( max=max, lag.max=lag.max, p.max=p.max,
         min=min, lag.min=lag.min, p.min=p.min )
  }

  ## have to get SpaceTime data into a list where each row is an entry.
  result <- if (nprocessors==1) lapply( as.list(as.data.frame(t(a.data(spaceTime)))), cor.func ) else {
    ## mc.preschedule is faster for shorter calls. mc.set.seed=TRUE is **essential**, otherwise the same seed is used on all cores.
    mclapply( as.list(as.data.frame(t(a.data(spaceTime)))), cor.func, mc.set.seed=TRUE, mc.cores=nprocessors, mc.preschedule=TRUE )
  }

  ## build a ccfMap  object
  p.max <- as.numeric(sapply( result, '[[', 'p.max' ))  ## used in calculating is.signif.max/min
  p.min <- as.numeric(sapply( result, '[[', 'p.min' ))
  
  return( new('ccfMap', 
              max.corr=as.numeric(sapply( result, '[[', 'max' )),
              min.corr=as.numeric(sapply( result, '[[', 'min' )),
              lag.max=as.numeric(sapply( result, '[[', 'lag.max' )),
              lag.min=as.numeric(sapply( result, '[[', 'lag.min' )),
              p.max=p.max,
              p.min=p.min,
              is.signif.max= p.max <= (1-confidence.level),
              is.signif.min= p.min <= (1-confidence.level),
              confidence.level=confidence.level,
              lon=a.lon(spaceTime),
              lat=a.lat(spaceTime),
              POSIXct=a.POSIXct(spaceTime),
              time.range=paste(format(a.POSIXct(spaceTime)[c(1,length(a.POSIXct(spaceTime)))],'%d/%m/%Y'),collapse='-'),
              timeSeries.name=a.data.name(timeSeries),
              spaceTime.name=a.data.name(spaceTime)
              )
         )
  
}
)

## a method for subsetting a correlation map
subset.ccfMap <- function( ccfMap, inds ) {

  a.max.corr(ccfMap) <- a.max.corr(ccfMap[inds])
  a.min.corr(ccfMap) <- a.min.corr(ccfMap[inds])
  a.lag.max(ccfMap) <- a.lag.max(ccfMap[inds])      
  a.lag.min(ccfMap) <- a.lag.min(ccfMap[inds])
  a.p.max(ccfMap) <- a.p.max(ccfMap[inds])
  a.p.min(ccfMap) <- a.p.min(ccfMap[inds])
  a.is.signif.max(ccfMap) <- a.is.signif.max(ccfMap[inds])
  a.is.signif.min(ccfMap) <- a.is.signif.min(ccfMap[inds])
  a.lon(ccfMap) <- a.lon(ccfMap[inds])
  a.lat(ccfMap) <- a.lat(ccfMap[inds]) 

  ccfMap

}
