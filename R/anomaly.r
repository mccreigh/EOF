######################################################################
## class: anomaly
setClass("anomaly",
         representation(anomaly = "spaceTime",
                        anom.dim="character",
                        period.fmt="character",
                        standardize="logical",
                        allow.missing="logical",
                        wh.keep= "vector",                        
                        input= "spaceTime")
         )

set.accessors('anomaly')

######################################################################
## method: calculate anomalies
## returns an anomaly class from a spaceTime input.
## this method on spaceTime returns an anomaly class object.
## is meant to be an internal, but could be used more generally.
if (!isGeneric("anomaly")) {  ## creates a generic function
  if (is.function("anomaly"))
    fun <- anomaly
  else fun <- function(st, ...) standardGeneric("anomaly")
  setGeneric("anomaly", fun)
}

## might eventually add a climatology of means and sd to return??
setMethod("anomaly",c("spaceTime"),
          function(st, anom.dim='time', period.fmt='n/a',
                   standardize=FALSE, allow.missing=FALSE ) {

            st.in=st  ## include input with the anomaly class output

            nspace <- length(a.lon(st))
            ntime <- length(a.POSIXct(st))

            ## anomalies in space or in time?
            if (anom.dim!='time' & anom.dim != 'space')
              warning("anomaly: anom.dim neither space nor time, using time.")

            ## only allow 1 div in space.
            divs <- if (anom.dim=='time') {
              if (period.fmt!='n/a') as.numeric(format(a.POSIXct(st), period.fmt )) else (1:ntime)*0+1
            } else (1:space)*0+1
            uniq.divs=unique(divs)
            nmax.cycles=max(table(divs))
            
            ## loop over divs in the dimension anom.dim
            for (dd in uniq.divs) {
              div <- which(divs==dd)             
              ncycles=length(div)  ## when anom.dim=='time', ncycles=nspace
              div.arr=if (anom.dim=='time') a.data(st)[,div] else a.data(st)[div,]

              anom.dim.num <- if (anom.dim=='time') 1 else 2              
              ## mean
              mean <- apply( div.arr, anom.dim.num, 'mean', na.rm=allow.missing)
              if (anom.dim=='time') {
                a.data(st)[,div] <- a.data(st)[,div] - array( data=replicate( ncycles, mean ), dim=c(nspace,ncycles) )
              } else {
                a.data(st)[div,] <- a.data(st)[div,] - t(array( data=replicate( ncycles, mean ), dim=c(ntime,ncycles)))
              }               

              ## sdev - note that these function could be customized as they are passed to apply by string.
              if (standardize) {
                sd <- apply( div.arr, anom.dim.num, 'sd', na.rm=allow.missing)
                if (anom.dim=='time') {
                  a.data(st)[,div] <- a.data(st)[,div] / array( data=replicate( ncycles, sd ), dim=c(nspace,ncycles) )
                } else {
                  a.data(st)[div,] <- a.data(st)[div,] / t(array( data=replicate( ncycles, sd ), dim=c(ntime,ncycles)))
                }               
              }
            } #for dd

            ## Remove any spatial points which are constant or all missing in time.
            ## --
            ## if anom.dim==time & allow.missing==FALSE: this will leave only points w complete timeseries.
            ## if anom.dim==time & allow.missing==TRUE: points with some data (incomplete timeseries) are retained.
            ## if anom.dim==space & allow.missing==FALSE: you can have many missing times have no data left.
            ## if anom.dim==space & allow.missing==TRUE: you can have many missing times have no data left.
            ## --
            ## if standardize==TRUE constant points are automatically picked up b/c their sd=0, an Inf (not NA) value.
            keep.space <-
             if (allow.missing) which( rowSums(is.finite(a.data(st))) >0 ) else which( rowSums(is.finite(a.data(st))) == ntime )
##            keep.space <- which( rowSums(is.finite(a.data(st))) == ntime )
            a.data(st) <- a.data(st)[keep.space,]
            a.lon(st) <- a.lon(st)[keep.space]
            a.lat(st) <- a.lat(st)[keep.space]
            if (!standardize & !(anom.dim=='time' & !allow.missing) ) {# if anom.dim='time' & !allow.missing then redundant.
              keep.space2 <- which( apply( a.data(st), 1, 'sd', na.rm=allow.missing) != 0 )  
              if (length(keep.space2) != length(keep.space)) {                               
                a.data(st) <- a.data(st)[keep.space2,]
                a.lon(st) <- a.lon(st)[keep.space2]
                a.lat(st) <- a.lat(st)[keep.space2]
                keep.space <- keep.space[keep.space2]                
              }
            }
                        
            ## warn about times where all values are missing.            
            keep.time=which ( colSums(is.finite(a.data(st))) < 1 )
            if (length(keep.time) > 1 ) warning('anomaly: warning some times have no valid data.')

            ## if missing data are being allowed, these are filled with zero.
            if (allow.missing) a.data(st)[which(!is.finite(a.data(st)))] = 0

            invisible( new("anomaly", anomaly=st,
                           anom.dim=anom.dim, period.fmt=period.fmt,
                           standardize=standardize, allow.missing=allow.missing, 
                           wh.keep=keep.space,
                           input=st.in ) )
            
          }
          )
