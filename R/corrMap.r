

## Class definition of a correlation map
setClass("corrMap", 
         representation( p='vector',
                         est='vector',
                         lon='vector',
                         lat='vector',
                         POSIXct='POSIXct',
                         time.range='character',
                         is.signif='vector',
                         confidence.level='numeric',
                         timeSeries.name='character',
                         spaceTime.name='character',
                         method='character',
                         alternative='character',
                         parameter='numericNULL'
                        ) 
         )
set.accessors( 'corrMap' )

## Generic
if (!isGeneric("corrMap")) {  ## creates a generic function
  fun <- if (is.function("corrMap")) EOF else function( timeSeries, spaceTime, ...) standardGeneric("corrMap")
  setGeneric("corrMap", fun)
}

setMethod("corrMap", c("timeSeries", "spaceTime"),
         
corrMap <- function( timeSeries, spaceTime, confidence.level=.95, ... )
{  
  wh.ts <- which( a.POSIXct(timeSeries) %in% a.POSIXct(spaceTime) )
  wh.st <- which( a.POSIXct(spaceTime) %in% a.POSIXct(timeSeries) )
  
  a.data(timeSeries) <- a.data(timeSeries)[wh.ts]
  a.POSIXct(timeSeries) <- a.POSIXct(timeSeries)[wh.ts]
  
  spaceTime <- subset( spaceTime, wh.st, dim='time' )
  
  cor.func <- function (st.time) { cor.test( a.data(timeSeries), st.time, confidence.level=confidence.level, ... ) }
  
  result <- apply( a.data(spaceTime), 1, 'cor.func' )
  
  ## build a corrMap  object
  p <- sapply( result, '[[', 'p.value' ) ## used in calculating is.signif

  return( new('corrMap', 
              p=p,
              est=sapply( result, '[[', 'estimate' ),
              lon=a.lon(spaceTime),
              lat=a.lat(spaceTime),
              POSIXct=a.POSIXct(spaceTime),
              time.range=paste(format(a.POSIXct(spaceTime)[c(1,length(a.POSIXct(spaceTime)))],'%d/%m/%Y'),collapse='-'),
              is.signif= p <= (1-confidence.level),
              confidence.level=confidence.level,
              timeSeries.name=a.data.name(timeSeries),
              spaceTime.name=a.data.name(spaceTime),
              method=result[[1]]$method,
              alternative=result[[1]]$alternative,
              parameter=result[[1]]$parameter
              )
         )
  
}
)


## a method for subsetting a correlation map
#subset.corrMap <-  function( corrMap, inds ) {
#
#}
