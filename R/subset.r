######################################################################
## method: subset on spaceTime
## spatial subsetting: use only one argument (start) of indices
##                     in the space dimension to specify the
##                     elements to keep.
## temporal subsetting: one or two arguments. with one argument,
##                      same as spatial, but in the time dim. if
##                      start and end times supplied in american
##                      mm/dd/yyyy format, keep those times and intervening.


if (!isGeneric("subset")) {  ## creates a generic function
  if (is.function("subset"))
    fun <- subset else fun <- function(x, ...) standardGeneric("subset")
  setGeneric("subset", fun)
}

setMethod("subset", "spaceTime", 
          function(x, start, end, dim='time') {

            if (dim != 'time' & dim !='space') warning('subset: dim is neither time or space, using time.')            
            if (dim=='time') {

              if (!missing('end')) {

                ## temporal range subsetting
                if (is.POSIXct(start)) start <- format(start,'%m/%d/%Y')
                if (is.POSIXct(end)) end <- format(end,'%m/%d/%Y')
                wh.start <- which(format(a.POSIXct(x),'%m/%d/%Y') == start)
                wh.end <- which(format(a.POSIXct(x),'%m/%d/%Y') == end)
                if (length(wh.start)==0 | length(wh.end)==0)
                  warning( 'Either start or end time not found in the timeSeries, check their values.', immediate.=TRUE)
                
                st <- a.data(x)[,wh.start:wh.end]
                if (length(wh.start:wh.end)==1) st <- matrix(st, ncol=1)
                a.data(x) <- st
                a.POSIXct(x) <- a.POSIXct(x)[wh.start:wh.end]
                return(x)

              } else { ## temporal index subsetting

                a.data(x) <- a.data(x)[,start]
                a.POSIXct(x) <- a.POSIXct(x)[start]
              }

            } else { #space index dimension

              a.data(x) <- a.data(x)[start,]
              a.lon(x) <- a.lon(x)[start]
              a.lat(x) <- a.lat(x)[start]
              
            }
            
            invisible(x)
            
          }
          )
