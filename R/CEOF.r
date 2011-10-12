######################################################################
## class: CEOF
## i dont know how to do partial inheritance, so this is may be a bit redundant
if (!isClassUnion("corrCovarOrNULL")) setClassUnion("corrCovarOrNULL", c("corrCovar","NULL") )
setClass("CEOF", representation(CEOF="array",
                                timeseries="array",
                                modes="vector",
                                lon="vector", lat="vector",
                                POSIXct="POSIXct",
                                eigen="vector",
                                eigen.uncert="vector",
                                variance.frac= "numeric",
                                variable="character",
                                corr.covar="character",
                                fast= "logical",
                                mtx="corrCovarOrNULL"
                                ) )

set.accessors("CEOF")

######################################################################
## method: CEOF
## Complex EOF 

if (!isGeneric("CEOF")) {  ## creates a generic function
  fun <- if (is.function("CEOF")) EOF else function( mtx, ...) standardGeneric("CEOF")
  setGeneric("CEOF", fun)
}

setMethod("CEOF", c("corrCovar"),
          function( mtx, ... ) {
            
            if (a.complex(mtx)==FALSE) warning("CEOF failure: input matrix must be complex.", immediate.=TRUE)

            eof <- EOF( mtx, ...)  ## same thing, really.... 

            invisible( new("CEOF",
                           CEOF=eof@EOF,
                           modes=eof@modes,
                           lon=eof@lon,
                           lat=eof@lat,
                           timeseries=eof@timeseries,
                           POSIXct=eof@POSIXct,
                           eigen=eof@eigen,
                           eigen.uncert=eof@eigen.uncert, 
                           variance.frac=eof@variance.frac,
                           variable=eof@variable,
                           corr.covar=eof@corr.covar,
                           fast=eof@fast,
                           mtx=eof@mtx ) )            

            ## should probably include phase and amplitude sequences
            ## should I separate out the Re and Im spatial components?
            
          }
          )

