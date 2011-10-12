######################################################################
## class: spectralEOF
## i dont know how to do partial inheritance, so this is may be a bit redundant
## with in the bands slot, a list contains spectralEOFband objects
setClass("spectralEOFband", representation(EOF="array",
                                            timeseries="array",
                                            modes="vector",
                                            eigen="vector",
                                            eigen.uncert="vector",
                                            variance.frac= "numeric"#,
#                                            mtx="corrCovar"
                                            ) )
set.accessors("spectralEOFband")

setClass("spectralEOF", representation(band="list",
                                lon="vector", lat="vector",
                                POSIXct="POSIXct",
                                variable="character",
                                bands="matrix",
                                corr.covar="character",
                                complex='logical',
                                fast= "logical"
                                ) )
set.accessors("spectralEOF")


if (!isGeneric("spectralEOF")) {  ## creates a generic function
  fun <- if (is.function("spectralEOF")) EOF else function( stb, ...) standardGeneric("spectralEOF")
  setGeneric("spectralEOF", fun)
}


#####################################################################
setMethod("spectralEOF", c("spaceTimeBands"),
function( stb, corr=!covar, covar=!corr, complex=FALSE, nprocessors=1, ... ) {
  
  if( missing(corr) & missing(covar) ) corr=TRUE
  nbands <- dim(a.bands( stb ))[1]
  
  EOF.band <- function( band ) {
    st <- new("spaceTime", data=band,
                    lon=a.lon(stb), lat=a.lat(stb), POSIXct=a.POSIXct(stb),
                    data.name=a.data.name(stb), data.units=a.data.units(stb) )
    mtx <- if (corr) corrMtx(st, complex=complex) else covarMtx(st, complex=complex)
    eof <- if (!complex) EOF( mtx ) else EOF( mtx )
    return( new("spectralEOFband", EOF=a.EOF(eof), timeseries=a.timeseries(eof),
                  modes=a.modes(eof), eigen=a.eigen(eof), eigen.uncert=a.eigen.uncert(eof),
                  variance.frac=a.variance.frac(eof)#, mtx=mtx
                 ) )
  }    
  
  band.list <- choose.lapply( stb@data, EOF.band, nprocessors=nprocessors, ...)

  corr.covar <- if (corr) 'corr' else 'covar'
  
  invisible( new( "spectralEOF", band=band.list,
              lon=a.lon(stb), lat=a.lat(stb),
              POSIXct=a.POSIXct(stb),
              variable=a.data.name(stb),
              bands=a.bands(stb),
              corr.covar=corr.covar,
              complex=complex,
              fast=FALSE
              ) )
  
}
)
