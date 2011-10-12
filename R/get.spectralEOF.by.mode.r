
## get the space an timeseries information out of a
## spectralEOF object by mode.

gather.real.space <-  function(b, mode=mode) Re( a.EOF(b)[,mode] )
gather.imaginary.space <-  function(b, mode=mode) Im( a.EOF(b)[,mode] )
gather.amplitude.space <- function(b, mode=mode) Mod( a.EOF(b)[,mode] )
gather.phase.space <- function(b, mode=mode) Arg( a.EOF(b)[,mode] )

gather.real.ts <-  function(b, mode=mode) Re( a.timeseries(b)[,mode] )
gather.imaginary.ts <-  function(b, mode=mode) Im( a.timeseries(b)[,mode] )
gather.amplitude.ts <- function(b, mode=mode) Mod( a.timeseries(b)[,mode] )
gather.phase.ts <- function(b, mode=mode) Arg( a.timeseries(b)[,mode] )

gather.nmodes <-  function(b) length(a.modes(b))

get.spectralEOF.by.mode <- function( seof, mode=mode ) {

  if (missing(mode)) {
    warning('A mode must be specified to get.spectralEOF.')
    return(invisible(NULL))
  }

  mode.max=min(unlist((lapply( a.band(seof), gather.nmodes ))))
  if (mode>mode.max) {
    warning(paste("At least one band does not have this many (",mode,") modes.",sep=''))
    return(invisible(NULL))
  }

  if (seof@complex) {
    invisible( list(
                    space= list(
                      real=lapply( a.band(seof), gather.real.space, mode=mode),
                      imag=lapply( a.band(seof), gather.imaginary.space, mode=mode),
                      amp=lapply( a.band(seof), gather.amplitude.space, mode=mode),
                      phase=lapply( a.band(seof), gather.phase.space , mode=mode),
                      lon=a.lon(seof),
                      lat=a.lat(seof)
                      ),

                    time= list(
                      real=lapply( a.band(seof), gather.real.ts, mode=mode),
                      imag=lapply( a.band(seof), gather.imaginary.ts, mode=mode),
                      amp=lapply( a.band(seof), gather.amplitude.ts, mode=mode),
                      phase=lapply( a.band(seof), gather.phase.ts , mode=mode),
                      POSIXct=a.POSIXct(seof)
                      )
                    )
              )
  } else {
    invisible( list(
                    space= list(
                      real=lapply( a.band(seof), gather.real.space, mode=mode ),
                      lon=a.lon(seof),
                      lat=a.lat(seof)
                      ),

                    time= list(
                      real=lapply( a.band(seof), gather.real.ts, mode=mode),
                      POSIXct=a.POSIXct(seof)
                      )
                    )
              )

  }
  
}
    
