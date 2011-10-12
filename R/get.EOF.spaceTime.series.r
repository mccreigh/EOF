## this gives all the modes by default.
get.spectralEOF.spaceTime.series <- function( seof, nprocessors=1,
                                               modes=1:max(unlist(lapply( (lapply(seof@band, a.modes)), max ))) )
  {
    get.spaceTime <- function( seof.band, modes=1:(min(c(max(modes),max(seof.band@modes)))) ) {
      print(modes)   
      ntime <- length(seof.band@timeseries[,1])
      nspace <- length(seof.band@EOF[,1])
      perm.space <- if (length(modes)>1) c(1,3,2) else c(1,2)
      perm.time <- if (length(modes)>1) c(3,1,2) else c(2,1)
      print( perm.space)
      print(nspace)
      print(ntime)
      
      resize.space <- function(space) aperm( array( space, c( dim(as.array(space)), ntime) ), perm.space )
      resize.time <- function(time) aperm( array( time, c( dim(as.array(time)), nspace) ), perm.time )
      ( resize.space(Re(seof.band@EOF[,modes])) * resize.time(cos(Arg(seof.band@timeseries[,modes]))) ) +
        ( resize.space(Im(seof.band@EOF[,modes])) * resize.time(sin(Arg(seof.band@timeseries[,modes]))) )
    }
      
    choose.lapply( seof@band, get.spaceTime, nprocessors=nprocessors, modes=modes )
  }
