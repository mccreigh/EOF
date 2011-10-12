######################################################################
## method: ggEOF
## takes an EOF object and returns a list (of ggplot objects)
## use ggplot to plot an EOF obj/class.


if (!isGeneric("ggEOF")) {  ## creates a generic function
  fun <- if (is.function("ggEOF")) EOF else function( eof, ...) standardGeneric("ggEOF")
  setGeneric("ggEOF", fun)
}

setMethod("ggEOF", c("EOF"),
          function( eof,  nmodes=4, ncol.space=ceiling(sqrt(nmodes)), worldmap=FALSE, cut.space=NULL, ...) {

            require(ggplot2)
            
            ## stopping criteria is whack for this example, bc so much temporal autocorrelation in pdsi?
            ##nmodes=max(which( (eof$W[1:(neof-1)]- eof$W[2:neof])/eof$dW[1:(neof-1)]  > .5 ))-1
            
            ## what if i just want to specify an individual mode which is not the first one??
            
            nmax.eof=length(a.modes(eof))
            nmodes=min(nmodes, nmax.eof)

            ##############################
            ## space
            ## build spatial data frame
            nspace=length(a.lon(eof))
            space <- data.frame( lon=a.lon(eof), lat=a.lat(eof) )
            names.pre <- names(space)
            for (mm in 1:nmodes) space <- cbind(space,a.EOF(eof)[,mm])
            names(space) <- c(names.pre, paste('Mode.',1:nmodes,sep=''))
            
            space <- melt(space, id.vars=c("lon","lat"))
            if (!is.null(cut.space)) {
              breaks=c(-1)* max(abs(space$value)) + ((0:cut.space) *diff( c(-1,1)*max(abs(space$value)))/cut.space)
              space$value=cut(space$value, breaks, include.lowest=TRUE)
            }
            names(space) <- c('lon','lat','modes','values.space')

            
            ## ggplot object in space
            gg.space <- ggplot( space )
            ## titles / facet titles 
            ggspace <- if (nmodes>1) {
              gg.space <- gg.space + facet_wrap(~ modes, ncol=ncol.space)
            } else {
              gg.space <- gg.space + opts(title='Mode 1')
            }

            ## include a world map outline?
            if (worldmap) gg.space <- gg.space + world.map()

            ## example
            ## gg.space + geom_tile( aes(x=lon,y=lat,fill=PDSI) ) + coord_map( orientation=c(90,0,0) ) + ylim(-60, 77)

            ##############################
            ## variance explaned / eigen values.
            variance <- data.frame( modes=a.modes(eof),
                                    variance.percent=a.variance.frac(eof) )
            gg.var <- ggplot( variance ) 
            
            ##############################
            ## time            
            time=data.frame( time=a.POSIXct(eof) )
            for (mm in 1:nmodes) time <- cbind(time,a.timeseries(eof)[,mm])
            names(time) <- c( 'time', paste("Mode.",1:nmodes,sep='') )

            time <- melt(time, id.vars=c('time'))
            names(time) <- c('time','modes','values.time')

            gg.time <- ggplot( time )
            ## ggplot( time, aes(x=time, y=PDSI, colour=modes) ) + geom_line() + facet_wrap(~modes)
          

            ## invisible( new( 'ggEOF', space=gg.space, time=gg.tim, var=gg.var ) )
            invisible( list( space=gg.space, time=gg.time, var=gg.var ) )
            
          })
            
