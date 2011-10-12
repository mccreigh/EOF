## returns a list of ggplot objects for each band in the spectralEOF object
## all arguments can also be passed to ggEOF including plot=c(mode.i,mode.k)

ggSpectralEOF <-  function( seof, modes=1,
                            plot=TRUE, plot.world.map=TRUE,
                            phase.mult=3, do.unwrap.phase=TRUE,
                            ... )
{
  ggseof <- list()
  for (mm in modes) {

    ## get the data for this mode
    mode.data <- get.spectralEOF.by.mode( seof, mode=mm )

    #####################################################################
    ## space
    ## pull out the real spatial part
    space.real <- melt(cbind(as.data.frame(mode.data$space$real),
                             lon=mode.data$space$lon,
                             lat=mode.data$space$lat),
                       id.vars=c("lon","lat") )
    space.real$variable <- factor( space.real$variable, labels=names(seof@band) )

    ## gg the real part
    gg.real.out <- ggplot( space.real )
    gg.real <- gg.real.out +
                 geom_tile( aes(x=lon,y=lat, fill=value) ) +
                 geom_contour( aes(x=lon, y=lat, z=value), color='grey' ) +
                 facet_wrap(~variable, ncol=1) +
                 scale_x_continuous(name='') + scale_y_continuous(name='') +
                 scale_fill_gradient2( low='blue',high='orange',midpoint=0 ) +
#                 xlim( min(space.real$lon), max(space.real$lon) ) +
#                 ylim( min(space.real$lat), max(space.real$lat) ) +
                 opts( plot.margin = unit(c(.5, .5, .25, .25), "lines") ) +
                 opts( legend.key.width = unit(.5, "lines") )
    gg.real <- if (seof@complex) gg.real + opts( title=paste('Real - Mode',mm) ) else
                                 gg.real + opts( title=paste('Mode',mm) )     
    if (plot.world.map) gg.real <-
      gg.real + world.map( alpha=.6, size=.3, lon=space.real$lon,
                          lat=space.real$lat, buffer=10, posi=all(space.real$lon>=0) )
    
    if (seof@complex) {
      ## pull out the imaginary spatial part
      space.imag <- melt(cbind(as.data.frame(mode.data$space$imag),
                               lon=mode.data$space$lon,
                               lat=mode.data$space$lat),
                         id.vars=c("lon","lat") )
      space.imag$variable <- factor( space.imag$variable, labels=names(seof@band) )
      
      ## pull out the phase components and amplitude
      space.phase.amp <- melt(cbind(as.data.frame(mode.data$space$amp),
                                    lon=mode.data$space$lon,
                                    lat=mode.data$space$lat),
                              id.vars=c("lon","lat") )    
      space.phase.amp$variable <- factor( space.amp$variable, labels=names(seof@band) )
      space.phase.amp$amp <- space.phase.amp$value
      space.phase.amp$phase.re <- space.real$value
      space.phase.amp$phase.im <- space.imag$value
      space.phase.amp <- space.phase.amp[-which(names(space.phase.amp)=='value')]

      ## gg the imaginary part
      gg.imag.out <- ggplot( space.imag )
      gg.imag <- gg.imag.out +
                   geom_tile( aes(x=lon,y=lat, fill=value) ) +
                   geom_contour( aes(x=lon, y=lat, z=value), color='grey' ) +
                   facet_wrap(~variable, ncol=1) +
                   scale_x_continuous(name='') + scale_y_continuous(name='') +
                   xlim( min(space.real$lon), max(space.real$lon) ) +
                   ylim( min(space.real$lat), max(space.real$lat) ) +
                   opts( title=paste('Imaginary - Mode',mm) ) +
                   opts( plot.margin = unit(c(.5, .5, .25, .25), "lines") ) +
                   opts( legend.key.width = unit(.5, "lines") )    
      if (plot.world.map) gg.imag <- gg.imag + world.map( alpha=.6, size=.3 )

      ## gg the amp and phase.
      gg.phase.amp.out <- ggplot( space.phase.amp )
      gg.phase.amp <- gg.phase.amp.out +
                        geom_tile(aes(x=lon, y=lat, fill=amp)) +
                        geom_contour(aes(x=lon, y=lat, z=amp), color='grey', alpha=.4) +
                        geom_segment( aes_string( x="lon", y="lat",
                                                 xend=paste("lon+",phase.mult,"*phase.im"),
                                                 yend=paste("lat+",phase.mult,"*phase.re" ),
                                                 arrow="arrow(length=unit(.1,'cm'))" ),
                                     color='white' ) +
                        facet_wrap(~variable, ncol=1) +
                        scale_x_continuous(name='') + scale_y_continuous(name='') +
                        xlim( min(space.real$lon), max(space.real$lon) ) +
                        ylim( min(space.real$lat), max(space.real$lat) ) +
                        opts( title=paste('Amp & Phase - Mode',mm) ) +
                        opts( plot.margin = unit(c(.5, .5, .25, .25), "lines") ) +
                        opts( legend.key.width = unit(.5, "lines") )
      if (plot.world.map) gg.phase.amp <- gg.phase.amp + world.map( alpha=.6, size=.3 )
    }

    #####################################################################
    ## the timeseries parts
    nbands=length(names(seof@band))
    unwrap.yn <- function(p,y) if (y) unwrap.phase(p) else p
    in.list <- as.list(1:nbands); names(in.list)=names(seof@band)
    get.time.funct <- if (seof@complex) {
      function(bb) melt( data.frame(real=mode.data$time$real[[bb]],
                                    imag=mode.data$time$imag[[bb]],
                                    phase=unwrap.yn(mode.data$time$phase[[bb]],
                                      do.unwrap.phase),
                                    amp=mode.data$time$amp[[bb]],
                                    POSIXct=mode.data$time$POSIXct), id.vars='POSIXct') 
    } else {
      function(bb) melt(data.frame(time=mode.data$time$real[[bb]],
                                   POSIXct=mode.data$time$POSIXct), id='POSIXct')
    }
    
    time.bands <- lapply( in.list, get.time.funct )
    gg.time.out <- lapply( time.bands, ggplot )

    make.gg.time <- function(gg)
      { gg$data + geom_line( aes(x=POSIXct, y=value) ) + facet_grid(variable~. , scales='free_y') +
          scale_x_datetime(name='', major='20 years') +
          scale_y_continuous(name='') +
          opts(title=gg$name,
               plot.title = theme_text(size = 6 * 1.2),
               plot.margin = unit(c(-.1, -.1, -.5, -.5), "lines") )
      }

    combine.list.w.name <- function( l ) {      
      out <- lapply( as.list(1:length(names(l))), function(ii) list( data=l[[ii]], name=names(l)[ii] ) )
      names(out) <- names(l)
      out
    }
    gg.time.name.combo <- combine.list.w.name(gg.time.out)
    
    gg.time <- lapply( gg.time.name.combo, make.gg.time )

    
    if (plot) {
      ## save default plots to output
      gg.real.out <- gg.real
      if (seof@complex) {
        gg.imag.out <- gg.imag
        gg.phase.amp.out <- gg.phase.amp
      }
      gg.time.out <- gg.time
      
      ## set up view ports
      vplayout <- function(x,y) viewport(layout.pos.row=x,layout.pos.col=y)
      check.new.dev()
      grid.newpage()

      pushViewport(viewport(layout=grid.layout(nbands, if (seof@complex) 4 else 2 )))
      
      ## plot
      print(gg.real, vp=vplayout(1:nbands,1))
      if (seof@complex) {
        print(gg.imag, vp=vplayout(1:nbands,2))
        print(gg.phase.amp, vp=vplayout(1:nbands,3))
      }
      for (bb in 1:nbands) print(gg.time[[bb]], vp=vplayout(bb, if (seof@complex) 4 else 2))
      
    }

    ## arrange data for output
    ggseof[[mm]] <- if (seof@complex) {
      list( space = list( real=gg.real.out, imag=gg.imag.out, phase.amp=gg.phase.amp.out ),
           time = gg.time.out )
    } else {
      list( space = gg.real.out, time = gg.time.out )
    }
    
  }

  invisible(ggseof)

}
