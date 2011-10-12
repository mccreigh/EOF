## ggplot method for spaceTime objects
## purpose: plot multiple times in a timeseries.

## field: a spaceTime field
## times: a numeric index vector in the time dimension or
##        a character vector (format mm/dd/yyyy) which can
##        be matched against times in the optional POSIXct
##        variable.
## breaks: =10, number of plotting levels
## ncol: =1, number of colums in facet_wrap
## world: include a world map

ggField <- function( st, times, breaks=10, ncol=1, world=TRUE, plot=TRUE, ... ) {

  wh.times=times
  times.char <- as.character(times)  ## facet titles
  if (is.character(times) & exists('POSIXct') ) {
    for (tt in 1:length(wh.times))
      { nparts=length(strsplit(times[tt],'/')[[1]])
        if (nparts==2) {y.ind=2; m.ind=1} else {y.ind=3; m.ind=1; d.ind=2}
        y <- as.numeric(format(POSIXct,'%Y'))==as.numeric(strsplit(times[tt],'/')[[1]][y.ind]) 
        wh.y=which(y)
        ym <- as.numeric(format(POSIXct,'%m'))==as.numeric(strsplit(times[tt],'/')[[1]][m.ind]) & y
        wh.ym <- which(ym)
        if (nparts>2)
          wh.ymd <- which(as.numeric(format(POSIXct,'%d'))==as.numeric(strsplit(times[tt],'/')[[1]][d.ind]) & ym) else wh.ymd=NULL

        print(wh.y)
        print(wh.ym)
        print(wh.ymd)
        
        wh.times[tt] <-
          if (length(wh.ymd)>0)
            { if (length(wh.ymd)>1) warning(paste("Multiple matches: ",times[tt],sep=''),immediate.=TRUE)
              wh.ymd
            } else {
              if (length(wh.ym)>0)
                { if (length(wh.ym)>1) warning(paste("Multiple matches: ",times[tt],sep=''),immediate.=TRUE)
                  wh.ym
                }  else {
                  if (length(wh.y)>1) warning(paste("Multiple matches: ",times[tt],sep=''),immediate.=TRUE)
                  wh.y
                }
            }
      }
  }
  wh.times <- as.numeric(wh.times)
  
  
  ## loop through the required times
  for (tt in 1:length(wh.times)) {
    wh.keep=which( !is.na(as.vector(a.data(st)[,wh.times[tt]])) )

    dat.dum <- (a.data(st)[,wh.times[tt]])[wh.keep]
    dat <- if (tt==1) dat.dum else c(dat, dat.dum)    

    lon.dum=a.lon(st)[wh.keep]
    lon <- if (tt==1) lon.dum else c(lon, lon.dum)
    
    lat.dum=a.lat(st)[wh.keep]
    lat <- if (tt==1) lat.dum else c(lat, lat.dum)

    date.dum=rep(times.char[tt], length(wh.keep))
    date <- if (tt==1) date.dum else c(date,date.dum)
  }

  gg.frame <- data.frame(lon=lon, lat=lat,
                         date=factor(date, levels=times.char),
                         value=cut(dat, breaks, inc=T))
  gg.out <- if (!plot) ggplot( gg.frame ) else ggplot()

  if (length(wh.times)>1) gg.out <- gg.out + facet_wrap(~ date, ncol=ncol)
  
  if (plot) print(gg.out + geom_tile(data=gg.frame, aes(x=lon,y=lat,fill=value),alpha=.7) +
#                           xlim(min(gg.frame$lon),max(gg.frame$lon)) +
#                           ylim(min(gg.frame$lat),max(gg.frame$lat))  +
                           labs(fill=a.data.name(st)) +
                           world.map(add=world, lon=gg.frame$lon, lat=gg.frame$lat, ...)
                  )
  
  invisible(gg.out)
  
}
