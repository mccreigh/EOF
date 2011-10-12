require(ggplot2) 

if (FALSE) {
  ## this is how you'd read the ncdf file if you didnt have the data already saved
  ## in the EOF package in r binary format.
  require(ncdf)
  kaplan.sst <- open.ncdf('~/methods/cluster_lecture/kaplan.sst.ncdf')
  str(kaplan.sst) ## details, unix equivalent: ncdump -h kaplan.sst.ncdf
  print(names(kaplan.sst$dim))
  print(names(kaplan.sst$var))
  
  sst <- get.var.ncdf( kaplan.sst, 'ssta' )
  lon <- (get.var.ncdf( kaplan.sst, 'X' ) +360 ) %% 360 ##??
  lon[ lon>180 ] <- lon[ lon > 180 ] - 360
  lat <- get.var.ncdf( kaplan.sst, 'Y' )
  time <- get.var.ncdf( kaplan.sst, 'T' )/12 + 1960 ## odd units, dont really need
  close.ncdf(kaplan.sst)
  
  nlat=length(lat); nlon=length(lon); ntime=length(time)
  time.kaplan <- evenPOSIXct( ntime, origin=1856+1/24, dY=1/12 )
  lat.kaplan <- t(matrix( rep(lat,each=nlon), ncol=nlon, nrow=nlat, byrow=T ))
  lon.kaplan <- t(matrix( rep(lon,each=nlat), ncol=nlon, nrow=nlat, byrow=F ))
  
  ## flatten, put time on cols, space on rows of df
  ## putting POSIXct on column names is kinda futile, keep separate.
  kaplan <- adply(sst, c(1,2) ) 
  kaplan[1] <- as.vector(lon.kaplan)
  kaplan[2] <- as.vector(lat.kaplan)
  names(kaplan)[c(1,2)] <- c('lon','lat') 
  save(kaplan, time.kaplan, file='~/R/jlm_lib/EOF/data/kaplan_sst.rsav')
  
} 

load("~/R/jlm_lib/EOF/data/kaplan_sst.rsav") ## load the data from the EOF "package"
## keep only 1900-last
wh.keep <- which( format(time.kaplan,'%Y') >= 1900 )
kaplan <- kaplan[, c(1,2,wh.keep+2)]
time.kaplan <- time.kaplan[wh.keep]

## one time slice, looking good
ggplot( kaplan ) + geom_tile( aes(x=lon,y=lat,fill=V529) ) + world.map(posi=F) 

## make a space time object, as defined in the EOF package
ksst.ST <- new("spaceTime", data=as.matrix(kaplan[,c(-1,-2)]),
                 lon=kaplan$lon, lat=kaplan$lat, POSIXct=time.kaplan,
                 data.name='Kaplan SST', data.units='deg C')

## compare w the previous plot, could give it more slices in the second position... 
dev.new()
ggField( ksst.ST, 1 )

## a function which will show the results, not built in cause it would be too inflexible...
plot.eof <- function(gg.eof, ...) {  ## ... is passed to world.map
  s <- gg.eof$space+ geom_tile( aes(x=lon,y=lat,fill=values.space) ) +
    facet_wrap(~modes, ncol=1) +
    world.map( ...  ) +
    scale_fill_gradient2( low='blue', high='orange', midpoint=0 )
  t <- gg.eof$time + geom_line( aes(x=time, y=values.time) ) +
    geom_point( aes(x=time, y=values.time), color='red',alpha=.2, size=.3) +
    facet_wrap( ~modes, ncol=1 )
  v <- gg.eof$var + geom_point( aes(x=modes, y=variance.percent) )  
  nrow=6; ncol=10
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol)))
  vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col=y)
  print(s, vp=vplayout(1:6, 1:6) )
  print(t, vp=vplayout(1:4, 7:10) )
  print(v, vp=vplayout(5:6, 7:10) )
  invisible(s) ## trick for later
}

## perform a correlation based EOF on the kaplan monthly SSTs.
corr.ksst <- corrMtx(ksst.ST) 
eof.corr.ksst <- EOF(corr.ksst)
gg.eof.corr.ksst <- ggEOF(eof.corr.ksst)  
plot.eof(gg.eof.corr.ksst, lon=ksst.ST@lon,lat=ksst.ST@lat, posi=FALSE, buffer=10)

## el nino is the second mode and has the following amplitudes in spaceTime
el.nino <- eof.corr.ksst@EOF[,2] %o% eof.corr.ksst@timeseries[,2]

## incomplete timeseries (space) points were removed from the correlation calc
## the indices kept are eof.corr.ksst@mtx@A@wh.keep
## but we also want to subtract el.nino, which is from the anomaly space,
## from the anomalies... not the original data.
## note the anomalies spaceTime obj with the incomplete timeseres removed in
## eof.corr.ksst@mtx@A@anomaly
no.nino <- eof.corr.ksst@mtx@A@anomaly
no.nino@data <- no.nino@data -el.nino

## examine the correlation based EOF with el nino removed.
corr.nonino <- corrMtx(no.nino) 
eof.corr.nonino <- EOF(corr.nonino)
gg.eof.corr.nonino <- ggEOF(eof.corr.nonino) 
dev.new()
plot.eof(gg.eof.corr.nonino)
## the modes are the same, some signs are reversed. not surprising.


## get JJAS for 20S-20N
wh.jjas=which( as.numeric(format(ksst.ST@POSIXct,'%m')) >=6 &
               as.numeric(format(ksst.ST@POSIXct,'%m')) <=9 )
ksst.jjas.pac <- subset( ksst.ST, wh.jjas, dim='time')
wh.trop.pac <- which( ksst.ST@lat >= -20 &ksst.ST@lat <= 20 &
                      ( ksst.ST@lon >125 | ksst.ST@lon < -75 ) )                    
ksst.jjas.pac <- subset(ksst.jjas.pac, wh.trop.pac, dim='space' )
ksst.jjas.pac@lon <- (ksst.jjas.pac@lon + 360) %% 360

## use plyr to collapse the time dimension
## there's no september in 2011... so drop the last 3 
ksst.jjas.pac@data <-
  aaply(ksst.jjas.pac@data, 1,
        function(x) aggregate(x,list(format(ksst.jjas.pac@POSIXct,'%Y')), mean )$x )[,1:111]
ksst.jjas.pac@POSIXct <- evenPOSIXct( 111, origin=1900+7/12+1/24, dY=1 )

## examine the above work
ggField(ksst.jjas.pac, 1, buffer=20)

## take after 1975 only
wh.last <- which( ksst.jjas.pac@POSIXct > 1975)
ksst.jjas.pac.last <- subset(ksst.jjas.pac, wh.last, dim='time')

## correlation based EOF 
corr.ksst.jjas.pac.last <- corrMtx(ksst.jjas.pac.last)
eof.corr.ksst.jjas.pac.last <- EOF(corr.ksst.jjas.pac.last)
gg.eof.corr.ksst.jjas.pac.last <- ggEOF(eof.corr.ksst.jjas.pac.last)  ## plot the first 3 modes
## plot the results, just to see them and then
##take what we need for plotting the 3 spatial modes
zz <- plot.eof(gg.eof.corr.ksst.jjas.pac.last,
               lon=ksst.jjas.pac.last@lon,lat=ksst.jjas.pac.last@lat,
               posi=TRUE, buffer=10)$data
zz <-zz[which(zz$modes!='Mode.4'),]

##  plotting the 3 spatial modes
gg.eof.sp <- ggplot( zz ) + geom_tile( aes(x=lon,y=lat,fill=values.space) ) +
  facet_wrap(~modes, ncol=1) +
  world.map(posi=T, lon=zz$lon, lat=zz$lat, buffer=15 ) +
  scale_fill_gradient2( low='blue', high='orange', midpoint=0 ) +
  opts(title='EOF spatial patterns')
print(gg.eof.sp)

## will compare EOF and cluster "modes"
## look at the EOF modes
get.mode <- function(mode.number) {
  dum <- as.data.frame(t(eof.corr.ksst.jjas.pac.last@EOF[,mode.number] %o%
                         eof.corr.ksst.jjas.pac.last@timeseries[,mode.number]))
  dum$modes <- paste('Mode.',mode.number,sep='')
  dum$POSIXct <- eof.corr.ksst.jjas.pac.last@POSIXct
  melt(dum,id=c('modes','POSIXct'))
}
mode1 <- get.mode(1); mode2 <- get.mode(2); mode3 <- get.mode(3)
eof.modes <- rbind( Mode.1=mode1, Mode.2=mode2, Mode.3=mode3)

time.comps <- gg.eof.corr.ksst.jjas.pac.last$time$data ## cuz want to overplot time components
time.comps <- time.comps[which(time.comps$modes!='Mode.4'),]
time.comps$modes <- factor(time.comps$modes,levels=levels(time.comps$modes)[-4])
## normalize the time component for plotting purposes.
for (i in 1:3) time.comps$values.time[which(as.numeric(time.comps$modes)==i)] <-
  time.comps$values.time[which(as.numeric(time.comps$modes)==i)] *
  max(eof.corr.ksst.jjas.pac.last@EOF[,i]) #ugly, i see another way but it's comparable

## plot the EOF modes         
gg.eof.modes <-
  ggplot( eof.modes ) +
  geom_line( aes(x=POSIXct,y=value, color=modes, group=variable) ) +
  geom_line( data=time.comps, aes(x=time,y=values.time), size=2) +
  facet_wrap(~modes, ncol=1) + opts(title='EOF modes')
print(gg.eof.modes)
      
## cluster the data in to 3 groups. the following discard the NA points in space.
wh.keep <- eof.corr.ksst.jjas.pac.last@mtx@A@wh.keep
good.data <- ksst.jjas.pac.last@data[wh.keep,]
kclust <- kmeans( good.data , 3 )  ## actually cluster.
clust.frame <- as.data.frame(t(good.data))  ## transform data
names(clust.frame) <- 1:length(clust.frame)  ## will use to map the clusters onto the data
clust.frame$POSIXct <- eof.corr.ksst.jjas.pac.last@POSIXct
clust.frame <- melt(clust.frame,'POSIXct')
clust.frame$cluster <-
  factor(paste('Clust.',kclust$cluster[as.numeric(clust.frame$variable)],sep='')) ## boo-yah

clust.centers <- as.data.frame(t(kclust$centers)) ## we'll plot the centroids as well.
names(clust.centers) <- paste("Clust.",1:3,sep='')
clust.centers$POSIXct <- eof.corr.ksst.jjas.pac.last@POSIXct
clust.centers <- melt(clust.centers,id='POSIXct'); names(clust.centers)[2] <- 'cluster'

gg.clust <-
  ggplot( clust.frame ) +
  geom_line( aes(x=POSIXct,y=value, color=cluster, group=variable) ) +
  geom_line( data=clust.centers, aes(x=POSIXct,y=value), size=2) +
  facet_wrap( ~cluster, ncol=1 ) + opts(title='Clusters')
print(gg.clust)

## the cluster spatial pattern
clust.df <-  data.frame( value=factor(kclust$cluster), lon=ksst.jjas.pac.last@lon[wh.keep],
                        lat=ksst.jjas.pac.last@lat[wh.keep] )
gg.clust.sp <- ggplot( clust.df) + geom_tile( aes(x=lon,y=lat,fill=value) ) +
  world.map(posi=TRUE, lon=clust.df$lon, lat=clust.df$lat, buffer=15) +
  opts(title='Cluster pattern')

pdf(w=11,h=8.5,file="hw2_cluster_1.pdf")
nrow=4; ncol=3
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow,ncol)))
vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col=y)
print(gg.eof.modes, vp=vplayout(1:4, 1) )
print(gg.clust, vp=vplayout(1:4, 2) )
print(gg.eof.sp, vp=vplayout(1:3, 3) )
print(gg.clust.sp, vp=vplayout(4, 3) )
dev.off()
