

ggCorrMap <- function( cm, only.signif=TRUE, world=TRUE, plot=TRUE, ... )
{

  df <- data.frame( p=a.p(cm),
                    est=a.est(cm),
                    lon=a.lon(cm),
                    lat=a.lat(cm) )
  if (only.signif) df <- df[which(a.is.signif(cm)),]
  
  gg.out <- ggplot( df )

  if (plot) print(gg.out + geom_tile(aes(x=lon,y=lat,fill=est)) +
                  scale_fill_gradient2(name=paste( attributes(a.est(cm)[1])$names,'at\n',a.confidence.level(cm),'confidence' ),
                                       midpoint=0, low='blue', high='red' ) +
                  xlim(min(a.lon(cm)),max(a.lon(cm))) + ylim(min(a.lat(cm)),max(a.lat(cm))) +
                  opts(title=paste(a.method(cm),' of\n',a.spaceTime.name(cm),' with\n',a.timeSeries.name(cm),'\n',
                                   a.time.range(cm), sep='') ) + world.map(world, ...)
                  )

}

