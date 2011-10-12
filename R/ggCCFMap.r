## i havent tested this with only.signif=FALSE

ggCCFMap <- function( ccfMap, only.signif=TRUE, world=TRUE, plot=TRUE, ... )
{

  map.length=length(a.lon(ccfMap))
  
  df.corr <-  data.frame( value=center.discrete( c(a.max.corr(ccfMap), -1*a.min.corr(ccfMap)), step=.1, include.lowest=TRUE),
                          variable=factor( c( rep('Correlation',map.length), rep('Anti-correlation',map.length) ), levels=c('Correlation','Anti-correlation') ), 
                          lon=c(a.lon(ccfMap),a.lon(ccfMap)),
                          lat=c(a.lat(ccfMap),a.lat(ccfMap)) )

  all.lags <-  c(a.lag.max(ccfMap), a.lag.min(ccfMap))
  df.lag <-  data.frame( value=center.discrete( all.lags, step=2, odd=FALSE, include.lowest=TRUE ),
                         variable=factor( c( rep('Lag.Corr',map.length), rep('Lag.Anti',map.length) ), levels=c('Lag.Corr','Lag.Anti') ), 
                         lon=c(a.lon(ccfMap),a.lon(ccfMap)),
                         lat=c(a.lat(ccfMap),a.lat(ccfMap)) )
                        
  df.signif <-  data.frame( value=c(a.is.signif.max(ccfMap), a.is.signif.min(ccfMap)),
                           variable=factor( c( rep('Signif.Corr',map.length), rep('Signif.Anti',map.length) ), levels=c('Signif.Corr','Signif.Anti') ), 
                           lon=c(a.lon(ccfMap),a.lon(ccfMap)),
                           lat=c(a.lat(ccfMap),a.lat(ccfMap)) )
  
  if (only.signif) {
    df.corr <- df.corr[which(df.signif$value),]
    df.lag <- df.lag[which(df.signif$value),]
    df.signif <- df.signif[which(df.signif$value),]
  }

  ncorrbreaks=length(levels(df.corr$value))
  nlagbreaks=length(levels(df.lag$value))
  
  if (plot) {

    gg.corr <- ggplot(df.corr) + geom_tile(aes(x=lon,y=lat,fill=value)) + world.map() + facet_wrap( ~variable) +
                 ##scale_fill_gradient2(name=paste('Correlation at\n',a.confidence.level(ccfMap),'confidence' ), midpoint=0, low='blue', high='red', mid='yellow' ) +
                 scale_fill_manual(values=rev(c( rev(brewer.pal(ncorrbreaks %/% 2 +1,'YlGnBu')[-1]), brewer.pal(ncorrbreaks %/% 2+1,'YlOrRd')  ) ),
                                   name='(Anti)Correlation') +
                 opts(title=paste('Correlation of ',a.spaceTime.name(ccfMap),' with ', a.timeSeries.name(ccfMap),' ', a.time.range(ccfMap), sep='') )  +
                 xlim(min(a.lon(ccfMap)),max(a.lon(ccfMap))) + ylim(min(a.lat(ccfMap)),max(a.lat(ccfMap))) 

    gg.lag <- ggplot(df.lag) + geom_tile(aes(x=lon,y=lat,fill=value)) + world.map() + facet_wrap( ~variable) +
                 ##scale_fill_gradient2(name=paste('Lag' ), midpoint=0, low='blue', high='red', mid='yellow' ) +
                 scale_fill_manual(values=rev(c( rev(brewer.pal(nlagbreaks %/% 2 +1,'YlGnBu')[-1]), brewer.pal(nlagbreaks %/% 2+1,'YlOrRd')[-1]  ) ),
                                   name='Lag') +
                 opts(title=paste('Shift of ',a.spaceTime.name(ccfMap),' w.r.t. ', a.timeSeries.name(ccfMap),' ', a.time.range(ccfMap), sep='') ) +
                 xlim(min(a.lon(ccfMap)),max(a.lon(ccfMap))) + ylim(min(a.lat(ccfMap)),max(a.lat(ccfMap))) 

    check.new.dev()

    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 1)))
    vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col=y)
    ## plot into view ports
    print(gg.corr,   vp=vplayout( 1 , 1 ) )
    print(gg.lag, vp=vplayout( 2 , 1 ) )
    
  }

}

