default.train.start <- function(temporal, EOF, CEOF) {
  get.time1 <- function(x) as.vector(x[1])  
  temp.t1 <- if (temporal)
    as.POSIXct(max(laply(llply(temporal.predictor.list,a.POSIXct), get.time1)),origin=origin) else NA
  EOF.t1 <- if (EOF)
    as.POSIXct(max(laply(llply(spatial.EOF.predictor.list,a.POSIXct), get.time1)),origin=origin) else NA
  CEOF.t1 <- if (CEOF)
    as.POSIXct(max(laply(llply(spatial.CEOF.predictor.list,a.POSIXct), get.time1)),origin=origin) else NA
  as.POSIXct( max( c(temp.t1, EOF.t1, CEOF.t1), na.rm=T) + origin )
}

default.forecast.end <- function(temporal, EOF, CEOF) {
  get.t.end <- function(x) as.vector(x[length(x)])
  temp.f.end <- if (temporal)
    as.POSIXct(min(laply(llply(temporal.predictor.list,a.POSIXct), get.t.end)), origin=origin) else NA
  EOF.f.end <- if (EOF)
    as.POSIXct(min(laply(llply(spatial.EOF.predictor.list,a.POSIXct), get.t.end)), origin=origin) else NA
  CEOF.f.end <- if (CEOF)
    as.POSIXct(min(laply(llply(spatial.CEOF.predictor.list,a.POSIXct), get.t.end)), origin=origin) else NA
  as.POSIXct( min( c(temp.f.end, EOF.f.end, CEOF.f.end), na.rm=T) + origin )
}

default.forecast.start <- function(temporal, EOF, CEOF) {
  temp.f1 <- if  (temporal)
  a.POSIXct(temporal.predictor.list[[1]])[round(mean(match(c(train.start, forecast.end),
                                                               temporal.predictor.list[[1]]@POSIXct)))] else NA
  EOF.f1 <- if  (EOF)
  a.POSIXct(spatial.EOF.predictor.list[[1]])[round(mean(match(c(train.start, forecast.end),
                                                               spatial.EOF.predictor.list[[1]]@POSIXct)))] else NA
  CEOF.f1 <- if  (CEOF)
  a.POSIXct(spatial.CEOF.predictor.list[[1]])[round(mean(match(c(train.start, forecast.end),
                                                               spatial.CEOF.predictor.list[[1]]@POSIXct)))] else NA
  as.POSIXct( max( c(temp.f1, EOF.f1, CEOF.f1), na.rm=T) + origin )  
}

mk.lead.time.list <- function(lead.times) {
  lead.time.list <- as.list(lead.times)
  ntime=length(target$POSIXct);  dY <- get.dY(target$POSIXct)
  ## putting the dates in becomes a bit of work since the forecast lead times may not be in the predictor set at all.
  ## months numerically in [0,11], so using the forecast.time month actually gives the next month numerically.
  names(lead.time.list) <- evenPOSIXct( length(lead.times),
                                 origin=as.numeric(format(target$POSIXct[ntime],'%Y')) +
                                 (as.numeric(format(target$POSIXct[ntime],'%m')))*dY+(dY/2), dY)
  lead.time.list
}

mk.band.list <- function() { ## looks at target, temporal.preds, EOF.preds, and CEOF.preds in the calling frame.
  ## there could be disparity in the number of bands between target, temporal, EOF, and EOF depending if
  ## any are global or not.
  nbands.target <-length(target$filter$filter[,1])
  nbands.temporal <- if (e.temporal.preds) length( temporal.preds$filter$filter[,1] ) else NA
  nbands.EOF <- if (e.EOF.preds) length( EOF.preds$ts ) else NA
  nbands.CEOF <- if (e.CEOF.preds) length( CEOF.preds$ts ) else NA
  nbands.all <- c( nbands.target, nbands.temporal, nbands.EOF, nbands.CEOF)
  nbands <- min(nbands.all, na.rm=TRUE )
  which.min <- match(nbands, nbands.all)[1]
  band.list <- as.list(1:nbands)
  names(band.list) <- eval(parse(text= c('target$filter$window.labels.out',
                                         'temporal.preds$filter$window.labels.out',
                                         'names(EOF.preds$ts)', 'names(CEOF.preds$ts)')[which.min]))[1:nbands]
  band.list
}


review.all.gws <- function( all.gws ) {
  all.gws <- ldply(all.gws) ## data.frame
  all.gws$.id <- as.numeric(factor(all.gws$.id))
  TransRevLog2 <<- Trans$new("revlog2", function(x) (log(x, 2)), function(x) 2^(x), function(x) bquote(2^.(x)))  
  breaks=2^c(-2:8); lab.breaks=format(breaks, digits=1, nsmall=2)
  print(ggplot( all.gws, aes( x=period, y=gws, color=.id, group=.id ) ) + geom_line() +
        scale_x_continuous(breaks=breaks, labels=lab.breaks, name="Period (years)",trans="RevLog2") +
        opts(legend.position='none'))
  print("Set spectral bands over all forecasts? (y/n):")
  if (tolower(substr(readLines(n=1),1,1))=='y') {
    windows <-  
      unlist(choose.lapply( forecast.times[length(forecast.times)], setup.target, review.all.gws=TRUE, nproc=1 ))
  }
  windows
}

#######################################################################
## target/predictand
## subset the data to the training set.
## subset the target.ts to the current training set including the desired lead times.
setup.target <- function(forecast.time,  just.gws=FALSE, review.all.gws.first=FALSE ) {

  train.end <- forecast.time #a.POSIXct(target.ts)[which( a.POSIXct(target.ts) == forecast.time ) -1]
  
  if (!exists('target')) {  # target exists in a parent environ if global.target was TRUE
    ## set up the training period on the target   
    target.subset <- ts.subset( target.ts, format(train.start,format='%m/%d/%Y'),format(train.end,format='%m/%d/%Y') )
    
    ## wavelet transform the target to understand our problem spectrally
    wt <- waveletTransform(target.subset, pad=TRUE, mother="morlet", recon=TRUE, verbose=FALSE)
    
    ## just return period and gws if looking at gws temporal non-stationarity over the forecast period
    if (just.gws) return(data.frame(period=a.period(wt),gws=rowMeans(abs(a.transform(wt))^2)))
    
    ## Set up the spectral windows on the target a) at each time, or b) globally.
    ## interactively decide on spectral bands/wavelet filter of the target.ts    
    accept.bands=TRUE
    if (review.all.gws.first) accept.bands=FALSE
    if (!accept.bands) ggFilter( wt, use.current.window=TRUE )
    while (!accept.bands) {
      print('Enter windows as a vector or an (n.window X 2) matrix:')
      windows=eval(parse(text= readLines(n=1)))
      target.filter <- waveletFilter( wt, window=windows )     
      ggFilter( wt, target.filter, use.current.window=TRUE )
      print(windows)
      print('accept? y/n:')
      if (tolower(substr(readLines(n=1),1,1))=='y') accept.bands=TRUE      
    }
  
  ## if choosen to set the gws globally, windows are set while looking at the last gws
  if (review.all.gws.first) return(windows)
  
  ## if set globally, calculate the filter now
  if (!exists('target')) target.filter <- waveletFilter( wt, window=windows )

  list( wt=wt, filter=target.filter, POSIXct=a.POSIXct(a.input(wt)) )
  } else {
    stop()
    ## have to subset the existing global target to the training period
#    train.start.ind <- which( == train.start)
#    train.end.ind <- which( == train.end)
    ##subset... both wt and wt.filter?
  }
  
}



#######################################################################
## predictor.list.CEOF
c2r <- function(c) { if (is.complex(c)) data.frame(r=Re(c),i=Im(c),a=Mod(c),p=Arg(c)) else data.frame(ts=c) }

setup.spatial <- function( forecast.time ) {

  is.EOF <- if () TRUE else FALSE
  is.CEOF <- if () TRUE else FALSE
  EOF.string <- if (is.EOF) 'EOF' else 'CEOF'
  save.path <- get(paste("global.",EOF.string,".save.path",sep=''))
  
  train.end <- forecast.time #a.POSIXct(target.ts)[which( a.POSIXct(target.ts) == forecast.time ) -1]

  if ( (is.EOF & !exists('EOF.preds')) | (is.CEOF & !exists('CEOF.preds')) ) {
    
    ## if save path and the file exists, load it
    if ((save.path)!='' & file_test('-f',save.path)) {        
      
      print(paste('loading:', save.path))
      load(save.path)
      preds.specData <- get(paste('preds.spec',EOF.string,sep=''))[wh.preds]
      ## check the restored variables against those passed in, keep only what was passed and in that order
      in.pred.names <- names(get(paste('spatial.',EOF.string,'.predictor.list'sep='')))
      load.pred.names <- preds.specData
      wh.preds <- match( in.pred.names, load.pred.names )
      if (any(is.na(wh.preds)))
          warning(paste('Some spatial.', EOF.string, 'predictor.list vars no in the restored global.',
                        , EOF.string, '.save.path file. stopping.',sep=''), .immediate=TRUE)
      ## temporal subsetting happens when the timeseries are extracted.
      
    } else {
      
      print(paste('setting up ',EOF.string,' predictors.'))
      ## subset to the training period
      make.pred.subset <- function(pred) {
        subset(pred, format(train.start,format='%m/%d/%Y'), format(train.end,format='%m/%d/%Y') ) }
      preds.subset <- llply( get(paste('spatial.',EOF.string,'.predictor.list',sep='')), make.pred.subset )
      
      ## filter the predictor fields.
      ## this can use a large number of processors, up to the number of spatial points in each field.
      do.filter.preds <- function(pred) {
        filterField( pred, target$wt, target$filter, nprocessors=nprocessors/length(preds.subset) ) }
      preds.filter <- choose.lapply( preds.subset, do.filter.preds, nprocessors=min(nprocessors,length(preds.subset) ))
      
      ## calculate the EOF/CEOF of each band.
      ## spectralCEOF uses upto a processor for each band and we have length(preds.filter) number of
      ## predictor variables, their product is the ideal # of processors.
      do.specEOF <- function(pred)
            eval(parse(text=paste('spectral',EOF.string,'(pred, corr=FALSE, nproc=nbands)',sep='')))
      preds.specData <- choose.lapply(preds.filter, do.specEOF, nproc=min(nprocessors, nbands*length(preds.filter)))
      
      ## save if the save path was specified and this is not an overwrite.
      if ((save.path)!='' & !file_test('-f',save.path)) save(preds.specData,file=save.path)

    }

    ## extract the CEOF timeseries over all bands and in the specified modes
    nbands=length(preds.specData[[1]]@band)
    ## create a list of CEOF modes if none was specified.
    ## else force the names to be ordered in the wh.modes list as in the predictor list, so we can mapply
    wh.modes <-  get(paste('wh.',EOF.string,'.modes',sep=''))
    if (!is.list(wh.modes))  {
      wh.modes <- as.list( rep(wh.modes,nbands) )
      names(wh.modes) <- names(preds.specData)
    } else { 
      wh.names <- match( names( preds.specData ), names(wh.modes) )
      if (any(is.na(wh.names)))
        warning('Improper names supplied for wh.',EOF.string,'.modes. stopping.', .immediate=TRUE)
      wh.modes <- wh.modes[wh.names]
    }
    
    get.bands <- function( bb ) {
      get.modes <- function( vv ) {
        ##check desired modes really exist, take those which do.
        wh.wh.modes = wh.modes[[vv]][which(wh.modes[[vv]] %in% 1:length(EOFpred[[vv]][1,]))]
        dum <- apply(as.matrix(EOFpred[[vv]][,wh.wh.modes]),2,c2r) ## real to real if not complex
        names(dum) <- paste("mode.", wh.wh.modes,sep='')
        as.data.frame(dum)
      }
      ## the the same band over all vars
      EOFpred <- lapply(preds.specData,function(l) a.timeseries(l@band[names(preds.specData[[1]]@band)[[bb]]][[1]]))
      names.list <- names(EOFpred); names(names.list) <- names(EOFpred)
      as.data.frame( llply( names.list , get.modes ) )
    }
    band.num.list <- as.list(1:nbands); names(band.num.list) <- names(preds.specData[[1]]@band)
    EOF.preds <- llply( band.num.list, get.bands )
                                
    EOF.POSIXct <- a.POSIXct(preds.specData[[1]])
    ## this is not redundant in the case of predictors loaded from file
    keep.inds <- which(EOF.POSIXct == train.start):which(EOF.POSIXct == train.end)
    EOF.preds <-llply( EOF.preds, function(b) b[keep.inds,])
    EOF.POSIXct <- a.POSIXct(preds.specData[[1]])[ keep.inds ]
    return(list( ts=EOF.preds, POSIXct=EOF.POSIXct ))
    
  } else {
    ## subset the global EOF.preds to the current train.end
    EOF.POSIXct <- EOF.preds$POSIXct
    keep.inds <- which(EOF.POSIXct == train.start):which(EOF.POSIXct == train.end)
    EOF.preds$ts <-llply( EOF.preds$ts, function(b) b[keep.inds,])
    EOF.preds$POSIXct <- EOF.POSIXct[ keep.inds ]
    return( EOF.preds )

  }
  
}


#####################################################################

forecast.leads <- function( lead ) {
  
  forecast.bands <- function( band ) {
    ## shifting the predictand timeseries into the past by lead wrt to predictors.
    ntime <- length( target$filter$filter[1,] )    
    model.frame <- data.frame( target= target$filter$filter[ band, (1+lead):ntime ] )

    print('experimental lag hard-coded')
    model.frame <- cbind( model.frame, data.frame( target= target$filter$filter[ band, 1:(ntime-lead) ] ) )

    pred.frame <- model.frame[1,] * NA
    if (e.temporal.preds) { model.frame <- cbind( model.frame, target$filter[ band, 1:(ntime-lead) ] )
                            pred.frame <- cbind( pred.frame, target$filter[ band, ntime ] ) }
    if (e.EOF.preds) { model.frame <- cbind( model.frame, EOF.preds$ts[[band]][1:(ntime-lead),] )
                       pred.frame <- cbind( pred.frame, EOF.preds$ts[[band]][ntime,] ) }
    if (e.CEOF.preds) { model.frame <- cbind( model.frame, CEOF.preds$ts[[band]][1:(ntime-lead),] )
                        pred.frame <- cbind( pred.frame, CEOF.preds$ts[[band]][ntime,] ) }
    band.forecast <- eval(call(model, model.frame, pred.frame,
                               confidence.level=confidence.level, return.quantiles=return.quantiles,
                               nsamp.conf=nsamp.conf, nprocessors=nprocessors ))
  }
  band.preds <- choose.lapply( band.list, forecast.bands, nprocessors=1)

  ## combine resampled confidence to get an ensemble, NOTE this is recursive.
  comb.ens <- function(l) { n=length(l); if (n==1) return(l[[1]]); outer(l[[1]], comb.ens(l[2:n]), '+') }
  full.ens <- quantile(as.vector(comb.ens( band.preds )), return.quantiles)
  
}

spectral.regression <- function(model.frame, pred.frame,
                               confidence.level=.95, return.quantiles=c(.05, .25, .5, .75, .95),
                               nsamp.conf=5, nprocessors=1)
{
  the.lm <- lm( target ~ ., model.frame )
  ## calculate the prediction confidence.
  ntime=length(model.frame[,1])
  pred.lm <- predict( the.lm, pred.frame, interval='prediction', level=confidence.level )  
  ## resample the confidence interval in each to yield an ensemble of estimates??
  ens <- if (nsamp.conf>2) seq(pred.lm[2],pred.lm[3],length.out=nsamp.conf) else pred.lm[1] 
}
