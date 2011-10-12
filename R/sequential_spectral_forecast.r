## forecast times are the times of the "last available data", when the forecast is made on the
## latest possible data, the same as train.end. lead times extend forward in time from here by
## interger multiples of the timestep.

## Spectral params are selected and the model fit sqeuentially,
## starting with some training period and advancing in time.
## Before fitting a spectral model (e.g. regression), some prelim model params and its consequences:
## 1. spectral bands (params) are selected for wavelet filtering (either once and for all or at each timestep)
## 2. the target/predictand is filtered on the selected bands at each time step,
## 2. the spaceTime predictor fields are filtered on the same bands (at each point),
## 3. a CEOF is computed on each spaceTime fields spectral bands.
## Then these are passed to a spectral model which relates the predictors and predictands.

## option to lag predictors via a list name=lag.vector ??

sequential.spectral.forecast <-
  function( model='spectral.regression',
           confidence.level=.95,
           nsamp.conf=3,
           return.quantiles=c(.05,.25,.5,.75,.95),
           target.ts,
           temporal.predictor.list,
           spatial.EOF.predictor.list,
           spatial.CEOF.predictor.list,
           train.start=default.train.start(e.temporal.preds,e.EOF.preds,e.CEOF.preds),
           forecast.end=default.forecast.end(e.temporal.preds, e.EOF.preds, e.CEOF.preds),
           forecast.start=default.forecast.start(e.temporal.preds, e.EOF.preds, e.CEOF.preds),
           forecast.lead.times=1,
           forecast.step=1,
           review.all.gws.first=TRUE,
           spectral.windows,
           global.windows=TRUE,           
           global.target=TRUE,
           global.temporal.preds=TRUE,
           global.EOF.preds=TRUE,
           global.CEOF.preds=TRUE,
           global.EOF.save.path='',
           global.CEOF.save.path='',
           wh.EOF.modes=1,
           wh.CEOF.modes=1,
           nprocessors=1
           )
{

  ## note which predictor sets we are dealing with here.
  e.temporal.preds=!missing(temporal.predictor.list)
  e.EOF.preds=!missing(spatial.EOF.predictor.list)
  e.CEOF.preds=!missing(spatial.CEOF.predictor.list)
  e.spatial.preds <- e.EOF.preds | e.CEOF.preds
  
  source("/u/wk/jmccreig/R/jlm_lib/EOF/sequential_spectral_forecast_routines.r", local=TRUE)  

  ## this is the routine eventually called to run forecast times in parallel.
  spectral.forecast <- function( forecast.time ) {

    #######################################################################
    ## Set up the data at the current forecast.time

    ## training set runs from start of time to 1 step before the forecast time, all the data prior to "present"
    ##train.end = a.POSIXct(target.ts)[which( a.POSIXct(target.ts) == forecast.time ) -1]
    ## the maximum number of spectral bands in use.
    ##nbands = dim( as.matrix(windows) )[1]

    ## these functions have access to the variables defined in the parent env to the current function.
    ## The target: band filter
    target <- setup.target( forecast.time )    
    ## Temporal predictors: band filter each
    if (e.temporal.preds) temporal.preds <- setup.temporal( forecast.time )    
    ## Spatial EOF predictors: filter, EOF, collect timeseries
    if (e.EOF.preds | e.CEOF.preds) spatial.preds <- setup.spatial( forecast.time )    
    
    ## #####################################################################
    ## call the spectral model on each lag for the current forecast time.    
    ## first need to bring these functions into the current environment/frame.
    environment(mk.lead.time.list) <- environment()
    environment(forecast.leads) <- environment()
    environment(mk.band.list) <- environment()
    
    lead.time.list <- mk.lead.time.list(forecast.lead.times)
    band.list <- mk.band.list() ## this is available inside forecast.leads
    choose.lapply( lead.time.list, forecast.leads, nprocessors=1 )
  }

  ## Begin the main routine.
  ## #####################################################################
  ## Set up the desired forecast.times which represent the "present" - when the forecast is made,
  ## they allow the fullest possible training period not using data beyond 
  ## this time. The end of training set and beginning of the response are cut to accomodate lags.
  wh.target.start <- match( forecast.start, a.POSIXct(target.ts))
  wh.target.end <- match( forecast.end, a.POSIXct(target.ts))
  forecast.times <- as.list( a.POSIXct(target.ts)[seq(wh.target.start, wh.target.end,by=forecast.step)] )
  names(forecast.times) <- a.POSIXct(target.ts)[seq(wh.target.start, wh.target.end,by=forecast.step)]

  ## #####################################################################
  ## Global settings.

  ## Target  
  ## If global.windows are not supplied then the non-stationarity of the gws can be examined
  ## and then global windows set by insepection of the gws at the final time.
  if (missing(spectral.windows) & global.windows) {
    all.gws <- choose.lapply( forecast.times, setup.target, just.gws=TRUE, nproc=nprocessors)
    windows <- review.all.gws( all.gws )
  } else if (global.windows) windows <- spectral.windows
  if (global.target)
    target <- setup.target( forecast.times[length(forecast.times)] )
  ## temporal preds
  if (e.temporal.preds & global.temporal.preds)
    temporal.preds <- setup.temporal( forecast.times[length(forecast.times)] )
  ## EOF/CEOF preds
  if (e.spatial.preds) {
    spatial.preds <- list()
    if (e.EOF.preds & global.EOF.preds) spatial.preds[1] <- 'EOF'
    if (e.CEOF.preds & global.CEOF.preds) spatial.preds[length(spatial.preds)] <- 'CEOF'
    if (length(spatial.preds) > 0) setup.spatial( forecast.times[length(forecast.times)] )
  }

  ## #####################################################################
  ## forecast! nprocessors here is the number of parallelizations of the above
  print('calling forecast model')
  choose.lapply( forecast.times , spectral.forecast, nprocessors=nprocessors, mc.preschedule=TRUE )
  
}
  
