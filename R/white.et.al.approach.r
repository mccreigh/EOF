## this is barely started and incomplete, their methodology is difficult to make sense of completely,
## though it seems good on the surface.

#####################################################################
## CEOF phase space matching
#

## Each band of the target time series is cross-correlated against a real time series of the corresponding  band of a
## predictor's spectral CEOF. In each band, the spectral CEOF timeseries are then shifted in time (within 1/2 cycle) to
## give maximum correlation/anticorrelation with the target. Thi moves all components (real, complex, phase, amp) which
## requires a compensating phase shift of each spectral CEOF phase series to keep phase at its original location in
## time. The result of this temproal and phase shift yields target variable maximum/minimum at 0/180 deg of CEOF
## timeseries phase.

## a general routine should accept the filtered target, the spectralCEOF, number of modes, option to use analytic signal? vs only use the real part against the filter target.

#####################################################################
## inputs
filter <- rn.filter
specCEOF <- slp.specCEOF
modes=1
dY=12 ## can be obtained from the wavelet transform @dt, alternatively, i could include this with the filter... 
## ...

################################### main function body
filter.hilbert <- analyticFilter( rn.filter )
nbands <- length(rn.filter$window.out[,1])

for (mm in modes) {
  
  specCEOF.mode.ts <- get.spectralCEOF.by.mode(specCEOF, mode=mm)$time

  ## do over all bands via lapply:

  get.ccf.trans.shift <- function( bb, imaginary=FALSE, ... ) {
    mean.period.years=mean(filter.hilbert$window.out[bb,])
    band.ccf <- ccf( eval(call( if (!imaginary) "Im" else "Re", filter.hilbert$filter[bb,])),
                    specCEOF.mode.ts[[ if (!imaginary) "imag" else "real" ]][[bb]],
                    lag.max = mean.period.years*.5*(1/dY), plot = TRUE, ...)
    wh.max <- which(abs(band.ccf$acf) == max(abs(band.ccf$acf)))
    ccf.max <- band.ccf$acf[wh.max]
    ccf.lag.max <- band.ccf$lag[wh.max]
    phase.shift.radians <- (ccf.lag.max*dY)/(mean.period.years) * 360 #2 *pi
    list( ccf.max=ccf.max, ccf.lag.max=ccf.lag.max, phase.shift.radians=phase.shift.radians )
  }

  ## ccf the real parts, translate and phase shift
  trans.shift.real <- lapply( as.list(1:nbands), get.ccf.trans.shift )
  trans.shift.imag <- lapply( as.list(1:nbands), get.ccf.trans.shift, imaginary=TRUE )
  
  ## ccf the imaginary parts

  
  
