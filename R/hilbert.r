## discrete time "analytic signal" via hilbert transform (via fft)

## input: real timeseries or matrix with time along columns (nspaceXntime).
##        any imaginary component is discarded.
## output: the analytic signal = input + i*HILBERT(input) for each point/timeseries,
## a complex spaceXtime result.

## The hilbert transform gives a quarter phase shift to the input timeseries.

## license: GPL license
## caveat emptor: no warranty whatsoever, use this code at your own risk.
## author: james mccreight # gmail o com

## reference:
## S. Lawrence Marple, Jr., Computing the discrete-time analytic 
##     signal via FFT, IEEE Transactions on Signal Processing, Vol. 47, 
##     No. 9, September 1999, pp.2600--2603.

## (this can be verified against hilbert() for univariate ts in the EMD package,
##   if that's any solace.)

hilbert <- function (input)  
{

  if (!is.real(input)) {
    print("*** Warning: discaring imaginary part of input.")
    input <- Re(input)
  }

  if (is.array(input)) {
    nspace <- dim(as.matrix(input))[1]
    ntime <- dim(as.matrix(input))[2]
  } else {
    ntime <- length(input)
    nspace <- 1
    input=t(input)
  }
  
  mult <- if (ntime%%2==0) {                        ## as per Marple's standard transform.
    c(1, rep(2, ntime/2-1), 1, rep(0, ntime/2-1) )  ## even length time series
  } else {
    c(1, rep(2, (ntime-1)/2 ), rep(0, (ntime-1)/2 ) )  ## odd length time series
  }

  mult=as.array(replicate(nspace, mult), dim=c(nspace,ntime))

  ## mvfft computs the fft of each column rather than a 2-d fft.
  ## but it's annoying/i'm stupid on a single vector TS for some reason, hence the if.
  if (nspace==1)
    return( fft( t(mult) * fft(input), inverse=TRUE) / ntime )
  else 
    return( t( mvfft( mult * mvfft(t(input)) , inverse=TRUE) / ntime ) )  ## inverse fft is unscaled in R.

}
