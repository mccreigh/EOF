## Take a wavelet filter (not yet a formal object) and comupte the hilbert transform on
## each of the filtered timeseries, return the input list  where the filtered
## timeseries are replaced these by their "analytic" counterparts.

analyticFilter <- function( wf ) {
  nbands=length( wf$window.labels.out )
  ## the as real suppresses the warnings...
  for (bb in 1:nbands) wf$filter[bb,] <- hilbert(as.real(wf$filter[bb,]))
  wf
}
