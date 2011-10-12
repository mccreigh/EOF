######################################################################
## classes: corrMtx , covarMtx
## Correlation Matrix
## Because the B field is optional have an optional slot (B used in SVD / EEOF analysis).
setClassUnion("anomalyNULLmissing", c("anomaly", "NULL", "missing"))
setClassUnion("vectorNULLmissingLogical", c("vector","NULL","missing","logical"))
## corrMtx
setClass("corrMtx",
         representation(corrMtx="array",
                        A="anomaly", weight.type.A='character', weights.A='vectorNULLmissingLogical',
                        B="anomalyNULLmissing", weight.type.B='character', weights.B='vectorNULLmissingLogical',
                        complex="logical",
                        lag="numeric", fast="logical") )
set.accessors( 'corrMtx')
## corrMtx
setClass("covarMtx",
         representation(covarMtx="array",
                        A="anomaly", weight.type.A='character', weights.A='vectorNULLmissingLogical',
                        B="anomalyNULLmissing", weight.type.B='character', weights.B='vectorNULLmissingLogical',
                        complex='logical',
                        lag="numeric", fast="logical") )
set.accessors( 'covarMtx')

######################################################################
## methods corrMtx and covarMtx
## Calculate correlation and covariance matrices. can accept 1 or 2 matrices, A or B.

## fast:
## The fast options only applies when a single matrix is supplied. In this case, if the
## time dimension (2) is smaller than the space dimension, then the matrix is transposed
## prior to A%*%t(A), because t(A)%*%A is much smaller. If used this is noted by EOF
## which will recover the spatial patterns which are identical to the slower way.
## For notes on the "fast" method, see http://www.idlcoyote.com/code_tips/eof_analysis.html
## (I know that one of corr/covar could call the other while only altering a key word. They
## are kept separate with their own for clarity purposes later on.)

## weights (geographical):
## because many grids are not equal area, points at northern latitudes (closer together
## and representing smaller areas) can exhibit inflated correlation/covariance. weighting
## points by sqrt( abs( cos(pi * a.lat(get(wh.weight))/180) ) ) is the common cure.
## One can optionally apply weighting before corr/covar calculation. 
## Weights can be calculated on latitude (weights=TRUE) or be user
## supplied in an array with dimensions of space (ie if precise
## areas are known). I still need to make sure weighting preserves total variance...
## Note that the "weights" variable becomes "weight.type" and the actual
## weights go into "weights" in the code.
## The weights are returned with the corr/covarMtx object and they are also embedded in
## the calculated matrix.

## complex:
## create an "analytical" time series at each point in space. this is done by hilbert
## transform applied to the calculated anomalies, which yields the timeseries ts+H(ts)*i
## at each point. correlation and covariance matrices are calculated using the conjugate
## transpose. the complex timeseries is not returned, only the complex corr/covar matrix.
## (hilbert() of the anomalies or original fields can get one the "analytic" timeseries.)
## The compelex corr/covar matrices are used in CEOF to reveal information on propigating
## signals, including their spatial amplitude and phase for each mode. Phase sequences
## for the timeseries are also handled by the CEOF routine.

## lag:
## if two matrices are supplied, their correlation can be time lagged.
## B is lagged behind A using positive integers (timestep inferred from POSIXct),
## and A behind B by negative integers.

## anom.dim='time' (default) or 'space'. mean (and standard deviation if standardize=TRUE)
## are calculated over this dimension for all points in the other dimension and then removed.
## this is typically, though not always time. using 'space' may yield some bugs as it is
## currently under-tested.

## allow.missing:
## instead of removing points from a spaceTime field, missing points can be
## withheld in calculation of anomalies

if (!isGeneric("corrMtx")) {  ## creates a generic function
  fun <- if (is.function("corrMtx")) {
    corrMtx
  } else {
    function( A, B, ...) standardGeneric("corrMtx")
  }
  setGeneric("corrMtx", fun)
}

if (!isGeneric("covarMtx")) {  ## creates a generic function
  fun <- if (is.function("covarMtx")) {
    covarMtx
  } else {
    function( A, B, ...) standardGeneric("covarMtx")
  }
  setGeneric("covarMtx", fun)
}

setClassUnion("spaceTimeNULLmissing", c("spaceTime", "NULL", "missing"))

setMethod("corrMtx", c("spaceTime", "spaceTimeNULLmissing"),
          function( A, B, lag=0, fast=TRUE,  
                   weights=FALSE, allow.missing=FALSE, anom.dim='time', period.fmt='n/a', complex=FALSE,
                   weights.B=weights, allow.missing.B=allow.missing,
                   anom.dim.B=anom.dim, period.fmt.B=period.fmt, complex.B=complex, ... ) {
            invisible( corrCovarMtx('corrMtx', A, B, lag, fast, 
                                    weights, allow.missing, anom.dim, period.fmt, complex, 
                                    weights.B, allow.missing.B, anom.dim.B, period.fmt.B, complex.B ) )
          })

setMethod("covarMtx", c("spaceTime", "spaceTimeNULLmissing"),
          function( A, B, lag=0, fast=TRUE, 
                   weights=FALSE, allow.missing=FALSE, anom.dim='time', period.fmt='n/a', complex=FALSE,
                   weights.B=weights, allow.missing.B=allow.missing,
                   anom.dim.B=anom.dim, period.fmt.B=period.fmt, complex.B=complex, ... ) {
            invisible( corrCovarMtx('covarMtx', A, B, lag, fast, 
                                    weights, allow.missing, anom.dim, period.fmt, complex,
                                    weights.B, allow.missing.B, anom.dim.B, period.fmt.B, complex.B ) )
          })

## this routine is meant to be internal
corrCovarMtx <- function(mtxType, A, B, lag, fast,  
                         weights.A, allow.missing.A, anom.dim.A, period.fmt.A, complex.A,  
                         weights.B, allow.missing.B, anom.dim.B, period.fmt.B, complex.B ) {
  
  if (!missing(B)) {
    if ( dim(a.data(A))[2] == dim(a.data(B))[2] )
      warning(paste(mtxType,'failure: \nsecond dimensions of A and B are required to match'), immediate.=TRUE)
    if ( (complex.A-.5)*(complex.B-.5) < 1 )
      warning(paste(mtxType,'failure: \nBoth or neither of A and B must be complex.', immediate=TRUE))    
  }
  
  ###############
  ## Calculate the anomalies of the input fields.
  anom.A <- anomaly(A, standardize=(mtxType=='corrMtx'),
                    allow.missing=allow.missing.A, anom.dim=anom.dim.A, period.fmt=period.fmt.A)
  anom.B <- if (!missing(B)) anomaly(B, standardize=(mtxType=='corrMtx'),
                                     allow.missing=allow.missing.B, anom.dim=anom.dim.B, period.fmt=period.fmt.B) else NULL
  
  ## get the dimensions from the anomalies
  nspace.A=length(a.lon(a.anomaly(anom.A))); ntime.A=length(a.POSIXct(a.anomaly(anom.A)))
  if (!missing(B)) { nspace.B=length(a.lon(a.anomaly(anom.B))); ntime.B=length(a.POSIXct(a.anomaly(anom.B))) }

  ###############
  ## Geographical weighting
  for (wh.weight in c('A','B')) {

    weights <- eval(parse(text=paste('weights.',wh.weight,sep='')))
    
    weight.type <- if (typeof(weights) == 'logical') { if(weights) 'lat weights' else 'none'} else 'user supplied'
    if (wh.weight=='B' & missing(B)) weight.type <- 'none'
    
    if (weight.type!='none') {
      
      if (weight.type=='lat weights')
        { weights <- sqrt( abs( cos(pi * a.lat(get(wh.weight))/180) ) )
          weights[which(weights < .01)] <- NA }

      if (weight.type=='user supplied')
        { if (length(a.lat(get(wh.weight)))!=length(weights))
            { warning(paste(mtxType,': length of weights does not match spatial',
                            'dimension of ',wh.weight,'. No weighting used.',sep=''))
              weight.type <- 'none' }
        }
      
      ## set weights.A or weights.B
      if (weight.type!='none') weights <- weights[ a.wh.keep( get(paste('anom.',wh.weight,sep='')) ) ]
      print("warning: need to conserve the total field variance!!! ")  ## for now a todo reminder
      
    }

    if (weight.type=='none') weights <- NULL ## if none

    assign(paste("weight.type.",wh.weight,sep=''), weight.type)
    assign(paste("weights.",wh.weight,sep=''), weights)

  } ## for wh.weight

  ###############
  ## Time lagging
  if (lag==0){
    if (missing(B)) {
      ntime <- ntime.A
    } else {
      if ( all(a.POSIXct(a.anomaly(anom.A)) == a.POSIXct(a.anomaly(anom.B))) ) {  ## timeseries line up
        ntime <- ntime.A
      } else {
        warning("Timeseries dont have the same abscissa.", immediate.=TRUE)
      }
    }
  } else {
    warning("No lag routines have been written yet!", immediate.=T)
  }
    
  ###############
  
  if (missing(B) & fast) ## is fast option faster or not?
    if (length(a.POSIXct(a.anomaly(anom.A))) > length(a.lon(a.anomaly(anom.A)))) fast <- FALSE
  if (complex.A) fast <- FALSE  ## not allowing for CEOF since I'm not getting exactm atching results right now.
  
  ## calculate the A matrix considering geographical weighting and complex/hilbert transform.
  anom.weight.A <- if (weight.type.A=='none') {
    if (!exists('hilbert')) src("hilbert.r")
    if (!complex.A) a.data(a.anomaly(anom.A)) else hilbert(a.data(a.anomaly(anom.A)))
  } else {
    if (!complex.A) {
      a.data(a.anomaly(anom.A)) * t(array(data=replicate( ntime, weights.A ), dim=c(ntime,nspace.A)))
    } else {
      if (!exists('hilbert')) src("hilbert.r")
      hilbert(a.data(a.anomaly(anom.A))) * t(array(data=replicate( ntime, weights.A ), dim=c(ntime,nspace.A)))
    }
  }

  mtx <- if (missing(B)) {

    if (!fast) {
      if (!complex.A) anom.weight.A %*% t(anom.weight.A) else anom.weight.A %*% t( Conj(anom.weight.A)) ## conjugate transpose
    } else {
      if (!complex.A)  t(anom.weight.A) %*% anom.weight.A else t( anom.weight.A ) %*% (Conj(anom.weight.A))  ##does this work???
    }
    
  } else {

    ## calculate the B matrix considering geographical weighting and complex/hilbert transform.
    anom.weight.B <- if (weight.type.B=='none') {
      if (!complex.B) a.data(a.anomaly(anom.B)) else hilbert(a.data(a.anomaly(anom.B)))
    } else {
      if (!complex.B) {
        a.data(a.anomaly(anom.B)) * t(array(data=replicate( ntime, weights.B ), dim=c(ntime,nspace.B)))
      } else {
        hilbert(a.data(a.anomaly(anom.B))) * t(array(data=replicate( ntime, weights.B ), dim=c(ntime,nspace.B)))
      }
    }

    ## no fast calculation for cross covar/corr, but a cross complex conjugate transpose?? :)
    if (!complex.B) anom.weight.A %*% t(anom.weight.B) else anom.weight.A %*% t( Conj(anom.weight.B) ) 

  }

  ## return
  invisible( eval(parse(text=paste('new(mtxType, ',mtxType,'=mtx, ',
                          'A=anom.A, weight.type.A=weight.type.A, weights.A=weights.A, ',
                          'B=anom.B, weight.type.B=weight.type.B, weights.B=weights.B, ',
                          'lag=lag, fast=fast, complex=complex.A)',sep='') ) ) )
  
}

