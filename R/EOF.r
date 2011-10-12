## Some classes and methods use definitions from my pgwa code.
## any routines which do this are found in this file and not in pgwa.r.
######################################################################
## OOP EOF anlaysis (and other related things). structure is kinda a mess right now.
##
## classes: 
##          spaceTime
##          anomaly
##          corrMtx
##          covarMtx
##          EOF
##
## methods: input class(es) -> output class
##          subset: spaceTime -> spaceTime
##          anomaly: spaceTime -> anomaly
##          corrMtx: c("spaceTime", "spaceTimeNULLmissing") -> corrMtx
##          covarMtx: c("spaceTime", "spaceTimeNULLmissing") -> covarMtx
##          EOF: corrCovar -> EOF  (corrCovar is either corrMtx or covarMtx)

## james mccreight - mccreigh ^ gmail * com

######################################################################
## Note: all my slot accessor functions begin with "a." This means
## "access slot" and will be followed by the slot name and a function
## on an object, e.g. a.slot(obj). 

## spaceTime.r
#src("center.discrete.r") ## simple factorizing (as in cut()) may be used by gg routines.
#src("ggField.r")  ## plots spaceTime variables
#src("subset.r")   ## subsets spaceTime varibles in either space or time
#src("anomaly.r")  ## calculates anomalies on spaceTime variables
#src("corrCovarMtx.r")  ## calculates the correlation or covariance mtx of spaceTime variables.

## EOF.class.r
######################################################################

#src("ggEOF.r")
#src("unwrap.phase.r")
#src("hilbert.r")
#src("CEOF.r")
#src("analyticFilter.r")
#src("ggCEOF.r")
#src("spectralEOF.r")
#src("get.spectralEOF.by.mode.r")
#src("ggSpectralEOF.r")
#src("sequential_spectral_forecast.r")
#src("corrMap.r")
#src("ggCorrMap.r")
#src("ccfMap.r")
#src("ggCCFMap.r")
#src("world.map.r")
#rm(set.accessors)

## todo
## add lag capability to corr/covar
## recreate the example from clim.pact, mostly for eigen uncertainty
## can fast be done with complex EOFs?
## rotated EOFs
## corr/covar: *test*: allow missing values in calculation of eigen uncertainty
## need to preserve the total variance in each data set when geographical weighting
## what about a null corr/covarMtx: what if EOF analysis is desired on the original field?

######################################################################
## method: EOF
## benefits from recognizing the fast corrMtx/covarMtx by default
## in those methods / classes.
## The default nuber of modes is that which gives at least 95%
## of the total variance.

if (!isGeneric("EOF")) {  ## creates a generic function
  fun <- if (is.function("EOF")) EOF else function( mtx, ...) standardGeneric("EOF")
  setGeneric("EOF", fun)
}

setMethod("EOF", c("corrCovar"),
          function( mtx, nmodes=which(cumsum(variance.frac)>95)[1] ) {
           
            if (!is.null(a.B(mtx)))
              warning('EOF failure: analysis matrix should not be a cross-matrix.', immediate.=T)
            
            ## mtx = U D t(V)
            svd <- if (class(mtx)[[1]] == 'corrMtx') svd(a.corrMtx(mtx)) else svd(a.covarMtx(mtx))
            
            nspace=length(a.lon(a.anomaly(a.A(mtx))))
            ntime=length(a.POSIXct(a.anomaly(a.A(mtx))))

            eigen <- svd$d
            variance.frac <- eigen / sum(eigen) *100
            
            ## Uncertainty in the eigen values: adapted from clim.pact, which sites North.
            ## I have yet to verify this independently.
            ## Estimate the effective degrees of freedom at each point over its timeseries.

            ## I assume this only works when anomalies are calculated in time, which is consistent with
            ## the temporal lag correlation which is the basis of this method. All bets off for spatial
            ## correlations.
            if (a.anom.dim(a.A(mtx))!='time') {
              warning("EOF: Uncertainty in eigenvalues currently calculated on time anomalies only, set to zero here")
              eigen.uncert <- eigen*0
            } else {
              anom.dim.num <- if (a.anom.dim(a.A(mtx))=='space') 2 else 1
              l1corr <- function( ts ) acf( ts, lag.max=1, plot=FALSE )$acf[2,1,1]
              l1max <- max( apply( a.data(a.anomaly(a.A(mtx))), anom.dim.num, 'l1corr'), na.rm=TRUE)           
              n.eff <- round( ntime * (1 - l1max)/(1 + l1max) )
              eigen.uncert <- eigen * sqrt(2/n.eff)
            }           

            ## extract EOF and timeseries
            if (!a.fast(mtx)) {

              eof <- svd$u  ## spaceXmodes
              timeseries <- t(a.data(a.anomaly(a.A(mtx)))) %*% eof  ##  timeXspace %*% spaceXmodes 

            } else {

              eof <- array(NA, dim=c(nspace, ntime))  ## ntime could be replaced with nmode... or not.
              for (mm in 1:ntime) {
                tmp <- a.data(a.anomaly(a.A(mtx))) %*% svd$u[,mm]  ## spaceXtime %*% timeX1 = spaceX1
                eof[,mm] <- tmp / sqrt(sum(tmp^2))                 ## spaceX1
              }
              timeseries <- t(a.data(a.anomaly(a.A(mtx)))) %*% eof  ##  timeXspace %*% spaceXmodes 

            }

            ## truncate to nmodes
            ## for now this is the first nmodes to get to at least 95% of total variance.
            nmodes <- min(nmodes, dim(eof)[2])  ## in case nmodes is set to something else.
            nmodes <- max(nmodes,2)  ## needs to be at least 2
            eof <- eof[,1:nmodes]
            timeseries <- timeseries[,1:nmodes]
            eigen <- eigen[1:nmodes]
            eigen.uncert <- eigen.uncert[1:nmodes]
            variance.frac <- variance.frac[1:nmodes]
            
            invisible( new("EOF",
                           EOF=eof,
                           modes=1:nmodes,
                           lon=a.lon(a.anomaly(a.A(mtx))),
                           lat=a.lat(a.anomaly(a.A(mtx))),
                           timeseries=timeseries,
                           POSIXct=a.POSIXct(a.anomaly(a.A(mtx))),
                           eigen=eigen,
                           eigen.uncert=eigen.uncert, 
                           variance.frac=variance.frac,
                           variable=a.data.name(a.anomaly(a.A(mtx))),
                           corr.covar= if (a.standardize(a.A(mtx))) 'corr' else 'covar',
                           fast= a.fast(mtx),
                           mtx=mtx ) )
          }
          )
  


    
