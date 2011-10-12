## examples using the EOF package (soon to be a pacakge... ok maybe not that soon)

## some examples inspired/taken from the text you should buy/check out:
## Navarra. A., and V. Simoncini, 2010, A Guide to Empirical Orthogonal Functions for Climate Data Analysis,
##   DOI 10.1007 / 978 - 90 - 481 - 3702 - 2 5, c Springer Science+Business Media BV 

## examples:
## a time dominated EOF (non-fast)
## a space dominated EOF (fast)
## a allow.missing ex
## a lat weighting ex
## a user-supplied weights ex
## anomalys calculated in space ex
## examples of seasonal anomalies vs other.
ceof.prop=TRUE ## ceof of a propigating wave
## ceof of some geophysical phenom
## SVD
## lagged EOF/SVD



if (ceof.prop) { ## inspired by Navarra and Simoncini figure 5.6

  src("EOF.r")
  n=34; m=20     ## a synthetic propigating wave with 34 timesteps on a 20x20 spatial domain.
  xnum=2; ynum=0; omega=4
  u=seq(0,2*pi,length.out=m)
  v=seq(0,2*pi,length.out=n)
  y=u %o% rep(1,length(u))
  x=t(y)
  a=array(data=NA, dim=c(m,m,n))
  for (j in 1:n) a[,,j] = 10*sin(ynum*y+ xnum*x - omega*v[j]) + 10*rnorm(length(u));
  z=array(data=a, dim=c(m^2,n))
  propwave <- new("spaceTime", data=z, lon=as.vector(x), lat=as.vector(y),
                  POSIXct=as.POSIXct(1:34,origin='2000/1/1:00:00'),
                  data.name='wave', data.units='m') 
  corr.prop <- corrMtx(propwave, complex=TRUE) # disabled fast for CEOF, since results werent close enough.
  ceof.prop <- CEOF(corr.prop)
  gg.prop <- ggCEOF(ceof.prop, plot=1)  ## plot the first mode.
  
}
