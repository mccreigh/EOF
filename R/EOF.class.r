######################################################################
## class: EOF
## Because the B field is optional have an optional slot (B used in SVD / EEOF analysis).
setClassUnion("corrCovar", c("corrMtx", "covarMtx"))
## corrMtx
setClass("EOF", representation(EOF="array", modes="vector", lon="vector", lat="vector",
                               timeseries="array", POSIXct="POSIXct",
                               eigen="vector",
                               eigen.uncert="vector",
                               variance.frac= "numeric",
                               variable="character",
                               corr.covar="character",
                               fast= "logical",
                               mtx="corrCovar"
                               ) )
set.accessors("EOF")
