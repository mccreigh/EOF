require(waveprac)
######################################################################
## class: spaceTime
setClass("spaceTime",
         representation(data = "array",
                        data.name = "character",
                        data.units = "character",
                        lon = "vector", 
                        lat = "vector", 
                        POSIXct = "POSIXct")
         )

## spaceTime accessor functions
set.accessors( 'spaceTime' )
