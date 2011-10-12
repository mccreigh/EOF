
world.map <- function( add=TRUE, lon=NULL, lat=NULL, buffer=0,
                      alpha=.75, color='black', size=.5, positive.lon=FALSE, ... ) {
    require(ggplot2)
    if (add) {
      require(maps)
      world <- map_data("world")

      if (!is.null(lon)) if (any(lon >180)) positive.lon=TRUE      
      if (positive.lon) world$long <- world$long %% 360
      
      if (!is.null(lon)) world <- world[ world$long <= max(lon)+buffer &
                                        world$long >= min(lon)-buffer, ]
      if (!is.null(lat)) world <- world[ world$lat <= max(lat)+buffer &
                                        world$lat >= min(lat)-buffer, ]
      invisible(geom_path( data=world, aes(x=long,y=lat,group=group),
                          alpha=alpha, color=color, size=size))
    } else invisible(geom_blank())
    
}
